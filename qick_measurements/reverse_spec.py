# -*- coding: utf-8 -*-

from qick.averager_program import AveragerProgram
import numpy as np
import os
import time
import matplotlib.pyplot as plt
import userfuncs
from utility.plotting_tools import simplescan_plot
from utility.measurement_helpers import estimate_time

try:
    from scipy.optimize import curve_fit
except Exception:
    curve_fit = None


def lorentzian_peak(x, A, f0, gamma_fwhm, offset):
    return offset + A * (0.5 * gamma_fwhm) ** 2 / ((x - f0) ** 2 + (0.5 * gamma_fwhm) ** 2)

def lorentzian_dip(x, A, f0, gamma_fwhm, offset):
    return offset - A * (0.5 * gamma_fwhm) ** 2 / ((x - f0) ** 2 + (0.5 * gamma_fwhm) ** 2)


class ReverseCavitySpec(AveragerProgram):
    """
    Quasi_CW timing: qubit pulse then cavity readout pulse.
    Cavity frequency is set from Python sweep (re-instantiating program per point like pulsed_trans).
    """
    def initialize(self):
        cfg = self.cfg
        cav_ch = cfg["cav_channel"]
        qub_ch = cfg["qub_channel"]

        self.declare_gen(ch=cav_ch, nqz=cfg["nqz_c"])
        self.declare_gen(ch=qub_ch, nqz=cfg["nqz_q"])

        # readout
        readout = self.us2cycles(cfg["readout_length"], ro_ch=cfg["ro_channels"][0])
        for ch in cfg["ro_channels"]:
            self.declare_readout(ch=ch, length=readout, freq=cfg["cav_freq"], gen_ch=cav_ch)

        # cavity pulse
        freq_c  = self.freq2reg(cfg["cav_freq"], gen_ch=cav_ch, ro_ch=cfg["ro_channels"][0])
        phase_c = self.deg2reg(cfg["cav_phase"], gen_ch=cav_ch)
        gain_c  = int(cfg["meas_gain"])

        self.default_pulse_registers(ch=cav_ch, freq=freq_c, phase=phase_c, gain=gain_c)
        self.set_pulse_registers(ch=cav_ch, style="const",
                                 length=self.us2cycles(cfg["meas_window"], gen_ch=cav_ch))

        # qubit pulse (flat_top + gaussian edges)
        freq_q  = self.freq2reg(cfg["qub_freq"], gen_ch=qub_ch)
        phase_q = self.deg2reg(cfg["qub_phase"], gen_ch=qub_ch)
        gain_q  = int(cfg["qub_gain"])

        self.default_pulse_registers(ch=qub_ch, freq=freq_q, phase=phase_q, gain=gain_q)

        sigma = self.us2cycles(cfg["qub_sigma"], gen_ch=qub_ch)
        num_sigma = int(cfg["num_sigma"])
        self.add_gauss(ch=qub_ch, name="ex", sigma=sigma, length=int(sigma * num_sigma))
        self.set_pulse_registers(ch=qub_ch, style="flat_top", waveform="ex",
                                 length=self.us2cycles(cfg["quasi_CW_len"], gen_ch=qub_ch))

        self.synci(200)

    def body(self):
        cfg = self.cfg
        cav_ch = cfg["cav_channel"]
        qub_ch = cfg["qub_channel"]

        sigma = self.us2cycles(cfg["qub_sigma"], gen_ch=qub_ch)
        num_sigma = int(cfg["num_sigma"])
        pulse_len = self.us2cycles(cfg["quasi_CW_len"], gen_ch=qub_ch) + int(num_sigma * sigma)

        offset = self.us2cycles(cfg["adc_trig_offset"], gen_ch=cav_ch)
        meas_time = self.us2cycles(cfg["meas_time"], gen_ch=cav_ch)

        ex_time = meas_time - self.us2cycles(cfg["qub_delay"], gen_ch=qub_ch) - pulse_len

        self.trigger(adcs=self.ro_chs, pins=[0], adc_trig_offset=offset)
        self.pulse(ch=qub_ch, t=ex_time)
        self.pulse(ch=cav_ch, t=meas_time)
        self.wait_all()
        self.sync_all(self.us2cycles(cfg["relax_delay"]))


def get_reverse_cav_spec_settings():
    s = {}
    s["scanname"] = "reverse_cav_spec"
    s["meas_type"] = "ReverseCavSpec"

    # cavity frequency sweep (Hz)
    s["freq_start"]  = 7e9
    s["freq_step"]   = 1e6
    s["freq_points"] = 401

    # qubit gain sweep (DAC a.u.)  <-- THIS is the outer axis now
    s["qub_gain_start"]  = 1000
    s["qub_gain_step"]   = 200
    s["qub_gain_points"] = 1   # set >1 to sweep qubit power

    # fixed knobs for qubit tone
    s["qub_freq"]      = 4.5e9   # Hz
    s["quasi_CW_len"]  = 10.0    # us

    # fixed cavity gain
    s["meas_gain"] = 500         # DAC a.u. (fixed cavity/readout gain)

    # options
    s["fit_last_trace"] = True

    # card
    s["reps"] = 1
    s["soft_avgs"] = 2000
    return s


def reverse_cavity_spec(soc, soccfg, instruments, settings):
    exp_globals  = settings["exp_globals"]
    exp_settings = settings["exp_settings"]
    m_pulse = exp_globals["measurement_pulse"]
    q_pulse = exp_globals["qubit_pulse"]

    lo_freq = exp_globals.get("LO_freq", 0.0) * exp_globals.get("LO", 0)

    stamp    = userfuncs.timestamp()
    saveDir  = userfuncs.saveDir(settings)
    filename = exp_settings["scanname"] + "_" + stamp

    # sweep axes
    fpts = exp_settings["freq_start"] + exp_settings["freq_step"] * np.arange(exp_settings["freq_points"])
    qgpts = exp_settings["qub_gain_start"] + exp_settings["qub_gain_step"] * np.arange(exp_settings["qub_gain_points"])

    mags   = np.zeros((len(qgpts), len(fpts)))
    phases = np.zeros((len(qgpts), len(fpts)))
    Is     = np.zeros((len(qgpts), len(fpts)))
    Qs     = np.zeros((len(qgpts), len(fpts)))

    tstart = time.time()

    for ig, qub_gain in enumerate(qgpts):
        print(f"Current qubit gain: {qub_gain} (row {ig+1}/{len(qgpts)})")

        for jf, cav_f_hz in enumerate(fpts):
            cav_board_MHz = (cav_f_hz - lo_freq) / 1e6

            config = {
                "cav_channel": exp_globals["cav_channel"],
                "qub_channel": exp_globals["qub_channel"],
                "ro_channels": exp_globals["ro_channels"],

                "nqz_c": exp_globals.get("nqz_c", 1),
                "nqz_q": exp_globals.get("nqz_q", 2),

                # cavity/readout pulse (frequency swept, gain fixed)
                "cav_freq":   cav_board_MHz,
                "cav_phase":  m_pulse["cav_phase"],
                "meas_gain":  int(exp_settings["meas_gain"]),
                "meas_window": m_pulse["meas_window"],
                "meas_time":   m_pulse["meas_pos"],

                # qubit pulse (gain swept, freq fixed)
                "qub_freq":   exp_settings["qub_freq"] / 1e6,
                "qub_phase":  q_pulse["qub_phase"],
                "qub_gain":   int(qub_gain),
                "qub_sigma":  q_pulse["sigma"],
                "qub_delay":  q_pulse["delay"],
                "num_sigma":  q_pulse["num_sigma"],
                "quasi_CW_len": float(exp_settings["quasi_CW_len"]),

                # readout timing
                "readout_length":  m_pulse["init_buffer"] + m_pulse["meas_window"]
                                  if "init_buffer" in m_pulse else m_pulse["meas_window"],
                "adc_trig_offset": m_pulse["emp_delay"] + m_pulse["meas_pos"],

                "relax_delay": exp_globals["relax_delay"],
                "reps": int(exp_settings["reps"]),
                "soft_avgs": int(exp_settings["soft_avgs"]),
            }

            prog = ReverseCavitySpec(soccfg, config)
            I_buf, Q_buf = prog.acquire(soc, load_pulses=True, progress=False)

            I = float(I_buf[0][0])
            Q = float(Q_buf[0][0])

            mags[ig, jf]   = np.sqrt(I**2 + Q**2)
            phases[ig, jf] = np.arctan2(Q, I) * 180 / np.pi
            Is[ig, jf] = I
            Qs[ig, jf] = Q

        # estimate after first row (like pulsed_trans)
        if ig == 0 and len(qgpts) > 1:
            estimate_time(tstart, time.time(), len(qgpts))

        # live plot (yaxis is qubit gain now)
        full_data = {
            "xaxis":  fpts / 1e9,
            "mags":   mags[:ig+1],
            "phases": phases[:ig+1],
            "Is":     Is[:ig+1],
            "Qs":     Qs[:ig+1],
        }
        plot_data = {"xaxis": fpts / 1e9, "mags": mags[:ig+1], "phases": phases[:ig+1]}
        single_data = {"xaxis": fpts / 1e9, "mag": mags[ig], "phase": phases[ig]}

        yaxis  = qgpts[:ig+1]
        labels = ["Cavity Freq (GHz)", "Qubit Gain (DAC a.u.)"]

        simplescan_plot(plot_data, single_data, yaxis, filename, labels, identifier="", fig_num=1, IQdata=False)
        plt.savefig(os.path.join(saveDir, filename + "_fullColorPlot.png"), dpi=150)

        userfuncs.SaveFull(
            saveDir, filename,
            ["qgpts", "fpts", "mags", "phases", "Is", "Qs", "full_data", "single_data"],
            locals(), expsettings=settings, instruments={}
        )

    # Optional Lorentzian fit on last qubit-gain trace
    if exp_settings.get("fit_last_trace", True) and curve_fit is not None and len(fpts) >= 4:
        x = fpts / 1e9
        y = mags[-1].copy()

        hanger = bool(exp_globals.get("hanger", False))  # True => dip, False => peak
        model = lorentzian_dip if hanger else lorentzian_peak

        offset0 = float(np.median(y))
        idx0 = int(np.argmin(y) if hanger else np.argmax(y))
        f0_0 = float(x[idx0])

        step = float(np.mean(np.diff(x))) if len(x) > 1 else 0.001
        gamma0 = max(5 * step, 0.001)  # gamma is FWHM in GHz for this model
        A0 = float((offset0 - y[idx0]) if hanger else (y[idx0] - offset0))
        A0 = max(A0, 1e-9)

        p0 = [A0, f0_0, gamma0, offset0]
        bounds = ([0, x.min(), 0, -np.inf], [np.inf, x.max(), np.inf, np.inf])

        try:
            popt, _ = curve_fit(model, x, y, p0=p0, bounds=bounds, maxfev=20000)
        except Exception:
            popt = p0

        A, f0, gamma_fwhm, offset = [float(v) for v in popt]
        FWHM_MHz = gamma_fwhm * 1e3

        fig2, ax2 = plt.subplots(figsize=(7, 4.5), num=2, clear=True)
        fig2.suptitle(filename, fontsize=11, y=0.98)
        ax2.plot(x, y, "o", ms=4, label="Data (last qub gain)")
        xf = np.linspace(x.min(), x.max(), 1000)
        ax2.plot(xf, model(xf, *popt), "-", lw=2, label="Lorentzian fit")
        ax2.set_xlabel("Frequency (GHz)")
        ax2.set_ylabel("|S21| (a.u.)")
        ax2.legend()
        ax2.grid(alpha=0.3)
        kind = "dip" if hanger else "peak"
        ax2.set_title(f"Lorentzian {kind}: f0 = {f0:.7f} GHz, FWHM = {FWHM_MHz:.2f} MHz", fontsize=10)
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        plt.savefig(os.path.join(saveDir, filename + "_lasttrace_lorentz_fit.png"), dpi=150)

    full_data = {
        "xaxis":  fpts / 1e9,
        "yaxis":  qgpts,
        "mags":   mags,
        "phases": phases,
        "Is":     Is,
        "Qs":     Qs,
    }
    return full_data
