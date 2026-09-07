from .SCPIinst import SCPIinst
import numpy as np
import time

class Keysight33622A(SCPIinst):
    errcmds = {}
    errcmds['error'] = 'SYST:ERR?'

    commandlist = {}
    commandlist['core'] = {}

    core = {}

    core['phase_unit'] = 'UNIT:ANGLe'  # could be DEGree, RADian, SECond, DEFault

    core['Ch1_waveform'] = 'SOURce1:FUNCtion'  # could be SINusoid, DC, or others
    core['Ch1_frequency'] = 'SOURce1:FREQuency'
    core['Ch1_voltage'] = 'SOURce1:VOLTage'  # by default the unit is Vpp
    core['Ch1_phase'] = 'SOURce1:PHASe'      # by default the unit is DEGree
    core['Ch1_offset'] = 'SOURce1:VOLTage:OFFSet'  # by default the unit is V
    core['Ch1_output'] = 'OUTPut1'           # takes ON/1, OFF/0

    core['Ch2_waveform'] = 'SOURce2:FUNCtion'
    core['Ch2_frequency'] = 'SOURce2:FREQuency'
    core['Ch2_voltage'] = 'SOURce2:VOLTage'
    core['Ch2_phase'] = 'SOURce2:PHASe'
    core['Ch2_offset'] = 'SOURce2:VOLTage:OFFSet'
    core['Ch2_output'] = 'OUTPut2'

    commandlist['core'] = core

    def __init__(self, address, reset=False):
        """
        Initialize connection with instrument, reset it if specified, and clear all errors
        """
        self.instrument_type = 'Keysight33622A'
        super().__init__(address, self.commandlist, self.errcmds, reset)
        self.Ch1_output = 0
        self.Ch2_output = 0

        # for live ARB preview figures
        object.__setattr__(self, "_arb_preview", {})

        if reset:
            self.reset()

    def reset(self):
        self.inst.write('*RST; *CLS')
        self.inst.write('OUTP1 OFF')
        self.inst.write('OUTP2 OFF')

    def Ch1_dc_voltage_ramp(self, newV, step_size=0.005, step_time=0.001):
        self.Ch1_waveform = 'DC'
        deltaV = newV - self.Ch1_offset
        numSteps = max(2, int(np.abs(np.ceil(deltaV / step_size))))
        vsteps = np.linspace(self.Ch1_offset, newV, numSteps)
        for vstep in vsteps:
            self.Ch1_offset = np.round(vstep, 4)
            time.sleep(step_time)

    def Ch2_dc_voltage_ramp(self, newV, step_size=0.005, step_time=0.001):
        self.Ch2_waveform = 'DC'
        deltaV = newV - self.Ch2_offset
        numSteps = max(2, int(np.abs(np.ceil(deltaV / step_size))))
        vsteps = np.linspace(self.Ch2_offset, newV, numSteps)
        for vstep in vsteps:
            self.Ch2_offset = np.round(vstep, 4)
            time.sleep(step_time)

    def Ch1_sin_gen(self, amplitude, frequency, phase=0, offset=0, step_size=0.005, step_time=0.001):
        self.Ch1_waveform = 'SIN'
        self.Ch1_frequency = frequency
        self.Ch1_voltage = str(amplitude) + ' Vpp'
        self.Ch1_phase = phase

        deltaV = offset - self.Ch1_offset
        numSteps = max(2, int(np.abs(np.ceil(deltaV / step_size))))
        vsteps = np.linspace(self.Ch1_offset, offset, numSteps)
        for vstep in vsteps:
            self.Ch1_offset = np.round(vstep, 4)
            time.sleep(step_time)

    def Ch2_sin_gen(self, amplitude, frequency, phase=0, offset=0, step_size=0.005, step_time=0.001):
        self.Ch2_waveform = 'SIN'
        self.Ch2_frequency = frequency
        self.Ch2_voltage = str(amplitude) + ' Vpp'
        self.Ch2_phase = phase

        deltaV = offset - self.Ch2_offset
        numSteps = max(2, int(np.abs(np.ceil(deltaV / step_size))))
        vsteps = np.linspace(self.Ch2_offset, offset, numSteps)
        for vstep in vsteps:
            self.Ch2_offset = np.round(vstep, 4)
            time.sleep(step_time)

    def phase_sync(self):
        # source1 or source2 mean nothing here, simply give two channels a common reference point
        self.inst.write('SOURce1:PHASe:SYNChronize')

    # -------------------------------------------------------------------------
    # ARB helpers
    # -------------------------------------------------------------------------
    @staticmethod
    def _deg_to_rad(x):
        return np.deg2rad(x)

    @staticmethod
    def _pick_arb_length(Fs, f1, f2, n_min=2000, n_max=65536, tol_cycles=1e-9):
        """
        Pick an ARB length N such that:
            (f1 * N / Fs) and (f2 * N / Fs) are (close to) integers.
        This ensures seamless repetition when looping the record.
        """
        bestN = None
        bestErr = float('inf')

        if f1 <= 0 or f2 <= 0 or Fs <= 0:
            raise ValueError("Fs, f1, f2 must be positive.")

        for N in range(int(n_min), int(n_max) + 1):
            c1 = f1 * N / Fs
            c2 = f2 * N / Fs
            e1 = abs(c1 - round(c1))
            e2 = abs(c2 - round(c2))
            err = max(e1, e2)

            if err < bestErr:
                bestErr = err
                bestN = N
                if bestErr <= tol_cycles:
                    break

        return bestN, bestErr

    def _arb_live_preview_refresh(
        self,
        ch: int,
        t: np.ndarray,
        y_norm: np.ndarray,
        Fs: float,
        f1: float,
        f2: float,
        name: str,
        max_time_points: int = 2000,
        max_fft_mhz: float = 200.0,
    ):
        """
        Live-refresh a single preview figure per channel.
        Robust to SCPIinst.__getattr__ interception (doesn't rely on hasattr()).
        """
        import matplotlib.pyplot as plt

        # ---- robust attribute init: bypass SCPIinst.__getattr__ ----
        try:
            arb_preview = object.__getattribute__(self, "_arb_preview")
        except AttributeError:
            object.__setattr__(self, "_arb_preview", {})
            arb_preview = object.__getattribute__(self, "_arb_preview")

        key = f"ch{ch}"
        N = len(y_norm)

        # Prepare plotted time subset
        nplot = min(N, int(max_time_points))
        x_time = t[:nplot] * 1e9
        y_time = y_norm[:nplot]

        # FFT for preview
        w = np.hanning(N)
        Y = np.fft.rfft(y_norm * w)
        f = np.fft.rfftfreq(N, d=1.0 / Fs)
        mag = 20.0 * np.log10(np.maximum(np.abs(Y), 1e-15))

        # Limit FFT display range
        f_mhz = f / 1e6
        if max_fft_mhz is not None:
            m = f_mhz <= float(max_fft_mhz)
            f_mhz = f_mhz[m]
            mag = mag[m]

        # ---------- Create figure first time ----------
        if key not in arb_preview:
            plt.ion()
            fig = plt.figure(num=100 + ch, figsize=(10, 7))
            fig.clf()

            try:
                fig.canvas.manager.set_window_title(f"ARB Preview CH{ch}")
            except Exception:
                pass

            ax1 = fig.add_subplot(2, 1, 1)
            ax2 = fig.add_subplot(2, 1, 2)

            (line_time,) = ax1.plot(x_time, y_time, marker=".", linestyle="-")
            ax1.set_xlabel("Time (ns)")
            ax1.set_ylabel("y_norm (arb)")
            ax1.grid(True)

            (line_fft,) = ax2.plot(f_mhz, mag, linestyle="-")
            ax2.set_xlabel("Frequency (MHz)")
            ax2.set_ylabel("Magnitude (dB, arb)")
            ax2.grid(True)

            arb_preview[key] = {
                "fig": fig,
                "ax1": ax1,
                "ax2": ax2,
                "line_time": line_time,
                "line_fft": line_fft,
            }

        # ---------- Update existing ----------
        d = arb_preview[key]

        d["line_time"].set_data(x_time, y_time)
        d["ax1"].relim()
        d["ax1"].autoscale_view(scalex=True, scaley=True)

        d["line_fft"].set_data(f_mhz, mag)
        d["ax2"].relim()
        d["ax2"].autoscale_view(scalex=True, scaley=True)

        title = (
            f"CH{ch} '{name}' | N={N}, Fs={Fs/1e9:.3f} GS/s | "
            f"f1={f1/1e6:.3f} MHz, f2={f2/1e6:.3f} MHz | rep={Fs/N/1e6:.3f} MHz"
        )
        d["fig"].suptitle(title)

        d["fig"].canvas.draw()
        d["fig"].canvas.flush_events()


    def _upload_arb_dual_sine(
        self,
        ch: int,
        offset: float,
        amplitude_1: float,
        frequency_1: float,
        phase_1: float,
        amplitude_2: float,
        frequency_2: float,
        phase_2: float,
        Fs: float = 1e9,
        name: str = None,
        n_min: int = 2000,
        n_max: int = 65536,
        tol_cycles: float = 1e-9,
        plot: bool = False,
        plot_fft_mhz: float = 200.0,
        arb_filter: str = "NORM",
    ):
        """
        Build + upload a dual-sine ARB to channel ch (1 or 2) using binary DAC download.
        Sets sample rate Fs, channel Vpp=(A1+A2), offset, selects ARB function.

        plot=True: live-refreshes a preview figure (Spyder-friendly).
        """
        if ch not in (1, 2):
            raise ValueError("ch must be 1 or 2")

        # IMPORTANT: arb name length on 33600-series is limited; keep <= 12 chars
        if name is None:
            name = f"DUALSINE{ch}"  # <=12 chars

        A1 = float(amplitude_1)
        A2 = float(amplitude_2)
        f1 = float(frequency_1)
        f2 = float(frequency_2)
        ph1 = self._deg_to_rad(float(phase_1))
        ph2 = self._deg_to_rad(float(phase_2))

        if A1 < 0 or A2 < 0:
            raise ValueError("amplitude_1 and amplitude_2 must be >= 0")
        if f1 <= 0 or f2 <= 0:
            raise ValueError("frequency_1 and frequency_2 must be > 0")

        src = f"SOURce{ch}"

        # Choose record length N so the record repeats seamlessly
        N, err = self._pick_arb_length(Fs, f1, f2, n_min=n_min, n_max=n_max, tol_cycles=tol_cycles)

        # Time axis
        t = np.arange(N, dtype=np.float64) / Fs

        # Channel Vpp scaling
        Vpp_total = A1 + A2
        if Vpp_total <= 0:
            Vpp_total = 0.001
            y_norm = np.zeros(N, dtype=np.float64)
        else:
            y_norm = (A1 / Vpp_total) * np.sin(2 * np.pi * f1 * t + ph1) \
                   + (A2 / Vpp_total) * np.sin(2 * np.pi * f2 * t + ph2)

            # protect against clipping
            maxabs = float(np.max(np.abs(y_norm)))
            if maxabs > 1.0:
                y_norm = y_norm / maxabs

        # Live preview plot (refreshes the same figure each call)
        if plot:
            self._arb_live_preview_refresh(
                ch=ch, t=t, y_norm=y_norm, Fs=Fs, f1=f1, f2=f2, name=name,
                max_time_points=2000,
                max_fft_mhz=plot_fft_mhz,
            )

        # Set ARB reconstruction filter explicitly (recommended)
        # Use NORM for smoother sine-like output; STEP for better step response.
        self.inst.write(f"{src}:FUNCtion:ARBitrary:FILTer {arb_filter}")

        # Convert normalized [-1,1] -> int16 DAC codes [-32767,32767]
        dac = np.clip(y_norm, -1.0, 1.0)
        dac = np.round(dac * 32767).astype(np.int16)

        # IEEE definite-length binary block: #<n><len><data>
        payload = dac.tobytes(order="C")
        len_str = str(len(payload))
        header = f"#{len(len_str)}{len_str}".encode("ascii")

        # Clear volatile for this source, then send binary block
        self.inst.write(f"{src}:DATA:VOLatile:CLEar")

        # Binary DAC download (no quotes around name)
        cmd = f"{src}:DATA:ARBitrary:DAC {name},".encode("ascii")
        self.inst.write_raw(cmd + header + payload)

        # Select ARB and sample rate
        self.inst.write(f'{src}:FUNCtion:ARB "{name}"')
        self.inst.write(f"{src}:FUNCtion ARB")
        self.inst.write(f"{src}:FUNCtion:ARB:SRATe {Fs}")

        # Set amplitude + offset at the channel level
        self.inst.write(f"{src}:VOLTage {Vpp_total}")
        self.inst.write(f"{src}:VOLTage:OFFSet {float(offset)}")

        return {
            "channel": ch,
            "name": name,
            "Fs": Fs,
            "N": N,
            "repeat_freq": Fs / N,
            "cycle_err_max": err,
            "Vpp_set": Vpp_total,
            "offset_set": float(offset),
        }

    # -------------------------------------------------------------------------
    # Public ARB dual-sine methods (plot arguments supported)
    # -------------------------------------------------------------------------
    def Ch1_arb_dual_sine(
        self, offset, amplitude_1, frequency_1, phase_1, amplitude_2, frequency_2, phase_2,
        plot=False, plot_fft_mhz=200.0, arb_filter="NORM", n_min=2000
    ):
        return self._upload_arb_dual_sine(
            ch=1,
            offset=offset,
            amplitude_1=amplitude_1, frequency_1=frequency_1, phase_1=phase_1,
            amplitude_2=amplitude_2, frequency_2=frequency_2, phase_2=phase_2,
            Fs=1e9,
            name="DUALSINE1",     # <=12 chars
            plot=plot,
            plot_fft_mhz=plot_fft_mhz,
            arb_filter=arb_filter,
            n_min=n_min,
        )

    def Ch2_arb_dual_sine(
        self, offset, amplitude_1, frequency_1, phase_1, amplitude_2, frequency_2, phase_2,
        plot=False, plot_fft_mhz=200.0, arb_filter="NORM", n_min=2000
    ):
        return self._upload_arb_dual_sine(
            ch=2,
            offset=offset,
            amplitude_1=amplitude_1, frequency_1=frequency_1, phase_1=phase_1,
            amplitude_2=amplitude_2, frequency_2=frequency_2, phase_2=phase_2,
            Fs=1e9,
            name="DUALSINE2",     # <=12 chars
            plot=plot,
            plot_fft_mhz=plot_fft_mhz,
            arb_filter=arb_filter,
            n_min=n_min,
        )

    def arb_sync(self):
        """
        Sync both channels' ARB playback to start at the first sample together.
        """
        self.inst.write("FUNCtion:ARB:SYNChronize")
