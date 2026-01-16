from .SCPIinst import SCPIinst
import numpy as np
import time

class Keysight33622A(SCPIinst):
    errcmds = {}
    errcmds['error'] = 'SYST:ERR?'

    commandlist = {}
    commandlist['core'] = {}

    core = {}

    core['phase_unit'] = 'UNIT:ANGLe' #could be DEGree, RADian, SECond, DEFault

    core['Ch1_waveform'] = 'SOURce1:FUNCtion' #could be SINusoid, DC, or others
    core['Ch1_frequency'] = 'SOURce1:FREQuency' 
    core['Ch1_voltage'] = 'SOURce1:VOLTage' #by default the unit is Vpp
    core['Ch1_phase'] = 'SOURce1:PHASe' #by default the unit is DEGree
    core['Ch1_offset'] = 'SOURce1:VOLTage:OFFSet'  #by default the unit is V
    core['Ch1_output'] = 'OUTPut1' #takes ON/1, OFF/0

    core['Ch2_waveform'] = 'SOURce2:FUNCtion'
    core['Ch2_frequency'] = 'SOURce2:FREQuency'
    core['Ch2_voltage'] = 'SOURce2:VOLTage'
    core['Ch2_phase'] = 'SOURce2:PHASe'
    core['Ch2_offset'] = 'SOURce2:VOLTage:OFFSet' 
    core['Ch2_output'] = 'OUTPut2'

    commandlist['core'] = core

    def __init__(self, address, reset = False):
        '''
        Initialize connection with instrument, reset it if specified, and clear all errors
        '''
        self.instrument_type = 'Keysight33622A'

        super().__init__(address, self.commandlist, self.errcmds, reset)
        self.Ch1_output = 0
        self.Ch2_output = 0

        if reset:
            self.reset()

    def reset(self):
        self.inst.write('*RST; *CLS')
        self.inst.write('OUTP1 OFF')
        self.inst.write('OUTP2 OFF')

    def Ch1_dc_voltage_ramp(self, newV, step_size = 0.005, step_time = 0.001):
        self.Ch1_waveform = 'DC'
        #self.Ch1_output = 1
        deltaV = newV - self.Ch1_offset
        numSteps = max(2, int(np.abs(np.ceil(deltaV/step_size))))
        vsteps = np.linspace(self.Ch1_offset, newV, numSteps)
        for vstep in vsteps:
            self.Ch1_offset = np.round(vstep,4)
            time.sleep(step_time)

    def Ch2_dc_voltage_ramp(self, newV, step_size = 0.005, step_time = 0.001):
        self.Ch2_waveform = 'DC'
        #self.Ch2_output = 1
        deltaV = newV - self.Ch2_offset
        numSteps = max(2, int(np.abs(np.ceil(deltaV/step_size))))
        vsteps = np.linspace(self.Ch2_offset, newV, numSteps)
        for vstep in vsteps:
            self.Ch2_offset = np.round(vstep,4)
            time.sleep(step_time)

    def Ch1_sin_gen(self, amplitude, frequency, phase=0, offset=0, step_size = 0.005, step_time = 0.001):
        self.Ch1_waveform = 'SIN'
        self.Ch1_frequency = frequency
        self.Ch1_voltage = str(amplitude) + ' Vpp'
        self.Ch1_phase = phase
        #self.Ch1_output = 1

        deltaV = offset - self.Ch1_offset
        numSteps = max(2, int(np.abs(np.ceil(deltaV/step_size))))
        vsteps = np.linspace(self.Ch1_offset, offset, numSteps)
        for vstep in vsteps:
            self.Ch1_offset = np.round(vstep,4)
            time.sleep(step_time)

    def Ch2_sin_gen(self, amplitude, frequency, phase=0, offset=0, step_size = 0.005, step_time = 0.001):
        self.Ch2_waveform = 'SIN'
        self.Ch2_frequency = frequency
        self.Ch2_voltage = str(amplitude) + ' Vpp'
        self.Ch2_phase = phase
        #self.Ch2_output = 1

        deltaV = offset - self.Ch2_offset
        numSteps = max(2, int(np.abs(np.ceil(deltaV/step_size))))
        vsteps = np.linspace(self.Ch2_offset, offset, numSteps)
        for vstep in vsteps:
            self.Ch2_offset = np.round(vstep,4)
            time.sleep(step_time)

    def phase_sync(self):
        self.inst.write('SOURce1:PHASe:SYNChronize') #source1 or source2 mean nothing here, simply give two channels a common reference point


    @staticmethod
    def _deg_to_rad(x):
        return np.deg2rad(x)

    @staticmethod
    def _pick_arb_length(Fs, f1, f2, n_min=32, n_max=65536, tol_cycles=1e-9):
        """
        Pick an ARB length N such that:
            (f1 * N / Fs) and (f2 * N / Fs) are (close to) integers.
        This ensures seamless repetition when looping the record.

        Returns: (N_best, err_best)
        where err_best is max fractional-cycle error across both tones.
        """
        bestN = None
        bestErr = float('inf')

        # Quick sanity
        if f1 <= 0 or f2 <= 0 or Fs <= 0:
            raise ValueError("Fs, f1, f2 must be positive.")

        for N in range(n_min, n_max + 1):
            c1 = f1 * N / Fs
            c2 = f2 * N / Fs
            e1 = abs(c1 - round(c1))
            e2 = abs(c2 - round(c2))
            err = max(e1, e2)

            if err < bestErr:
                bestErr = err
                bestN = N
                if bestErr <= tol_cycles:
                    break  # good enough / exact within tolerance

        return bestN, bestErr

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
        n_max: int = 65536,
        tol_cycles: float = 1e-9,
    ):
        """
        Internal method: build and upload a dual-sine ARB to channel ch (1 or 2),
        set sample rate Fs, set Vpp to (A1 + A2), set offset, select ARB function.
        """
        if ch not in (1, 2):
            raise ValueError("ch must be 1 or 2")

        # Default waveform name (unique-ish per channel)
        if name is None:
            name = f"DUAL_SINE_CH{ch}"

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

        # Choose ARB length N so the record repeats seamlessly
        N, err = self._pick_arb_length(Fs, f1, f2, n_min=8, n_max=n_max, tol_cycles=tol_cycles)

        # Time axis
        t = np.arange(N, dtype=np.float64) / Fs

        # Choose channel Vpp scaling
        Vpp_total = A1 + A2
        if Vpp_total <= 0:
            # Edge case: both amplitudes 0 => just output DC offset
            # We'll upload zeros and set Vpp to something small but nonzero (instrument may reject 0 Vpp)
            Vpp_total = 0.001  # 1 mVpp as a safe tiny amplitude
            y_norm = np.zeros(N, dtype=np.float64)
        else:
            # Normalize so that setting channel Vpp to (A1+A2) gives correct Vpp components
            y_norm = (A1 / Vpp_total) * np.sin(2 * np.pi * f1 * t + ph1) \
                   + (A2 / Vpp_total) * np.sin(2 * np.pi * f2 * t + ph2)

            # Guard against tiny numerical overshoots
            maxabs = np.max(np.abs(y_norm))
            if maxabs > 1.0:
                y_norm = y_norm / maxabs  # avoid clipping in the ARB DAC
                # NOTE: This renormalization will slightly change effective A1/A2 ratio if clipping would occur.

        # Format points for SCPI ASCII list
        # Keep a reasonable precision; too many digits makes transfers slower.
        points_str = ",".join(f"{v:.8f}" for v in y_norm)

        src = f"SOURce{ch}"

        # ARB reconstruction filter: keep ON (STEP or NORM) if you want 1 GS/s
        self.inst.write(f"{src}:FUNCtion:ARBitrary:FILTer NORM")  # or STEP


        # Clear volatile ARB memory for that source (optional but helps avoid name collisions)
        # (Command is global DATA:VOLatile:CLEar; many setups accept per-source prefix too.)
        self.inst.write(f"{src}:DATA:VOLatile:CLEar")

        # Download data to volatile memory
        # Keysight supports: DATA:ARB <name>,<p1>,<p2>,...,<pN>
        self.inst.write(f"{src}:DATA:ARB {name},{points_str}")

        # Select ARB and sample rate
        self.inst.write(f'{src}:FUNCtion:ARB "{name}"')
        self.inst.write(f"{src}:FUNCtion ARB")
        self.inst.write(f"{src}:FUNCtion:ARB:SRATe {Fs}")

        # Set amplitude + offset at the channel level
        self.inst.write(f"{src}:VOLTage {Vpp_total}")
        self.inst.write(f"{src}:VOLTage:OFFSet {float(offset)}")

        # Optional: you can query/print err if you want visibility.
        # return useful info for logs
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

    # ---------------------------
    # New public methods
    # ---------------------------

    def Ch1_arb_dual_sine(self, offset, amplitude_1, frequency_1, phase_1, amplitude_2, frequency_2, phase_2):
        """
        Upload and play a dual-sine ARB on CH1 at 1 GS/s.
        Inputs:
          offset: V
          amplitude_1, amplitude_2: Vpp
          frequency_1, frequency_2: Hz
          phase_1, phase_2: degrees
        """
        return self._upload_arb_dual_sine(
            ch=1,
            offset=offset,
            amplitude_1=amplitude_1, frequency_1=frequency_1, phase_1=phase_1,
            amplitude_2=amplitude_2, frequency_2=frequency_2, phase_2=phase_2,
            Fs=1e9,
            name="DUAL_SINE_CH1",
        )

    def Ch2_arb_dual_sine(self, offset, amplitude_1, frequency_1, phase_1, amplitude_2, frequency_2, phase_2):
        """
        Upload and play a dual-sine ARB on CH2 at 1 GS/s.
        Inputs:
          offset: V
          amplitude_1, amplitude_2: Vpp
          frequency_1, frequency_2: Hz
          phase_1, phase_2: degrees
        """
        return self._upload_arb_dual_sine(
            ch=2,
            offset=offset,
            amplitude_1=amplitude_1, frequency_1=frequency_1, phase_1=phase_1,
            amplitude_2=amplitude_2, frequency_2=frequency_2, phase_2=phase_2,
            Fs=1e9,
            name="DUAL_SINE_CH2",
        )

    def arb_sync(self):
        """
        Sync both channels' ARB playback to start at the first sample together.
        This is the ARB equivalent of phase_sync().
        """
        # This is a global command on 33600-series: it stops/restarts both ARBs aligned.
        self.inst.write("FUNCtion:ARB:SYNChronize")