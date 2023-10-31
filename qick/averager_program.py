"""
Several helper classes for writing qubit experiments.
"""
from typing import List, Union
import numpy as np
from .qick_asm import QickProgram, obtain, QickRegister, QickRegisterManagerMixin

class AveragerProgram(QickProgram):
    """
    AveragerProgram class is an abstract base class for programs which do loops over experiments in hardware.
    It consists of a template program which takes care of the loop and acquire methods that talk to the processor to stream single shot data in real-time and then reshape and average it appropriately.

    :param soccfg: This can be either a QickSOc object (if the program is running on the QICK) or a QickCOnfig (if running remotely).
    :type soccfg: QickConfig
    :param cfg: Configuration dictionary
    :type cfg: dict
    """

    def __init__(self, soccfg, cfg):
        """
        Constructor for the AveragerProgram, calls make program at the end.
        For classes that inherit from this, if you want it to do something before the program is made and compiled:
        either do it before calling this __init__ or put it in the initialize method.
        """
        super().__init__(soccfg)
        self.cfg = cfg
        self.make_program()
        if "soft_avgs" in cfg:
            self.rounds = cfg['soft_avgs']
        if "rounds" in cfg:
            self.rounds = cfg['rounds']
        self.reps = cfg['reps']

    def initialize(self):
        """
        Abstract method for initializing the program and can include any instructions that are executed once at the beginning of the program.
        """
        pass

    def body(self):
        """
        Abstract method for the body of the program
        """
        pass

    def make_program(self):
        """
        A template program which repeats the instructions defined in the body() method the number of times specified in self.cfg["reps"].
        """
        p = self

        rjj = 14
        rcount = 15
        p.initialize()
        p.regwi(0, rcount, 0)
        p.regwi(0, rjj, self.cfg['reps']-1)
        p.label("LOOP_J")

        p.body()

        p.mathi(0, rcount, rcount, "+", 1)

        p.memwi(0, rcount, self.counter_addr)

        p.loopnz(0, rjj, 'LOOP_J')

        p.end()


    def acquire(self, soc, threshold=None, angle=None, readouts_per_experiment=1, save_experiments=None, load_pulses=True, start_src="internal", progress=False, debug=False):
        """
        This method optionally loads pulses on to the SoC, configures the ADC readouts, loads the machine code representation of the AveragerProgram onto the SoC, starts the program and streams the data into the Python, returning it as a set of numpy arrays.
        config requirements:
        "reps" = number of repetitions;

        :param soc: Qick object
        :type soc: Qick object
        :param threshold: threshold
        :type threshold: int
        :param angle: rotation angle
        :type angle: list
        :param readouts_per_experiment: readouts per experiment
        :type readouts_per_experiment: int
        :param save_experiments: saved readouts (by default, save all readouts)
        :type save_experiments: list
        :param load_pulses: If true, loads pulses into the tProc
        :type load_pulses: bool
        :param start_src: "internal" (tProc starts immediately) or "external" (each round waits for an external trigger)
        :type start_src: string
        :param progress: If true, displays progress bar
        :type progress: bool
        :param debug: If true, displays assembly code for tProc program
        :type debug: bool
        :returns:
            - expt_pts (:py:class:`list`) - list of experiment points
            - avg_di (:py:class:`list`) - list of lists of averaged accumulated I data for ADCs 0 and 1
            - avg_dq (:py:class:`list`) - list of lists of averaged accumulated Q data for ADCs 0 and 1
        """

        self.shot_angle = angle
        self.shot_threshold = threshold

        if save_experiments is None:
            save_experiments = range(readouts_per_experiment)

        d_buf, avg_d, shots = super().acquire(soc, reads_per_rep=readouts_per_experiment, load_pulses=load_pulses, start_src=start_src, progress=progress, debug=debug)

        # reformat the data into separate I and Q arrays
        # save results to class in case you want to look at it later or for analysis
        self.di_buf = d_buf[:,:,0]
        self.dq_buf = d_buf[:,:,1]

        if threshold is not None:
            self.shots = shots

        n_ro = len(self.ro_chs)
        avg_di = np.zeros((n_ro, len(save_experiments)))
        avg_dq = np.zeros((n_ro, len(save_experiments)))

        for nn, ii in enumerate(save_experiments):
            for i_ch, (ch, ro) in enumerate(self.ro_chs.items()):
                avg_di[i_ch][nn] = avg_d[i_ch, ii, 0]
                avg_dq[i_ch][nn] = avg_d[i_ch, ii, 1]

        return avg_di, avg_dq


    def acquire_decimated(self, soc, load_pulses=True, readouts_per_experiment=1, start_src="internal", progress=True, debug=False):
        """
        This method acquires the raw (downconverted and decimated) data sampled by the ADC. This method is slow and mostly useful for lining up pulses or doing loopback tests.

        config requirements:
        "reps" = number of tProc loop repetitions;
        "soft_avgs" = number of Python loop repetitions;

        The data is returned as a list of ndarrays (one ndarray per readout channel).
        There are two possible array formats.
        reps = 1:
        2D array with dimensions (2, length), indices (I/Q, sample)
        reps > 1:
        3D array with dimensions (reps, 2, length), indices (rep, I/Q, sample)
        readouts_per_experiment>1:
        3D array with dimensions (reps, expts, 2, length), indices (rep, expt, I/Q, sample)

        :param soc: Qick object
        :type soc: Qick object
        :param load_pulses: If true, loads pulses into the tProc
        :type load_pulses: bool
        :param readouts_per_experiment: readouts per experiment (all will be saved)
        :type readouts_per_experiment: int
        :param start_src: "internal" (tProc starts immediately) or "external" (each soft_avg waits for an external trigger)
        :type start_src: string
        :param progress: If true, displays progress bar
        :type progress: bool
        :param debug: If true, displays assembly code for tProc program
        :type debug: bool
        :returns:
            - iq_list (:py:class:`list`) - list of lists of averaged decimated I and Q data
        """

        buf = super().acquire_decimated(soc, reads_per_rep=readouts_per_experiment, load_pulses=load_pulses, start_src=start_src, progress=progress, debug=debug)
        # move the I/Q axis from last to second-last
        return np.moveaxis(buf, -1, -2)

class RAveragerProgram(QickProgram):
    """
    RAveragerProgram class, for qubit experiments that sweep over a variable (whose value is stored in expt_pts).
    It is an abstract base class similar to the AveragerProgram, except has an outer loop which allows one to sweep a parameter in the real-time program rather than looping over it in software.  This can be more efficient for short duty cycles.
    Acquire gathers data from both ADCs 0 and 1.

    :param cfg: Configuration dictionary
    :type cfg: dict
    """

    def __init__(self, soccfg, cfg):
        """
        Constructor for the RAveragerProgram, calls make program at the end so for classes that inherit from this if you want it to do something before the program is made and compiled either do it before calling this __init__ or put it in the initialize method.
        """
        super().__init__(soccfg)
        self.cfg = cfg
        self.make_program()
        self.reps = cfg['reps']
        self.expts = cfg['expts']
        if "rounds" in cfg:
            self.rounds = cfg['rounds']

    def initialize(self):
        """
        Abstract method for initializing the program and can include any instructions that are executed once at the beginning of the program.
        """
        pass

    def body(self):
        """
        Abstract method for the body of the program
        """
        pass

    def update(self):
        """
        Abstract method for updating the program
        """
        pass

    def make_program(self):
        """
        A template program which repeats the instructions defined in the body() method the number of times specified in self.cfg["reps"].
        """
        p = self

        rcount = 13
        rii = 14
        rjj = 15

        p.initialize()

        p.regwi(0, rcount, 0)

        p.regwi(0, rii, self.cfg['expts']-1)
        p.label("LOOP_I")

        p.regwi(0, rjj, self.cfg['reps']-1)
        p.label("LOOP_J")

        p.body()

        p.mathi(0, rcount, rcount, "+", 1)

        p.memwi(0, rcount, self.counter_addr)

        p.loopnz(0, rjj, 'LOOP_J')

        p.update()

        p.loopnz(0, rii, "LOOP_I")

        p.end()

    def get_expt_pts(self):
        """
        Method for calculating experiment points (for x-axis of plots) based on the config.

        :return: Numpy array of experiment points
        :rtype: array
        """
        return self.cfg["start"]+np.arange(self.expts)*self.cfg["step"]

    def acquire(self, soc, threshold=None, angle=None, load_pulses=True, readouts_per_experiment=1, save_experiments=None, start_src="internal", progress=False, debug=False):
        """
        This method optionally loads pulses on to the SoC, configures the ADC readouts, loads the machine code representation of the AveragerProgram onto the SoC, starts the program and streams the data into the Python, returning it as a set of numpy arrays.
        config requirements:
        "reps" = number of repetitions;

        :param soc: Qick object
        :type soc: Qick object
        :param threshold: threshold
        :type threshold: int
        :param angle: rotation angle
        :type angle: list
        :param readouts_per_experiment: readouts per experiment
        :type readouts_per_experiment: int
        :param save_experiments: saved readouts (by default, save all readouts)
        :type save_experiments: list
        :param load_pulses: If true, loads pulses into the tProc
        :type load_pulses: bool
        :param start_src: "internal" (tProc starts immediately) or "external" (each round waits for an external trigger)
        :type start_src: string
        :param progress: If true, displays progress bar
        :type progress: bool
        :param debug: If true, displays assembly code for tProc program
        :type debug: bool
        :returns:
            - expt_pts (:py:class:`list`) - list of experiment points
            - avg_di (:py:class:`list`) - list of lists of averaged accumulated I data for ADCs 0 and 1
            - avg_dq (:py:class:`list`) - list of lists of averaged accumulated Q data for ADCs 0 and 1
        """
        self.shot_angle = angle
        self.shot_threshold = threshold

        if save_experiments is None:
            save_experiments = range(readouts_per_experiment)

        d_buf, avg_d, shots = super().acquire(soc, reads_per_rep=readouts_per_experiment, load_pulses=load_pulses, start_src=start_src, progress=progress, debug=debug)

        # reformat the data into separate I and Q arrays
        # save results to class in case you want to look at it later or for analysis
        self.di_buf = d_buf[:,:,0]
        self.dq_buf = d_buf[:,:,1]

        if threshold is not None:
            self.shots = shots

        expt_pts = self.get_expt_pts()

        n_ro = len(self.ro_chs)
        avg_di = np.zeros((n_ro, len(save_experiments), self.expts))
        avg_dq = np.zeros((n_ro, len(save_experiments), self.expts))

        for nn, ii in enumerate(save_experiments):
            for i_ch, (ch, ro) in enumerate(self.ro_chs.items()):
                avg_di[i_ch][nn] = avg_d[i_ch, ii, :, 0]
                avg_dq[i_ch][nn] = avg_d[i_ch, ii, :, 1]

        return expt_pts, avg_di, avg_dq


class AbsQickSweep:
    """
    Abstract QickSweep class.
    """

    def __init__(self, prog: QickProgram, label=None):
        """
        :param prog: QickProgram in which the sweep will run.
        :param label: label to be used for the loop tag in qick asm program.
        """
        self.prog = prog
        self.label = label
        self.expts: int = None

    def get_sweep_pts(self) -> Union[List, np.array]:
        """
        abstract method for getting the sweep values
        """
        pass

    def update(self):
        """
        abstract method for updating the sweep value
        """
        pass

    def reset(self):
        """
        abstract method for resetting the sweep value at the beginning of each sweep.
        """
        pass


class QickSweep(AbsQickSweep):
    """
    QickSweep class, describes a sweeps over a qick register.
    """

    def __init__(self, prog: QickProgram, reg: QickRegister, start, stop, expts: int, label=None):
        """

        :param prog: QickProgram in which the sweep happens.
        :param reg: QickRegister object associated to the register to sweep.
        :param start: start value of the register to sweep, in physical units
        :param stop: stop value of the register to sweep, in physical units
        :param expts: number of experiment points from start to stop value.
        :param label: label to be used for the loop tag in qick asm program.
        """
        super().__init__(prog)
        self.reg = reg
        self.start = start
        self.stop = stop
        self.expts = expts
        step_val = (stop - start) / (expts - 1)
        self.step_val = step_val
        self.reg.init_val = start

        if label is None:
            self.label = self.reg.name
        else:
            self.label = label

    def get_sweep_pts(self):
        return np.linspace(self.start, self.stop, self.expts)

    def update(self):
        """
        update the register value. This will be called after finishing last register sweep.
        This function should be overwritten if more complicated update is needed.
        :return:
        """
        self.reg.set_to(self.reg, '+', self.step_val)

    def reset(self):
        """
        reset the register to the start value. will be called at the beginning of each sweep.
        This function should be overwritten if more complicated reset is needed.
        :return:
        """
        self.reg.reset()


def merge_sweeps(sweeps: List[QickSweep]) -> AbsQickSweep:
    """
    create a new QickSweep object that merges the update and reset functions of multiple QickSweeps into one. This is
    useful when multiple registers need to be updated at the same time in one sweep.
    :param sweeps: list of "QickSweep"s
    :return:
    """
    label = "-".join([swp.label for swp in sweeps])

    merged = AbsQickSweep(sweeps[0].prog, label)
    merged.get_sweep_pts = sweeps[0].get_sweep_pts
    expts_ = set([swp.expts for swp in sweeps])
    if len(expts_) != 1:
        raise ValueError(f"all sweeps for merging must have same number of expts, got{expts_}")
    merged.expts = sweeps[0].expts

    def _update():
        for swp in sweeps:
            swp.update()

    def _reset():
        for swp in sweeps:
            swp.reset()

    def _get_sweep_pts():
        sweep_pts = []
        for swp in sweeps:
            sweep_pts.append(swp.get_sweep_pts())
        sweep_pts = np.array(sweep_pts).T
        return sweep_pts

    merged.update = _update
    merged.reset = _reset
    merged.get_sweep_pts = _get_sweep_pts

    return merged


class NDAveragerProgram(QickRegisterManagerMixin, QickProgram):
    """
    NDAveragerProgram class, for experiments that sweep over multiple variables in qick. The order of experiment runs
    follow outer->inner: reps, sweep_n,... sweep_0.

    :param cfg: Configuration dictionary
    :type cfg: dict
    """

    def __init__(self, soccfg, cfg):
        """
        Constructor for the NDAveragerProgram. Make the ND sweep asm commands.
        """
        super().__init__(soccfg)
        self.cfg = cfg
        self.qick_sweeps: List[AbsQickSweep] = []
        self.expts = 1
        self.sweep_axes = []
        self.make_program()
        self.reps = cfg['reps']
        if "soft_avgs" in cfg:
            self.rounds = cfg['soft_avgs']
        if "rounds" in cfg:
            self.rounds = cfg['rounds']

    def initialize(self):
        """
        Abstract method for initializing the program. Should include the instructions that will be executed once at the
        beginning of the qick program.
        """
        pass

    def body(self):
        """
        Abstract method for the body of the program.
        """
        pass

    def add_sweep(self, sweep: AbsQickSweep):
        """
        Add a layer of register sweep to the qick asm program. The order of sweeping will follow first added first sweep.
        :param sweep:
        :return:
        """
        self.qick_sweeps.append(sweep)
        self.expts *= sweep.expts
        self.sweep_axes.append(sweep.expts)

    def make_program(self):
        """
        Make the N dimensional sweep program. The program will run initialize once at the beginning, then iterate over
        all the sweep parameters and run the body. The whole sweep will repeat for cfg["reps"] number of times.
        """
        p = self

        p.initialize()  # initialize only run once at the very beginning

        rcount = 13  # total run counter
        rep_count = 14  # repetition counter

        n_sweeps = len(self.qick_sweeps)
        if n_sweeps > 7:  # to be safe, only register 15-21 in page 0 can be used as sweep counters
            raise OverflowError(f"too many qick inner loops ({n_sweeps}), run out of counter registers")
        counter_regs = (np.arange(n_sweeps) + 15).tolist()  # not sure why this has to be a list (np.array doesn't work)

        p.regwi(0, rcount, 0)  # reset total run count

        # set repetition counter and tag
        p.regwi(0, rep_count, self.cfg["reps"] - 1)
        p.label("LOOP_rep")

        # add reset and start tags for each sweep
        for creg, swp in zip(counter_regs[::-1], self.qick_sweeps[::-1]):
            swp.reset()
            p.regwi(0, creg, swp.expts - 1)
            p.label(f"LOOP_{swp.label if swp.label is not None else creg}")

        # run body and total_run_counter++
        p.body()
        p.mathi(0, rcount, rcount, "+", 1)
        p.memwi(0, rcount, 1)

        # add update and stop condition for each sweep
        for creg, swp in zip(counter_regs, self.qick_sweeps):
            swp.update()
            p.loopnz(0, creg, f"LOOP_{swp.label if swp.label is not None else creg}")

        # stop condition for repetition
        p.loopnz(0, rep_count, 'LOOP_rep')

        p.end()

    def get_expt_pts(self):
        """
        :return:
        """
        sweep_pts = []
        for swp in self.qick_sweeps:
            sweep_pts.append(swp.get_sweep_pts())
        return sweep_pts

    def acquire(self, soc, threshold: int = None, angle: List = None, load_pulses=True, readouts_per_experiment=1,
                save_experiments: List = None, start_src: str = "internal", progress=False, debug=False):
        """
        This method optionally loads pulses on to the SoC, configures the ADC readouts, loads the machine code
        representation of the AveragerProgram onto the SoC, starts the program and streams the data into the Python,
        returning it as a set of numpy arrays.
        Note here the buf data has "reps" as the outermost axis, and the first swept parameter corresponds to the
        innermost axis.

        config requirements:
        "reps" = number of repetitions;

        :param soc: Qick object
        :param threshold: threshold
        :param angle: rotation angle
        :param readouts_per_experiment: readouts per experiment
        :param save_experiments: saved readouts (by default, save all readouts)
        :param load_pulses: If true, loads pulses into the tProc
        :param start_src: "internal" (tProc starts immediately) or "external" (each round waits for an external trigger)
        :param progress: If true, displays progress bar
        :param debug: If true, displays assembly code for tProc program
        :returns:
            - expt_pts (:py:class:`list`) - list of experiment points
            - avg_di (:py:class:`list`) - list of lists of averaged accumulated I data for ADCs 0 and 1
            - avg_dq (:py:class:`list`) - list of lists of averaged accumulated Q data for ADCs 0 and 1
        """

        self.shot_angle = angle
        self.shot_threshold = threshold

        if save_experiments is None:
            save_experiments = range(readouts_per_experiment)

        # avg_d calculated in QickProgram.acquire() assumes a different data shape, here we will recalculate based on
        # the d_buf returned.
        d_buf, avg_d, shots = super().acquire(soc, reads_per_rep=readouts_per_experiment, load_pulses=load_pulses,
                                              start_src=start_src, progress=progress, debug=debug)

        # reformat the data into separate I and Q arrays
        # save results to class in case you want to look at it later or for analysis
        self.di_buf = d_buf[:, :, 0]
        self.dq_buf = d_buf[:, :, 1]

        if threshold is not None:
            self.shots = shots

        expt_pts = self.get_expt_pts()

        n_ro = len(self.ro_chs)
        avg_di = np.zeros((n_ro, len(save_experiments), *self.sweep_axes[::-1]))
        avg_dq = np.zeros((n_ro, len(save_experiments), *self.sweep_axes[::-1]))

        for nn, ii in enumerate(save_experiments):
            for i_ch, (ch, ro) in enumerate(self.ro_chs.items()):
                avg_di[i_ch][nn] = avg_d[i_ch, ii, ..., 0]
                avg_dq[i_ch][nn] = avg_d[i_ch, ii, ..., 1]

        return expt_pts, avg_di, avg_dq

    def _average_buf(self, d_reps, reads_per_rep: int):
        """
        overwrites the default _average_buf method in QickProgram. Here "reps" is the outermost axis, and we reshape
        avg_d to the shape of the sweep axes.
        :param d_reps:
        :param reads_per_rep:
        :return:
        """
        avg_d = np.zeros((len(self.ro_chs), reads_per_rep, *self.sweep_axes[::-1], 2))
        for ii in range(reads_per_rep):
            for i_ch, (ch, ro) in enumerate(self.ro_chs.items()):
                avg_d[i_ch][ii] = np.sum(d_reps[i_ch][ii::reads_per_rep, :].reshape((self.reps, *self.sweep_axes[::-1], 2)),
                                         axis=0) / (self.reps * ro['length'])
        return avg_d
