import numpy as np
import warnings
from utility.scheduler import scheduler

def config_awg_schedule(awg, settings):
    """
    Helper function to initialize the scheduler object and load the buffer
    program to the AWG. Returns the scheduler object to the user

    Args:
        awg (HDAWG.AWG): awg core of the HDAWG
        settings (dict): list of settings for the measurement pulse (exp_globals['m_pulse'])

    Returns:
        scheduler: configured scheduler object
    """
    start_time = settings['meas_pos']
    window_time = settings['meas_window']
    awg_sched = scheduler(total_time=start_time+2*window_time, sample_rate=2.4e9)
    awg_sched.add_analog_channel(1, name='Qubit_I')
    awg_sched.add_analog_channel(2, name='Qubit_Q')

    awg_sched.add_digital_channel(1, name='Qubit_enable', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
    awg_sched.add_digital_channel(2, name='Cavity_enable', polarity='Pos', HW_offset_on=0, HW_offset_off=0)

    progFile = open(r"C:\Users\Kollarlab\Desktop\Kollar-Lab\pulsed_measurements\HDAWG_sequencer_codes\hdawg_placeholder.cpp",'r')
    rawprog  = progFile.read()
    loadprog = rawprog
    progFile.close()
    loadprog = loadprog.replace('_samples_', str(awg_sched.samples))
    awg.load_program(loadprog)

    return awg_sched

def load_pulse_sequence(awg, scheduler):
    """
    Upload programmed sequence to the AWG and start it

    Args:
        awg (HDAWG.AWG): awg core of the HDAWG
        scheduler (scheduler): pulse schedule object
    """

    awg.stop()
    [ch1, ch2, marker] = scheduler.compile_schedule('HDAWG', ['Qubit_I', 'Qubit_Q'], ['Qubit_enable', 'Cavity_enable'])
    awg.load_waveform('0', ch1, ch2, marker)
    awg.run_loop()

def add_gate_analog(scheduler, amp, angle, position, q_pulse):
    """
    Low level function that computes the IQ traces to generate the
    desired gate. It also adds an enable window "bracketing" the pulse.
    We added the capability to apply corrections to the IQ signal
    depending on the input angles (essentially removing remaining IQ
    imbalance in the system from HDAWG or SGS or something else)

    Args:
        scheduler (scheduler): pulse schedule object
        amp (float): scaling parameter for pulse amplitude (1 corresponds to a pi/2 pulse)
        angle (float): angle in XY plane of the drive (RADIANS)
        position (float): position of the leading edge of the gate
        q_pulse (dict): dictionary holding all the qubit pulse calibration params
    """
    imp_off = q_pulse['imp_off']
    imp_amp = q_pulse['imp_amp']
    imp_phi = q_pulse['imp_phi']
    p_amp   = q_pulse['p_amp']
    p_amp = amp*p_amp*imp_off/(imp_off+imp_amp*np.sin(angle+imp_phi))
    if p_amp>1:
        warnings.warn('HDAWG amplitude is > 1, check settings')
    I = p_amp*np.cos(angle)
    Q = p_amp*np.sin(angle)
    qubit_I = scheduler.analog_channels['Qubit_I']
    qubit_Q = scheduler.analog_channels['Qubit_Q']
    qubit_marker = scheduler.digital_channels['Qubit_enable']

    pulse_length = q_pulse['num_sigma']*q_pulse['sigma']+q_pulse['hold_time']
    buffer = q_pulse['delay']
    qubit_I.add_pulse('gaussian_square', position=position, 
                          amplitude=I, length = q_pulse['hold_time'], 
                          ramp_sigma=q_pulse['sigma'], num_sigma=q_pulse['num_sigma'])
    qubit_Q.add_pulse('gaussian_square', position=position, 
                          amplitude=Q, length = q_pulse['hold_time'], 
                          ramp_sigma=q_pulse['sigma'], num_sigma=q_pulse['num_sigma'])
    qubit_marker.add_window(position-buffer, position+pulse_length+buffer)

def add_gate(scheduler, gate, position, q_pulse):
    """
    Higher level function, takes in a string description of the gate 
    to apply and computes the appropriate IQ values to apply the gate
    Will need to be expanded as we want to try other gates 

    Args:
        scheduler (scheduler): pulse schedule object
        gate (string): string description of gate to apply
        position (float): position of the leading edge of the gate
        q_pulse (dict): dictionary holding all the qubit pulse calibration params
    """
    if gate=='X':
        amp = 1
        angle = 0
    if gate=='Y':
        amp = 1
        angle = np.pi/2
    if gate=='-X':
        amp = 1
        angle = np.pi
    if gate=='-Y':
        amp = 1
        angle = 3*np.pi/2
    if gate=='Z':
        amp = 0
        angle = 0
    if gate=='PI':
        amp = 2
        angle = 0
    else:
        try:
            data = eval(gate)
            amp = data[0]
            angle = data[1]
        except NameError:
            NameError('Unsupported gate')
    add_gate_analog(scheduler, amp, angle, position, q_pulse)

def add_meas_window(scheduler, m_pulse):
    """
    Measurement window function, creates the marker window 

    Args:
        scheduler (scheduler): pulse schedule object
        m_pulse (dict): dictionary holding the measurement window settings (length and position)
    """
    window_time = m_pulse['meas_window']
    window_pos  = m_pulse['meas_pos']
    cavity_marker = scheduler.digital_channels['Cavity_enable']
    cavity_marker.add_window(window_pos, window_pos+window_time+1e-6)

def state_prep_and_measure(scheduler, pulse_seq, meas_basis, settings):
    """
    General purpose function that creates the pulse schedule to 
    generate the desired state and measure it in the specified basis. 
    This function will reset all the channels (IQ, enable and cavity)

    Args:
        scheduler (scheduler): object holding the pulse schedule (lives in the utility folder)
        pulse_seq (list of strings): list of gates to apply to the qubit (e.g: 'X', '-Y', 'Z', etc.)
        meas_basis (string): string defining measurement basis ('X', 'Y', 'Z')
        settings (dict): dictionary of settings configuring the underlying pulses (sigma, hold length, meas time etc.)
    """
    scheduler.reset()

    q_pulse = settings['qubit_pulse']
    m_pulse = settings['measurement_pulse']

    meas_pos = m_pulse['meas_pos']
    num_pulse = len(pulse_seq)+1 #Accounts for the basis pulse at the end of the sequence
    pulse_length = q_pulse['num_sigma']*q_pulse['sigma']+q_pulse['hold_time']
    buffer = q_pulse['delay']

    position = meas_pos-(buffer+pulse_length)*num_pulse
    if position<buffer:
        warnings.warn('Pulse sequence too long for this measurement position')

    for gate in pulse_seq:
        add_gate(scheduler, gate, position, q_pulse)
        position+=buffer+pulse_length
    add_gate(scheduler, meas_basis, position, q_pulse)
    add_meas_window(scheduler, m_pulse)
