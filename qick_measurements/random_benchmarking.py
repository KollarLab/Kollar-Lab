# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 11:52:24 2024

@author: Ruthie Vogel
"""

from qick.averager_program import AveragerProgram

import numpy as np
import os
import matplotlib.pyplot as plt
import userfuncs
import cmath
import random
import time


class IQProgram_RB(AveragerProgram):
    
    '''
    This class is optimized for Ruthie's randomized benchmarking and state tomography thesis project. Once pulses are defined by the user (and given names) it can automatically play any string of pulses separated
    by an empirical buffer defined in the config file. 
    '''
    
    def initialize(self):
        cfg     = self.cfg   
        gen_ch  = cfg["cav_channel"]
        qub_ch  = cfg["qub_channel"]
        ro_chs  = cfg["ro_channels"]
        
        
        # set the nyquist zone
        self.declare_gen( ch = cfg["cav_channel"], nqz = cfg["nqz_c"])
        self.declare_gen( ch = cfg["qub_channel"], nqz = cfg["nqz_q"])
        
        
        # configure the readout lengths and downconversion frequencies (ensuring it is an available DAC frequency)
        for ch in cfg["ro_channels"]:
            self.declare_readout( ch = ch, length = self.us2cycles( cfg["readout_length"], ro_ch = ro_chs[0]),
                                 freq=self.cfg["cav_freq"], gen_ch = gen_ch)

        # convert frequency to DAC frequency (ensuring it is an available ADC frequency)
        freq_c  = self.freq2reg( cfg["cav_freq"], gen_ch = gen_ch, ro_ch = ro_chs[0])
        phase_c = self.deg2reg( cfg["cav_phase"], gen_ch = gen_ch)
        gain_c  = cfg["meas_gain"]

        self.default_pulse_registers(ch = gen_ch, freq = freq_c, phase = phase_c, gain = gain_c, mode = "oneshot")
        self.set_pulse_registers(ch=gen_ch, style="const", length=self.us2cycles(self.cfg["meas_window"],gen_ch=gen_ch))



        freq_q  = self.freq2reg( cfg["qub_freq"], gen_ch = qub_ch) #, ro_ch = ro_chs[0])
        phase_q = self.deg2reg( cfg["qub_phase"], gen_ch = gen_ch) #should this be gen_ch? 
        gain_q  = cfg["qub_gain"]
        
        self.default_pulse_registers(ch = qub_ch, freq = freq_q, phase = phase_q, gain = gain_q, mode = "oneshot")

        # loads each unique pulse from 'pulse_data' into memory
        for pulse in cfg['pulse_data']:
            self.add_pulse( ch = qub_ch, name = pulse[0], idata = pulse[1]*int(self.soccfg['gens'][qub_ch]['maxv']), 
                           qdata = pulse[2]*int(self.soccfg['gens'][qub_ch]['maxv']) )
        
        
        self.synci(200)  # give processor some time to configure pulses
    
    
  
    def body(self):
       
        '''
        The body function sets pulse registers for each pulse and plays the pulses automatically, then increments the time
        to play the next pulse by the length of this pulse plus an empirical buffer defined in the config dicionary
        '''
        cfg     = self.cfg
        gen_ch  = cfg["cav_channel"]
        qub_ch  = cfg["qub_channel"]
        print(cfg["pulse_schedule"])
        
        sigma = self.us2cycles(self.cfg["qub_sigma"])
        num_sigma = self.cfg["num_sigma"]
        
        offset = self.us2cycles( cfg["adc_trig_offset"], gen_ch = gen_ch )
        meas_time = self.us2cycles(self.cfg["meas_time"],gen_ch=gen_ch)


        self.trigger(adcs = self.ro_chs, 
                     pins = [0], 
                     adc_trig_offset = offset )
      
        flip_pulses = np.flip(np.array(cfg["pulse_schedule"]))

        pulse_len = int(len(cfg['pulse_data'][0][1])/16)
        buffer = self.us2cycles( cfg["buffer"], gen_ch = qub_ch) 
        playtime = meas_time - ( self.us2cycles( cfg["qub_delay"], gen_ch = qub_ch) 
                                + pulse_len*( len( cfg['pulse_schedule'] ) + 1 ) 
                                + buffer*len( cfg['pulse_schedule'] ) )

        plot_pulse_I = []
        plot_pulse_Q = []
        
        '''
        for item in flip_pulses:
            
            self.set_pulse_registers( qub_ch, style = 'arb', waveform = item )
            self.pulse( ch = qub_ch, t = playtime )
            #print( "we have played the pulse {} at t = {}".format( item, playtime ) )
            
            #print(self.cycles2us( (38320-32128)/16, gen_ch = qub_ch))

            for i in range( len ( cfg["pulse_data"] ) ):
                if item == cfg["pulse_data"][i][0]:
                    gate = cfg["pulse_data"][i]
                
                    addtime = buffer*16 + int( len( gate[1] )/16 )
                    playtime = playtime - addtime  
                    
                    buff_zero = [0]*(buffer*16)
                    
                    # Now let's add the pulse data to the arrays to plot
                    plot_pulse_I.extend(gate[1].tolist())
                    plot_pulse_I.extend(buff_zero)
                    plot_pulse_Q.extend(gate[2].tolist())
                    plot_pulse_Q.extend(buff_zero)
        '''

        for item in cfg['pulse_schedule']:
            self.set_pulse_registers( qub_ch, style = 'arb', waveform = item )
            self.pulse( ch = qub_ch, t = playtime )
            
            for i in range( len ( cfg["pulse_data"] ) ):
                if item == cfg["pulse_data"][i][0]:
                    gate = cfg["pulse_data"][i]
                
                    addtime = buffer + int( len( gate[1] )/16 )
                    playtime = playtime + addtime  
                    
                    buff_zero = [0]*(buffer*16)
                    
                    # Now let's add the pulse data to the arrays to plot
                    #print('adding gate{} to plot'.format(gate[0]))
                    plot_pulse_I.extend(gate[1].tolist())
                    plot_pulse_I.extend(buff_zero)
                    plot_pulse_Q.extend(gate[2].tolist())
                    plot_pulse_Q.extend(buff_zero)


        cfg['pulse_plot'] = [plot_pulse_I, plot_pulse_Q]
        
        
        self.pulse(ch=gen_ch,t=meas_time)

        self.wait_all()
        self.sync_all( self.us2cycles( self.cfg["relax_delay"] ) )
        
        
class gaussian_square():
    '''
    gaussian_square is the main class that creates the IQ data for the QICKboard based pulses.

    **Returns:**
        ``pulse.I``     : (list)  containing the In-Phase waveform for the pulse
        ``pulse.Q``     : (list)  containing the Quadrature waveform for the pulse
        ``pulse.rotate``: (array) this numpy array is the rotation matrix that this pulse applies to the qubit

    '''
    _defaults = {
        'sigma'      : 5,
        'num_sigma'  : 4,
        'hold_time'  : 0,     
        'amp'        : 1,
        'sample_rate': 430.080e6,
        'angle'      : 0,
        'rot_amount' : np.pi,
        'axis'       : 'x'
    }
    
    def __init__(self, pulse, settings):
        self._settings = {**self._defaults, **settings}
        self.set_settings()
        self.pulse = self.compile_pulse()
        self.rmatrix = self.rotation()

        
    def set_settings(self):
        for parameter in self._settings.keys():
            if parameter not in self._defaults.keys():
                print('Error or typo, {} is not a valid parameter for this class'.format(parameter))
            else: self.__setattr__(parameter, self._settings[parameter])
                
    @property
    def settings(self):
        full_settings = {}
        for param in self._defaults.keys():
            full_settings[param] = getattr(self, param)
        return full_settings
    
    def compile_pulse(self):
        '''
        This function is the main part of the pulse creation. It returns the I data and Q data that we need for our pulses. 
        '''
        
        amp         = self.amp
        sigma       = self.sigma/1e6
        length      = self.hold_time/1e6
        num_sigma   = self.num_sigma
        sample_rate = self.sample_rate # multiplied by 16 because each clock tick is broken up into 16 pieces.
        
        #create gaussian ramp
        
        samples = int(sigma*num_sigma*sample_rate)*16
        t = np.linspace(0, sigma*num_sigma, samples)
        t0 = num_sigma*sigma/2
        pulse = np.exp(-(t-t0)**2/(sigma**2)) # correct form for a gaussian has 2*(sigma**2) but the qick board is wack
        offset = pulse[0]
        init_amp = max(pulse)
        ramp = amp*(pulse-offset)/(init_amp-offset)
        #hold max amplitude
        square = np.ones(int(length*sample_rate))*amp
        ramp_up, ramp_down = np.split(ramp, 2)
        final = np.concatenate((ramp_up, square, ramp_down))
        
        self.pulse = final
        self.I = np.cos(self.angle)*final
        self.Q = np.sin(self.angle)*final
    
    def rotation(self):
        
        '''
        This function assesses the rotation that each pulse applies to the qubit, which then gets returned as {}.rotate
        This can also be used for any pulse regardless of the axis. 
        '''
        axis = self.axis
        theta = self.rot_amount
        
        accept_axes = ['x', 'y', 'z']
        
        if axis == 'x':
            rot_matrix = np.array([[np.cos(theta/2), -1j*np.sin(theta/2)], [-1j*np.sin(theta/2), np.cos(theta/2)]])
            
        elif axis == 'y':
            rot_matrix = np.array([[np.cos(theta/2), -np.sin(theta/2)], [np.sin(theta/2), np.cos(theta/2)]])

        elif axis not in accept_axes:
            print('Error or typo, {} is not a valid axis for this function'.format(axis))
            
        if theta == np.pi:
            rot_matrix = rot_matrix*1j

        
        
        # This section makes sure that everything gets rounded away nicely
        rot_matrix[np.abs(rot_matrix) <= np.finfo(np.float64).eps] = 0
        
        
        self.rmatrix = rot_matrix
        self.rotate = rot_matrix   
        

def random_seq_gen(config_dict):

    '''
    This function generates the RB schedule, both the total schedule and the truncated schedule for a given run.
    It also calculates and then adds the required pulse to return an ideal qubit back to the ground state at the end of the schedule.

    **Parameters:**
        ``config_dict``: the config dictionary for the run, composed of exp_settings, exp_globals, and the locally defined factors

        ``gen_num``: The number of new schedules to be created.

        ``num_gates``: the total number of gates in the RB sequence.

    **Returns:**
        ``full_seq_dict``: (dictionary) contains all the sequences for the run, with numbered keys.        
    '''

    cfg = config_dict
    gen_new = cfg["gen_new"]
    
    if gen_new == True:
        gen_num = cfg["gen_num"] 

        n = cfg["num_gates"]
        full_seq_dict = {}

        for seq in range(gen_num):
            
            rand_list = np.zeros(n)
            full_seq = list(map(str, np.zeros(n)))

            for i in range(n):
                rand_list[i] = random.randint(0,3)
                placeholder = int(rand_list[i])
                full_seq[i] = (cfg["pulse_data"][placeholder][0])
        
            full_seq_dict[seq] = full_seq

    else:
        full_seq_dict = config_dict['full_seq_dict']

    return full_seq_dict

def meas_pulse(config_dict):
        
    cfg = config_dict
    
    '''
    This section takes the density matrix of the initial state and operates on it with the rotation matrices of each 
    pulse applied in the benchmarking sequence defined in 'pulse_schedule'. It then calculates what rotation is necessary
    to return the qubit to the ground state.

    This is all done assuming that each gate and the initial state are theoretical and perfect.

    Once that calculation is done, it'll generate the new pulse and append it to the pulse schedule.
    '''


    # specify the state we expect our qubit to start in: for us, it's the ground state. all rotations are then applied
    init_state = [[1], [0]]
    state = init_state

    print("the current schedule is {}".format(cfg["pulse_schedule"]))
    
    
    for pulse in cfg["pulse_schedule"]:
        for i in range( len ( cfg["pulse_data"] ) ):
            if pulse == cfg["pulse_data"][i][0]:
                new_state = cfg["pulse_data"][i][3] @ state
                state = new_state

    '''
    Throughout this function, there are a number of times where we have to decompose the state or rotation matrix and
    zero out all the rounding errors that come from python doing trig or working with complex numbers because otherwise it
    very thoroughly messes everything else up.
    '''
    
    for count, element in enumerate(state):
        if np.isclose(element.real, [0]):
            state[count].real = element.real*0
        if np.isclose(element.imag, [0]):
            state[count].imag = element.imag*0

    '''
    The following section puts the state we're left with in the conventional format that quantum states are written out in.
    This ensures that the operations we're doing on it are doing the right thing!
    '''
    
    if np.imag(state[0]) != 0:
        conj = np.conjugate(state[0])
        for count, element in enumerate(state):
            state[count] = element*conj*np.sqrt(2)
            
    # Gotta zero out after taking the complex conjugate
    for count, element in enumerate(state):
        if np.isclose(element.real, [0]):
            state[count].real = element.real*0
        if np.isclose(element.imag, [0]):
            state[count].imag = element.imag*0

    if state[0] <0:
        state = state * -1        
    
    
    '''
    Every point on the bloch sphere is described by theta and phi. The following part uses our state to pull out those values.
    Right now phi is slightly hard coded and only works for the specific points that we can get to on the bloch sphere.
    TODO: calculate phi empirically!
    '''
    theta = cmath.acos(state[0].real)

    for angle in [0, np.pi, np.pi/2, -np.pi/2]:
        if np.isclose(theta, angle):
            theta = angle
    
    theta = -1 * theta #because we need to rotate the opposite way

    phi = 0

    if state[1].imag == 0:
        if state[1].real > 0:
            phi = 0
        elif state[1].real < 0:
            phi = np.pi

    elif state[1].imag > 0:
        phi = np.pi/2

    elif state[1].imag < 0:
        phi = -np.pi/2


    '''
    This section calculates the rotation matrix of the pulse needed to get us back to the ground state. It then sorts through
    the pulse dictionary and appends the right pulse to the pulse schedule
    '''
    
    rot_matrix = np.array([[np.cos(theta), -np.exp(-1j*phi)*np.sin(theta)], [np.exp(1j*phi)*np.sin(theta), np.cos(theta)]])

    # Gotta zero out again
    for i in range(len(rot_matrix)):
        sublist = np.array(rot_matrix[i])
        for count, element in enumerate(sublist):
            if np.isclose(element.real, [0]):
                real = 0
            else:
                real = element.real

            if np.isclose(element.imag, [0]):
                imag = 0
            else:
                imag = element.imag
            sublist[count] = real + imag*1j
        rot_matrix[i] = sublist


    if state[1] == 0:
        meas_pulse = np.array([[1, 0], [0, 1]])

    if state[0] == 0:
        meas_pulse = np.array([[0, 1], [1, 0]])
    else:
        meas_pulse = rot_matrix 

    for i in range( len ( cfg["pulse_data"] ) ):
        subtract = meas_pulse - cfg["pulse_data"][i][3]
        zero_array = np.zeros((2, 2))
        if np.allclose(subtract, zero_array):
            cfg["pulse_schedule"].append(cfg["pulse_data"][i][0])
            print('added {}'.format(cfg["pulse_data"][i][0]))
    
    return config_dict



def get_RB_settings():
    settings = {}
    
    settings['scanname']  = 'initial_RB_test'
    settings['meas_type'] = 'RandomizedBenchmarking'
    
    settings['cav_freq'] = 1e9
    settings['cav_gain'] = 1000
    
    # RB parameters
    settings['pulse_schedule'] = []
    settings['pulse_data']     = [] 
    settings['pulse_plot']     = []
    settings['full_seq_dict']  = {}
    settings['buffer']         = .1 #[us]
    settings['gen_new']        = False
    settings['run_full']       = False # only True for initial scans 
    settings['gen_num']        = 2

    settings['num_gates']      = 10
    settings['used_gates']     = 1
   
    #Card settings
    settings['reps'] = 1
    settings['soft_avgs'] = 1
    
    return settings

def random_bench(soc,soccfg,instruments,settings):
    exp_globals  = settings['exp_globals']
    exp_settings = settings['exp_settings'] 
    m_pulse      = exp_globals['measurement_pulse']
    q_pulse      = exp_globals['qubit_pulse']
    pauli_gates  = exp_globals['pauli_gates']
    
    soc.reset_gens()
    
    #if exp_globals['LO']:
    if False:
        logen = instruments['LO']
        
        logen.freq   = exp_globals['LO_freq']
        logen.power  = exp_globals['LO_power']
        logen.output = 1
    else:
        pass

    stamp    = userfuncs.timestamp()
    saveDir  = userfuncs.saveDir(settings)
    filename = exp_settings['scanname'] + '_' + stamp

    config = {
        'cav_channel'     : exp_globals['cav_channel'],
        'qub_channel'     : exp_globals['qub_channel'],
        'ro_channels'     : exp_globals['ro_channels'],

        'nqz_c'           : 1,
        'cav_phase'       : m_pulse['cav_phase'],
        'meas_window'     : m_pulse['meas_window'],
        'meas_time'       : m_pulse['meas_pos'],
        'meas_gain'       : exp_settings['cav_gain'],
        'cav_freq'        : (exp_settings['cav_freq']-exp_globals['LO_freq'])/1e6,
        
        'nqz_q'           : 2,
        'qub_phase'       : q_pulse['qub_phase'],
        'qub_freq'        : exp_settings['qub_freq']/1e6,
        'qub_gain'        : exp_settings['qub_gain'],
        'qub_sigma'       : q_pulse['sigma'],
        'qub_delay'       : q_pulse['delay'],
        'num_sigma'       : q_pulse['num_sigma'],
        'hold_time'       : q_pulse['Hold_time'],
        
        # RB specialized settings
        'pulse_data'      : exp_settings['pulse_data'],
        'pulse_schedule'  : exp_settings['pulse_schedule'],
        'full_seq_dict'   : exp_settings['full_seq_dict'],
        'buffer'          : exp_settings['buffer'],
        'pulse_plot'      : exp_settings['pulse_plot'],
        'gen_new'         : exp_settings['gen_new'],
        'run_full'        : exp_settings['run_full'],
        'gen_num'         : exp_settings['gen_num'],
        'num_gates'       : exp_settings['num_gates'],
        'used_gates'      : exp_settings['used_gates'],
        
        'readout_length'  : m_pulse['init_buffer'] + m_pulse['meas_window'] + m_pulse['post_buffer'],
        'adc_trig_offset' : m_pulse['emp_delay'] + m_pulse['meas_pos'] - m_pulse['init_buffer'],


        'relax_delay'     : exp_globals['relax_delay'],
        'reps'            : exp_settings['reps'],
        'soft_avgs'       : exp_settings['soft_avgs']
        }
    

    # The following section sets up the I and Q data for each of the standard pulses
    
    x_pulse = pauli_gates['x_pulse']
    y_pulse = pauli_gates['y_pulse']
    x2_pulse = pauli_gates['x2_pulse']
    y2_pulse = pauli_gates['y2_pulse']
    x2_neg = pauli_gates['x2_neg']
    y2_neg = pauli_gates['y2_neg']
    ident_pulse = pauli_gates['identity']

    pulse_set = [x_pulse, y_pulse, x2_pulse, y2_pulse, x2_neg, y2_neg, ident_pulse]

    pulse_list_RB = []

    for pulse in pulse_set:
        set_dict = {
            'sigma'      : q_pulse['sigma'], 
            'num_sigma'  : q_pulse['num_sigma'],
            'hold_time'  : q_pulse['Hold_time'],     
            'amp'        : pulse[1],
            'sample_rate': 430.080e6,
            'angle'      : pulse[3],
            'rot_amount' : pulse[2],
            'axis'       : pulse[4]
        }
        
        pulse_comp = gaussian_square(pulse = np.zeros(120), settings = set_dict)
        pulse_data = [pulse[0], pulse_comp.I, pulse_comp.Q, pulse_comp.rotate]

        pulse_list_RB.append(pulse_data)
 
    config['pulse_data'] = pulse_list_RB


    # This next section sets the sequences of pulses and iterates through them

    #full_seq_dict = random_seq_gen(config) TODO: Reinstate after debugging

    #config['full_seq_dict'] = full_seq_dict TODO: Reinstate after debugging

    full_seq_dict = config['full_seq_dict']

    power_data = np.zeros((len(full_seq_dict), len(config['used_gates'])))
    phase_data = np.zeros((len(full_seq_dict), len(config['used_gates'])))
    I_finals   = np.zeros((len(full_seq_dict), len(config['used_gates'])))
    Q_finals   = np.zeros((len(full_seq_dict), len(config['used_gates'])))
    # TODO: figure out if I want to save all of the IQ data from each run, that might be too much
    # If I wanted to make it possible I'd have to use dictionaries

    for seq in full_seq_dict:
        run_seq = full_seq_dict[seq]
        print('Running Sequence {}'.format(seq))                

        for i in range(len(config['used_gates'])):
            if i == 0:
                tstart = time.time()
            else:
                pass

            trunc_point = config['used_gates'][i]
            config['pulse_schedule'] = run_seq[: trunc_point]

            # Append the measurement pulse
            #config = meas_pulse(config)

            prog_RB = IQProgram_RB(soccfg, config)

            meas_start = prog_RB.us2cycles(m_pulse["init_buffer"],ro_ch=0)
            meas_end = meas_start+prog_RB.us2cycles(m_pulse["meas_window"],ro_ch=0)

            avg_di, avg_dq = prog_RB.acquire(soc, load_pulses=True, progress=False)

            # I and Q full are the entire data for the output traces of the cavity
            I_full = avg_di[0][0]
            Q_full = avg_dq[0][0]

            '''
            This section was taken out when the script switched to acquire from acquire_decimated

            # I and Q final are the average I and Q measurements for each point
            # because RB is largely immune to SPAM errors, it doesn't matter how we average
            I_final = np.mean(I_window)
            Q_final = np.mean(Q_window)
            '''
            # Save the mean IQ data for every sequence and truncate point to an array to look at later
            I_finals[seq, i] = I_full
            Q_finals[seq, i] = Q_full
            
            power_data[seq, i] = np.sqrt(I_full**2 + Q_full**2)
            phase_data[seq, i] = np.arctan2(Q_full, I_full)*180/np.pi
        
        if exp_settings['debug']:
            full_IQ_data = {}

            # This creates a 3D array for each sequence, with IQ x number of truncate points x the length of a measurement array
            IQ_per_seq = np.zeros((2, len(config['used_gates'], 829)))
            full_IQ_data['sequence{}'.format(seq)] = IQ_per_seq

            config['readout_length'] = 2.2
            config['adc_trig_offset'] = m_pulse['emp_delay'] + m_pulse['meas_pos'] + exp_settings['debug_time']

            I_finals_debug   = np.zeros((len(full_seq_dict), len(config['used_gates'])))
            Q_finals_debug   = np.zeros((len(full_seq_dict), len(config['used_gates'])))
            power_debug = np.zeros((len(full_seq_dict), len(config['used_gates'])))
            phase_debug = np.zeros((len(full_seq_dict), len(config['used_gates'])))

            iq_RB_debug = prog_RB.acquire_decimated(soc, load_pulses=True, progress=False, debug=False)

            # I and Q full are the entire data for the output traces of the cavity
            I_full_debug = iq_RB_debug[0][0]
            Q_full_debug = iq_RB_debug[0][1]

            # I and Q window select for only the time when the cavity pulse was happening
            I_window = I_full_debug[meas_start:meas_end]
            Q_window = Q_full_debug[meas_start:meas_end]


            IQ_per_seq[0][i] = I_full_debug
            IQ_per_seq[1][i] = Q_full_debug
            
            # I and Q final are the average I and Q measurements for each point
            # because RB is largely immune to SPAM errors, it doesn't matter how we average
            I_final = np.mean(I_window)
            Q_final = np.mean(Q_window)

            # Save the mean IQ data for every sequence and truncate point to an array to look at later
            I_finals_debug[seq, i] = I_final
            Q_finals_debug[seq, i] = Q_final

            power_debug[seq, i] = np.sqrt(I_final**2 + Q_final**2)
            phase_debug[seq, i] = np.arctan2(Q_final, I_final)*180/np.pi

            # Define an x axis the length of our IQ data and convert it to us instead of clock ticks    
            if seq == 0: 
                xaxis = np.linspace( 0, len( iq_RB_debug[0][0] ) -1, len( iq_RB_debug[0][0] ) )
                for x in range(0,len(xaxis)):
                    xaxis[x] = prog_RB.cycles2us(xaxis[x],ro_ch=1)
                
                x_axis_pulses = np.linspace(0, len(config['pulse_plot'][0])-1, len(config['pulse_plot'][0]))
                for x in range(0,len(x_axis_pulses)):
                    x_axis_pulses[x] = prog_RB.cycles2us(x_axis_pulses[x],gen_ch=6)
                x_axis_pulses = x_axis_pulses/16

            # Figure 1 plots single shot I and Q output cavity traces and Input traces
            fig3 = plt.figure(3, figsize=(13, 8))
            plt.clf()
            plt.subplot(121)
            plt.plot(xaxis, I_full_debug, label = 'I Data')
            plt.plot(xaxis, Q_full_debug, label = 'Q Data')
            plt.xlabel('us')
            plt.ylabel('a.u.')
            plt.legend()
            plt.subplot(122)
            plt.plot(x_axis_pulses, config['pulse_plot'][0], label = 'I Data')
            plt.plot(x_axis_pulses, config['pulse_plot'][1], label = 'Q Data')
            plt.xlabel('us')
            plt.ylabel('amplitude')
            plt.legend()
            plt.title('Single Shot Cavity Output/Input Traces')
            fig3.canvas.draw()
            fig3.canvas.flush_events()
            plt.savefig(os.path.join(saveDir, filename+'_Input_Output_Traces.png'), dpi = 150)

        '''
        Woo the plotting section!
        
        The debug plot (above) displays the output and input traces (cavity pulse data)
        Plot 2 displays the power v number of gates (will hopefully be updated from power)
        Plot 3 shows I/Q data for each one of the points
        '''

        # TODO: Find some way to determine fidelity and switch to using that instead of power

        # Figure 2 plots the average power v number of gates for each sequence
        fig1 = plt.figure(1, figsize=(13, 8))
        plt.clf()
        for i in range(0, seq +1):
            plt.scatter(config['used_gates'], power_data[i], label = 'Sequence {}'.format(i))
        plt.title('Power v Num Used Gates\n' + filename)
        plt.xlabel('Num Gates')
        plt.ylabel('Power (DAC a.u.)')
        #plt.legend()
        fig1.canvas.draw()
        fig1.canvas.flush_events()
        plt.savefig(os.path.join(saveDir, filename + '_fidelity_plot.png'))

        # Figure 3 plots IQ data for each of the runs on an IQ plane to help with fidelity measurements
        fig2 = plt.figure(2, figsize = (13, 8))
        plt.clf()
        plt.axhline(0)
        plt.axvline(0)
        plt.xlabel('I data (DAC a.u.)')
        plt.ylabel('Q data (DAC a.u.)')
        plt.title('I/Q distribution of measurements')
        for i in range(0, seq+1):
            plt.scatter(I_finals[i], Q_finals[i], label = 'Sequence {}'.format(i))
        plt.legend()
        fig2.canvas.draw()
        fig2.canvas.flush_events()
        plt.savefig(os.path.join(saveDir, filename+'_IQ_distribution.png'), dpi = 150)


    if exp_settings['debug']:
        userfuncs.SaveFull(saveDir, filename, ['power_data', 'phase_data', 'I_finals', 'Q_finals',
                            'I_full', 'Q_full','xaxis', 'full_IQ_data'],
                            locals(), expsettings=settings, instruments={})
        
        data = {'saveDir': saveDir, 'filename': filename, 'I_full': I_full, 'Q_full': Q_full, 
                'I_finals':I_finals_debug, 'Q_finals':Q_finals_debug, 'xaxis':xaxis, 'power_data':power_debug, 
                'phase_data': phase_debug, 'full_IQ_data':full_IQ_data}
    
    else:
        userfuncs.SaveFull(saveDir, filename, ['power_data', 'phase_data', 'I_finals', 'Q_finals'],
                                locals(), expsettings=settings, instruments={})
        
        data = {'saveDir': saveDir, 'filename': filename, 'I_full': I_full, 'Q_full': Q_full,
                'I_finals':I_finals, 'Q_finals':Q_finals, 'power_data':power_data, 'phase_data':phase_data }
    
    return data, power_data, phase_data, prog_RB
