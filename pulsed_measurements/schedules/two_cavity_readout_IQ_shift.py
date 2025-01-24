from .scheduler_v2 import base_schedule, pulse
import os
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

class two_cavity_readout_SSB(base_schedule):
    '''
    This is a simple program generating a qubit pulse (all in I atm) and a measurement tone.
    It currently has 4 sweepable parameters: hold time, tau, amplitude, and angle which can be accessed 
    as: schedule.tau = val. This should allow us to interface this with the sweeper class and treat it 
    as a meta instruments with SCPI like properties. It currently takes in a settings dictionary 
    for the qubit pulse (defining the main qubit params: sigma, hold time etc.) and another for 
    the measurement tone as well as a total schedule time (the time for which we should generate 
    traces). 
    It's currently configured to take in an AWG core instead of the full HDAWG object as that's 
    what we interface with when uploading stuff. 
    The idea is that if you want a new pulse schedule you can quickly define it by modifying the 
    "update" condition (here labelled as the "define schedule" function) as well as the initial 
    channel assignments. 
    '''
    def __init__(self, total_time, qpulse, mpulse1, mpulse2, SSB_freq=0, awg='sim', awg2='sim'):
        '''
        My current idea here is that we define the various channels used for the schedule and 
        store whatever extra variables we want. This will initialize a blank schedule object and
        populate whatever channels we want.
        '''
        super().__init__(total_time, sample_rate=2.4e9)
        self.add_analog_channel(1,'Qubit_I')
        self.add_analog_channel(2,'Qubit_Q')
        self.add_analog_channel(3,'Blank')
        self.add_analog_channel(4,'Blank2')
        self.add_digital_channel(1, 'Qubit_enable')
        self.add_digital_channel(2, 'Readout')
        self.add_digital_channel(1, 'Boost')
        self.add_digital_channel(2, 'Blank')
        
        self._tau   = 0

        self.awg_core = awg
        self.awg_core2 = awg2

        self.SSB_freq = SSB_freq
        self.qpulse = qpulse
        self.mpulse1 = mpulse1
        self.mpulse2 = mpulse2
        self.qpulse.compile_pulse()

        self.buffer = 100e-9
        self.delay = 150e-9

        if awg=='sim':
            print('Operating schedule in simulation mode, no waveforms will be uploaded')
            self.sim = True
        else:
            #Finding a standard way to point to the placeholder file where the HDAWG defaults are 
            directory = os.path.dirname(os.path.realpath(__file__))
            self.placeholder = os.path.join(str(Path(directory).parents[0]),'HDAWG_sequencer_codes','hdawg_placeholder_4channels.cpp')
            self.progFile = open(self.placeholder,'r')
            rawprog  = self.progFile.read()
            loadprog = rawprog
            loadprog = loadprog.replace('_samples_', str(self.samples))
            self.awg_core.load_program(loadprog)
            self.sim = False
    
    def define_schedule(self):
        '''
        This is the core of the function where we define the pulse list. It 
        is hopefully written in a simple way where you can add pulses/ measurements
        to whichever channel is desired. This also recompiles the schedule and 
        uploads it to the AWG (should be called everytime a property is modified so
        the changes actually happen in the code)
        '''
        self.compute_position()
        self.reset()
        qpulse = self.qpulse
        pulse_I = pulse(qpulse.position, qpulse.I)
        pulse_Q = pulse(qpulse.position, qpulse.Q)
        self.add_pulse(pulse_I, 'Qubit_I')
        self.add_pulse(pulse_Q, 'Qubit_Q')
        self.add_marker(self.qpulse, buffer=self.buffer, channel='Qubit_enable')
        self.add_measurement(self.mpulse1, 'Readout')
        self.add_measurement(self.mpulse2, 'Boost')
        self.compute_schedule()
        #self.plot_waveforms()
        if not self.sim:
            self.upload_schedule()

    def upload_schedule(self):
        ''' 
        Simple AWG uploader. I was initially considering adding it to the scheduler object
        but this allows for flexibility in defining which channels end up being used/ where
        things get uploaded
        '''
        self.awg_core.stop()
        [ch1, ch2, marker] = self.compile_schedule(analog_list=['Qubit_I', 'Qubit_Q'], digital_list=['Qubit_enable', 'Readout'])
        [ch3, ch4, marker2] = self.compile_schedule(analog_list=['Blank', 'Blank2'], digital_list=['Boost', 'Blank'])
        t_ax = np.linspace(0,len(ch1)-1, len(ch1))/2.4e9
        w_SSB = self.SSB_freq*2*np.pi
        self.ch1_up = ch1*np.cos(w_SSB*t_ax)+ch2*np.sin(w_SSB*t_ax)
        self.ch2_up = ch2*np.cos(w_SSB*t_ax)-ch1*np.sin(w_SSB*t_ax)    
        self.awg_core.load_waveform('0', self.ch1_up, self.ch2_up, marker)
        self.awg_core2.load_waveform('0', ch3, ch4, marker2)
        self.awg_core.run_loop()

    def compute_position(self):
        '''
        This is a QOL function to automatically compute the start position of the pulse 
        depending on the hold time, tau and qubit sigma. 
        '''
        qubit_time = len(self.qpulse.wave)/self.sample_rate
        delay = self._tau+self.delay
        self.qpulse.position = self.mpulse1.position-qubit_time-delay

    ''' 
    Here we define all the properties that are accessible to the user. For this simple program
    there are 4: tau (the delay between the end of the qubit pulse and the measurement tone), 
    hold_time (the time the qubit pulse stays high), amp (the amplitude of the pulse), and angle (
    angle of the pulse in IQ plane)
    '''
    @property
    def tau(self):
        return self._tau
    @tau.setter
    def tau(self, val):
        self._tau = val
        self.define_schedule()

    @property
    def hold_time(self):
        return self.qpulse.hold_time
    @hold_time.setter
    def hold_time(self, val):
        self.qpulse.hold_time = val
        self.qpulse.compile_pulse()
        self.define_schedule()
    
    @property
    def amp(self):
        return self.qpulse.amp
    @amp.setter
    def amp(self, val):
        self.qpulse.amp = val
        self.qpulse.compile_pulse()
        self.define_schedule()
    
    @property
    def angle(self):
        return self.qpulse.angle
    @angle.setter
    def angle(self, val):
        self.qpulse.angle = val
        self.qpulse.compile_pulse()
        self.define_schedule()

    @property
    def settings(self):
        ''' 
        Get a list of all the settings so we can store it in our final pickle file 
        '''
        full_settings = {}
        #Get a list of all the user defined properties to quickly add them to a dict
        user_props = [k for k, v in vars(type(self)).items() if isinstance(v, property)]
        for prop in user_props:
            if prop=='settings':
                continue
            else:
                full_settings[prop] = getattr(self,prop)
        full_settings['qpulse'] = self.qpulse.settings 
        full_settings['mpulse1'] = self.mpulse1.settings
        full_settings['mpulse2'] = self.mpulse2.settings
        return full_settings