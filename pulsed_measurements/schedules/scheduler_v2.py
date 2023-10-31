import numpy as np
import matplotlib.pyplot as plt

class pulse():
    def __init__(self, position, wave):
        self.position = position
        self.wave = wave

class flat_top(pulse):
    _defaults = {
            'sigma'       : 10e-9,
            'num_sigma'   : 4,
            'hold_time'   : 0e-9,
            'amp'         : 1, 
            'sample_rate' : 2.4e9,
            'angle'       : 0
                }
    def __init__(self, position, settings):
        self._settings = {**self._defaults, **settings}
        self.set_settings()
        self.position = position
        self.wave = self.compile_pulse()
        
    def set_settings(self):
        for parameter in self._settings.keys():
            if parameter not in self._defaults.keys():
                print('Error or typo, {} is not a valid parameter for this class'
                      .format(parameter))
            else:
                self.__setattr__(parameter, self._settings[parameter])
    
    @property
    def settings(self):
        full_settings = {}
        for param in self._defaults.keys():
            full_settings[param] = getattr(self, param)
        return full_settings

    def compile_pulse(self):
        
        amp         = self.amp
        sigma       = self.sigma
        length      = self.hold_time
        num_sigma   = self.num_sigma
        sample_rate = self.sample_rate
        
        # Create gaussian ramp
        samples = int(sigma*num_sigma*sample_rate)
        t = np.linspace(0,sigma*num_sigma, samples)
        t0 = num_sigma*sigma/2
        wave = np.exp(-(t-t0)**2/(2*sigma**2))
        offset = wave[0]
        init_amp = max(wave)
        ramp = amp*(wave-offset)/(init_amp-offset)
        #Hold max amplitude
        square = np.ones(int(length*sample_rate))*amp
        ramp_up, ramp_down = np.split(ramp, 2)
        final = np.concatenate((ramp_up, square, ramp_down))

        self.wave = final
        self.I = np.cos(self.angle)*final
        self.Q = np.sin(self.angle)*final

class marker():
    _defaults = {
        'position' : 50e-6,
        'length'   : 5e-6,
        'buffer'   : 100e-9
    }
    def __init__(self, settings):
        self._settings = {**self._defaults, **settings}
        self.set_settings()

    def set_settings(self):
        for parameter in self._settings.keys():
            if parameter not in self._defaults.keys():
                print('Error or typo, {} is not a valid parameter for this class'
                      .format(parameter))
            else:
                self.__setattr__(parameter, self._settings[parameter])
    
    @property
    def settings(self):
        full_settings = {}
        for param in self._defaults.keys():
            full_settings[param] = getattr(self, param)
        return full_settings

class AnalogChannel():
    def __init__(self, chanID, sample_rate, samples, name='waveform'):
        self.ID = chanID
        self.sample_rate = sample_rate
        self.samples = samples
        self.time_array = samples/sample_rate
        self.wave_array = np.zeros(samples)
        self.name = name
    
    def add_pulse(self, pulse):
        pos_index = int(pulse.position*self.sample_rate)
        #print('index:{}, len array:{}, array:{}'.format(pos_index, len(pulse.wave), len(self.wave_array)))
        self.wave_array[pos_index:pos_index+len(pulse.wave)] = pulse.wave

    def reset(self):
        self.wave_array = np.zeros(self.samples)
            
class DigitalChannel():
    def __init__(self, chanID, sample_rate, samples, name='marker', HW_offset_on=0, HW_offset_off=0, polarity='Pos'):
        self.ID = chanID
        self.polarity = polarity
        self.sample_rate = sample_rate
        self.samples = samples
        if polarity=='Pos':
            self.marker_array = np.zeros(samples)
        elif polarity=='Neg':
            self.marker_array = np.ones(samples)*chanID
        else:
            print("Acceptable polarity is 'Neg' or 'Pos'")
        self.polarity = polarity
        self.HW_offset_on  = HW_offset_on
        self.HW_offset_off = HW_offset_off
        self.name = name
        self.window_list = []

    def add_window(self, start, stop):
        self.window_list.append([start, stop])
    
    def compile_channel(self):
        for window in self.window_list:
            start_ind  = int(window[0]*self.sample_rate)
            stop_ind   = int(window[1]*self.sample_rate)
            offset_on  = int(self.HW_offset_on*self.sample_rate)
            offset_off = int(self.HW_offset_off*self.sample_rate)
            if self.polarity=='Pos':
                self.marker_array[start_ind-offset_on:stop_ind-offset_off] = self.ID
            if self.polarity=='Neg':
                self.marker_array[start_ind-offset_on:stop_ind-offset_off] = 0
    
    def reset(self, polarity=None):
        if not polarity:
            polarity = self.polarity
        if polarity=='Pos':
            self.marker_array = np.zeros(self.samples)
        elif polarity=='Neg':
            self.marker_array = np.ones(self.samples)*self.chanID
        self.window_list = []

class base_schedule():
    def __init__(self, total_time, num_analog_channels=0, num_digital_channels=0, sample_rate=2.4e9):
        ## Choose the right HDAWG file for the 2 vs 4 channel mode
        ## Create all the channels here?
        self.sample_rate = sample_rate
        self.total_time = total_time
        samples = int(total_time*sample_rate)
        if samples%32!=0:
            print('Warning, waveform lengths not a multiple of 32, padding to match AWG reqs')
            samples = np.ceil(samples/32)*32
        self.samples = int(samples)
        self.time_array = np.linspace(0, total_time, self.samples)
        self.analog_channels = {}
        self.digital_channels = {}
        for chan in range(num_analog_channels):
            self.add_analog_channel(chan+1)
        for chan in range(num_digital_channels):
            self.add_digital_channel(chan+1)

    def add_analog_channel(self, HW_ID, name='waveform{}'):
        sample_rate = self.sample_rate
        samples = self.samples
        name = name.format(HW_ID)
        self.analog_channels[name] = AnalogChannel(HW_ID, sample_rate, samples, name)
    
    def add_digital_channel(self, HW_ID, name='marker{}', polarity='Pos', HW_offset_on=0, HW_offset_off=0):
        sample_rate = self.sample_rate
        samples = self.samples
        name = name.format(HW_ID)
        self.digital_channels[name] = DigitalChannel(HW_ID, sample_rate, samples, name, HW_offset_on, HW_offset_off, polarity)

    def add_pulse(self, pulse, channel):
        self.pulse_list.append([pulse, channel])
    
    def add_marker(self, pulse, buffer, channel):
        position = pulse.position
        length   = len(pulse.wave)/self.sample_rate+2*buffer
        mark_channel = self.digital_channels[channel]
        mark_channel.add_window(position-buffer, position-buffer+length)

    def compute_schedule(self):
        ## Loop through all the channels and add their pulses to the arrays
        for [pulse, channel] in self.pulse_list:
            self.analog_channels[channel].add_pulse(pulse)
        for channel in self.digital_channels.keys():
            ch = self.digital_channels[channel]
            ch.compile_channel()
    
    def compile_schedule(self, AWG_type='HDAWG', analog_list=[], digital_list=[]):
        if AWG_type=='HDAWG':
            if digital_list:
                mark1 = self.digital_channels[digital_list[0]]
                mark2 = self.digital_channels[digital_list[1]]
                final_markers = mark1.marker_array + mark2.marker_array
            else:
                print('No markers provided, using blank array')
                final_markers = np.zeros(self.samples)
            if analog_list:
                analog1 = self.analog_channels[analog_list[0]].wave_array
                analog2 = self.analog_channels[analog_list[1]].wave_array
            else:
                print('No analog wf provided, using blank array')
                analog1 = np.zeros(self.samples)
                analog2 = np.zeros(self.samples)
            return [analog1, analog2, final_markers]
        else:
            print('Not implemented')

    def add_measurement(self, mpulse, channel):
        position = mpulse.position
        length   = mpulse.length
        cav_channel = self.digital_channels[channel]
        cav_channel.add_window(position, position+length)

    def plot_waveforms(self):
        analog_traces = len(self.analog_channels)
        digital_traces = len(self.digital_channels)
        num_traces = analog_traces+digital_traces
        fig = plt.figure(111, figsize=(10,8))
        plt.clf()
        ax_a = plt.subplot(211)
        ax_a.set_title('Analog Traces')
        ax_a.set_xlabel('Time (us)')
        ax_a.set_ylabel('Amp')
        ax_b = plt.subplot(212, sharex=ax_a)
        ax_b.set_title('Digital Traces')
        ax_b.set_xlabel('Time (us)')
        ax_b.set_ylabel('Amp')
        for a,b in zip(self.analog_channels.keys(), self.digital_channels.keys()):
            active_channel = self.analog_channels[a]
            wave_data = active_channel.wave_array
            time_data = self.time_array*1e6
            ax_a.plot(time_data, wave_data,label=a)
            active_channel = self.digital_channels[b]
            active_channel.compile_channel()
            wave_data = active_channel.marker_array
            time_data = self.time_array*1e6
            ax_b.plot(time_data, wave_data, label=b)
        ax_a.legend()
        ax_b.legend()
    
    def reset(self):
        self.pulse_list = []
        for chan in self.analog_channels.values():
            chan.reset()
        for chan in self.digital_channels.values():
            chan.reset()

