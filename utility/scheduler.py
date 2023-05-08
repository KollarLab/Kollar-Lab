import numpy as np
import matplotlib.pyplot as plt

def gaussian(sample_rate, amp, sigma, num_sigma):
    samples = int(sigma*num_sigma*sample_rate)
    t = np.linspace(0,sigma*num_sigma, samples)
    t0 = num_sigma*sigma/2
    wave = np.exp(-(t-t0)**2/(2*sigma**2))
    offset = wave[0]
    init_amp = max(wave)
    wave = (wave-offset)/(init_amp-offset)
    return wave*amp

def gaussian_square(sample_rate, amp, length, num_sigma, ramp_sigma):
    ramp = gaussian(sample_rate, amp, ramp_sigma, num_sigma)
    square = np.ones(int(length*sample_rate))*amp
    ramp_up, ramp_down = np.split(ramp, 2)
    final = np.concatenate((ramp_up, square, ramp_down))
    return final
       
class AnalogChannel():
    def __init__(self, chanID, sample_rate, samples, name='waveform'):
        self.ID = chanID
        self.sample_rate = sample_rate
        self.samples = samples
        self.time_array = samples/sample_rate
        self.wave_array = np.zeros(samples)
        self.name = name
    
    def add_pulse(self, type, position, amplitude=1, **keyargs):
        if type=='gaussian':
            sigma = keyargs['sigma']
            num_sigma = keyargs['num_sigma']
            pulse = gaussian(self.sample_rate, amplitude, sigma, num_sigma)
            pos_index = int(position*self.sample_rate)
            self.wave_array[pos_index:pos_index+len(pulse)] = pulse
        if type=='gaussian_square':
            ramp_sigma = keyargs['ramp_sigma']
            num_sigma  = keyargs['num_sigma']
            length = keyargs['length']
            pulse = gaussian_square(self.sample_rate, amplitude, length, num_sigma, ramp_sigma)
            pos_index = int(position*self.sample_rate)
            self.wave_array[pos_index:pos_index+len(pulse)] = pulse
        if type=='custom':
            pos_index = int(position*self.sample_rate)
            wave_data = keyargs['wave']
            self.wave_array[pos_index:pos_index+len(wave_data)] = amplitude*wave_data
            
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

class scheduler():
    def __init__(self, total_time, num_dig_channels=0, num_analog_channels=0, sample_rate=2.4e9):
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
        for chan in range(num_dig_channels):
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

    def compile_schedule(self, AWG_type='HDAWG', analog_list=[], digital_list=[]):
        if AWG_type=='HDAWG':
            if digital_list:
                mark1 = self.digital_channels[digital_list[0]]
                mark2 = self.digital_channels[digital_list[1]]
                mark1.compile_channel()
                mark2.compile_channel()
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

    def plot_waveforms(self):
        analog_traces = len(self.analog_channels)
        digital_traces = len(self.digital_channels)
        num_traces = analog_traces+digital_traces
        fig = plt.figure(111)
        plt.clf()
        for i, key in enumerate(self.analog_channels.keys()):
            if i==0:
                ax1 = plt.subplot(num_traces, 1, i+1)
                ax = ax1
            else:
                ax = plt.subplot(num_traces, 1, i+1, sharex=ax1)
            active_channel = self.analog_channels[key]
            wave_data = active_channel.wave_array
            time_data = self.time_array*1e6
            ax.plot(time_data, wave_data)
            ax.set_title(key)
            ax.set_xlabel('Time (us)')
            ax.set_ylabel('Amp')

        for i, key in enumerate(self.digital_channels.keys()):
            ax = plt.subplot(num_traces, 1, i+analog_traces+1, sharex=ax1)
            active_channel = self.digital_channels[key]
            active_channel.compile_channel()
            marker_data = active_channel.marker_array
            time_data = self.time_array*1e6
            ax.plot(time_data, marker_data)
            ax.set_title(key)
            ax.set_xlabel('Time (us)')
            ax.set_ylabel('Amp')

#        fig.canvas.draw()
#        fig.canvas.flush_events()
    
    def reset(self):
        for chan in self.analog_channels.values():
            chan.reset()
        for chan in self.digital_channels.values():
            chan.reset()
