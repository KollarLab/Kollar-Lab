import os
import time
import numpy as np

import zhinst.ziPython as ziPython
import zhinst.utils as ziUtils
import settingTools as sT

class HDAWG():
    apilevel = 6 # Sets the 'modernity' of the api used (bigger is more modern, 6 is current max version)
    #defaults_path='C:\\Users\\kollarlab\\AppData\\Roaming\\Zurich Instruments\\LabOne\\WebServer\\setting'
    #defaults_path=os.getcwd()
    defaults_path  = os.path.dirname(os.path.realpath(__file__))
    defaults_filename='default_settings.xml'

    def __init__(self, device):
        (self.daq,self.device,self.props) = ziUtils.create_api_session(device, self.apilevel) #connect to device specified by string
        ziUtils.disable_everything(self.daq,self.device) #disable all outputs of device
        self.load_default()

    #want method to load default settings (known state) to device
    def load_default(self):
        file = self.defaults_path + os.sep + self.defaults_filename
        ziUtils.load_settings(self.daq, self.device, file)

    #load settings from xml file
    def load_settings_xml(self,path,filename):
        file = path + os.sep + filename
        ziUtils.load_settings(self.daq,self.device, file)

    #save settings to xml file
    def save_settings_xml(self, path, filename):
        file = path + os.sep + filename
        ziUtils.save_settings(self.daq, self.device, file)

    #method to load arb waveform to device
    def load_waveform(self, index, wave1, wave2=None, markers=None):
        wave_awg = ziUtils.convert_awg_waveform(wave1, wave2, markers)
        path='/{:s}/awgs/0/waveform/waves/{:d}'.format(self.device, index)
        self.daq.setVector(path, wave_awg)

    #method to read waveform from device, requires it to be part of a 'playwave' command in the sequencer code
    def read_waveform(self,index, channels=1, markers_present=False):
        path='/{:s}/awgs/0/waveform/waves/{:d}'.format(self.device, index) #Path to wave data
        #prepare data to be read out, unclear what this does exactly but in combination with the poll command this will read arbitrarily long data
        checktime=0.1 #time to record data for (s)
        timeout=5 #timeout for poll command (ms)
        self.daq.getAsEvent(path) 
        pollReturn=self.daq.poll(checktime,timeout,flat=True)
        #poll command returns dictionary so we need to access the actual array data
        wave_awg=pollReturn[path][0]['vector']
        [ch1,ch2,markers]=ziUtils.parse_awg_waveform(wave_awg,channels,markers_present)
        return [ch1,ch2,markers]

    #load program to AWG, essentially copied from the example section for HDAWG
    def load_program(self, awg_program):
        awgModule = self.daq.awgModule()
        awgModule.set('device',self.device)
        awgModule.execute()

        awgModule.set('compiler/sourcestring',awg_program)
        # Note: when using an AWG program from a source file (and only then), the compiler needs to
        # be started explicitly with awgModule.set('compiler/start', 1)
        while awgModule.getInt('compiler/status') == -1:
            time.sleep(0.1)

        if awgModule.getInt('compiler/status') == 1:
            # compilation failed, raise an exception
            raise Exception(awgModule.getString('compiler/statusstring'))
            
        if awgModule.getInt('compiler/status') == 0:
            print("Compilation successful with no warnings, will upload the program to the instrument.")
        if awgModule.getInt('compiler/status') == 2:
            print("Compilation successful with warnings, will upload the program to the instrument.")
            print("Compiler warning: ", awgModule.getString('compiler/statusstring'))

        # Wait for the waveform upload to finish
        time.sleep(0.2)
        i = 0
        while (awgModule.getDouble('progress') < 1.0) and (awgModule.getInt('elf/status') != 1):
            #print("{} progress: {:.2f}".format(i, awgModule.getDouble('progress')))
            time.sleep(0.001)
            i += 1
        print("{} progress: {:.2f}".format(i, awgModule.getDouble('progress')))

        if awgModule.getInt('elf/status') == 0:
            print("Upload to the instrument successful.")
        if awgModule.getInt('elf/status') == 1:
            raise Exception("Upload to the instrument failed.")

    #Set amplitude of AWG
    def set_AWGamp(self,amps, channels):
        ranges = [0.2,0.4,0.6,0.8,1.0,2.0,3.0,4.0,5.0]
        for amp, channel in zip(amps, channels):
            for maxval in ranges:
                if amp <= maxval:
                    print('Channel: {}, range: {}'.format(channel, maxval))
                    self.daq.setDouble('/{}/sigouts/{}/range'.format(self.device, channel), maxval) #set the max range of output channel
                    scaled_amp = amp/maxval
                    AWGmod = int(np.floor(channel/2)) #get the AWG module number (for 8163, there are two built in modules)
                    AWGchannel = np.mod(channel, 2) #get channel number within module
                    self.daq.setDouble('/{}/awgs/{}/outputs/{}/amplitude'.format(self.device, AWGmod, AWGchannel), scaled_amp) #scale amplitude of AWG data to match desired voltage
                    while(self.daq.getInt('/{}/sigouts/{}/busy'.format(self.device, channel))>0): #wait until channel settings are completed 
                        time.sleep(0.1)
                    break
                else:
                    if amp > max(ranges):
                        raise Exception("Amplitude is too large, max amplitude is {}".format(max(ranges)))
                    continue

    #Functions to do settings stuff, there's got to be a better way but this is a start
    def enable_channels(self, channels):
        for channel in channels:
            self.daq.setDouble('/{}/sigouts/{}/on'.format(self.device,channel),1)
            while(self.daq.getInt('/{}/sigouts/{}/busy'.format(self.device, channel))>0): #wait until channel settings are completed 
                time.sleep(0.1)

    #Disable/ turn off channels 
    def disable_channels(self, channels={}):
        for channel in channels:
            self.daq.setDouble('/{}/sigouts/{}/on'.format(self.device,channel),0)
        if not channels:
            for channel in range(4): 
                self.daq.setDouble('/{}/sigouts/{}/on'.format(self.device,channel),0)
        
    #Enable Marker channels (by default AWG will use Triggers sequencer commands but they are less precise)
    def enable_markers(self,markers={}):
        for marker in markers:
            self.daq.setDouble('/{}/triggers/out/{}/source'.format(self.device,marker),marker+4)
        if not markers:
            for marker in range(4):
                self.daq.setDouble('/{}/triggers/out/{}/source'.format(self.device,marker),marker+4)

    #Use Triggers instead of markers for waveforms
    #TRIGGERS don't seem to work at the moment (can't get them in GUI...)
    def enable_triggers(self,markers={}):
        for marker in markers:
            self.daq.setDouble('/{}/triggers/out/{}/source'.format(self.device,marker),marker)
        if not markers:
            for marker in range(4):
                self.daq.setDouble('/{}/triggers/out/{}/source'.format(self.device,marker),marker)

    #Turn on AWG and run through program once
    def run(self,AWGcore=0):
        self.daq.setInt('/{}/awgs/{}/single'.format(self.device,AWGcore),1)

    #Run AWG continuously
    def run_loop(self,AWGcore=0):
        self.daq.setInt('/{}/awgs/{}/enable'.format(self.device,AWGcore),1)
    
    #disconnect device
    def done(self):
        self.daq.disconnectDevice(self.device)
