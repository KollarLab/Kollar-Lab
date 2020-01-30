import os
import time

import zhinst.ziPython as ziPython
import zhinst.utils as ziUtils

class HDAWG():
    apilevel = 6 # Sets the 'modernity' of the api used (bigger is more modern, 6 is current max version)
    #defaults_path='C:\\Users\\kollarlab\\AppData\\Roaming\\Zurich Instruments\\LabOne\\WebServer\\setting'
    defaults_path=os.getcwd()
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
    def load_settings(self,path,filename):
        file = path + os.sep + filename
        ziUtils.load_settings(self.daq,self.device, file)

    #save settings to xml file
    def save_settings(self, path, filename):
        file = path + os.sep + filename
        ziUtils.save_settings(self.daq, self.device, file)

    #method to load arb waveform to device
    def load_waveform(self, index, wave1, wave2=None, markers=None):
        wave_awg = ziUtils.convert_awg_waveform(wave1, wave2, markers)
        path='/{:s}/awgs/0/waveform/waves/{:d}'.format(self.device, index)
        self.daq.setVector(path, wave_awg)

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
            print("{} progress: {:.2f}".format(i, awgModule.getDouble('progress')))
            time.sleep(1.0)
            i += 1
        print("{} progress: {:.2f}".format(i, awgModule.getDouble('progress')))

        if awgModule.getInt('elf/status') == 0:
            print("Upload to the instrument successful.")
        if awgModule.getInt('elf/status') == 1:
            raise Exception("Upload to the instrument failed.")

    #Functions to do settings stuff, there's got to be a better way but this is a start
    def enable_channels(self, channels):
        for channel in channels:
            self.daq.setDouble('/{}/sigouts/{}/on'.format(self.device,channel),1)

    def disable_channels(self, channels={}):
        for channel in channels:
            self.daq.setDouble('/{}/sigouts/{}/on'.format(self.device,channel),0)
        if not channels:
            for channel in range(4): 
                self.daq.setDouble('/{}/sigouts/{}/on'.format(self.device,channel),0)
        
    #Turn on AWG
    def AWG_run(self):
        self.daq.setInt('/dev8163/awgs/0/enable',1)
    
    #disconnect device
    def done(self):
        self.daq.disconnectDevice(self.device)
