import os
import time
import numpy as np

import zhinst.ziPython as ziPython
import zhinst.utils as ziUtils
import settingTools as sT
from enum import Enum
from bidict import bidict

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
        self.nodepaths = self.fill_paths(self.device)
        self.Channels = []
        self.Triggers = []
        self.AWGs     = []
        for i in range(4):
            self.Channels.append(HDAWGchannel(self.daq,i))
        for i in range(2):
            self.Triggers.append(HDAWGtrigger(self.daq,i))
            self.AWGs.append(HDAWGawg(self.daq,i))

    def fill_paths(self,device):
        nodes={}
        nodes['channelgrouping']='/{}/system/awg/channelgrouping'.format(device)
        nodes['samplerate']='/{}/awgs/0/time'.format(device)
        nodes['referenceClock']='/{}/system/clocks/referenceclock/source'.format(device)
        return nodes

    _groupInt=bidict({
        '2x2':0,
        '1x4':1
    })
    _refInt=bidict({
        'Internal':0,
        'External':1
    })

    @property
    def channelgrouping(self):
        node = self.nodepaths['channelgrouping']
        val = self.daq.getInt(node)
        return self._groupInt.inverse[val]
    @channelgrouping.setter
    def channelgrouping(self,val):
        node = self.nodepaths['channelgrouping']
        if val not in self._groupInt.keys():
            print('Invalid channel grouping, must be {}'.format(list(self._groupInt)))
        else:
            self.daq.setInt(node,self._groupInt[val])
    @property
    def referenceClock(self):
        node = self.nodepaths['referenceClock']
        val = self.daq.getInt(node)
        return self._refInt.inverse[val]
    @referenceClock.setter
    def referenceClock(self,val):
        node = self.nodepaths['referenceClock']
        if val not in self._refInt:
            print('Invalid reference clock parameter, must be {}'.format(list(self._refInt)))
        else:
            self.daq.setInt(node,self._refInt[val])

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

    #disconnect device
    def done(self):
        self.daq.disconnectDevice(self.device)

class HDAWGawg:
    def __init__(self, daq, AWGcore, samplerate='2.4Ghz',device='dev8163', channelgrouping='1x4'):
        self.device          = device
        self.daq             = daq
        self.index           = AWGcore
        self.nodepaths       = self.fill_paths()
        self.channelgrouping = channelgrouping
        self.Triggers        = []
        for i in range (2):
            self.Triggers.append(HDAWGtrigger(self.daq,i,AWGcore=self.index))
        self.samplerate      = samplerate

    def fill_paths(self):
        nodes = {}
        nodes['waves']='/{}/awgs/{}/waveform/waves/index'.format(self.device, self.index)
        nodes['single']='/{}/awgs/{}/single'.format(self.device,self.index)
        nodes['enable']='/{}/awgs/{}/enable'.format(self.device,self.index)
        return nodes
    
    def getSettings(self):
        settings = {}
        settings['samplerate'] = self.samplerate
        settings['Triggers'] = {}
        for i in range(2):
            temp = 'Trigger{}'.format(i)
            settings['Triggers'][temp]=self.Triggers[i].getSettings
        return settings
    
    def setSettings(self, settings):
        self.samplerate = settings['samplerate']
        if 'Triggers' in settings.keys():
            for i in range(len(settings['Triggers']))


    
    _sampleInt=bidict({
        '2.4GHz'   : 0,
        '1.2GHz'   : 1,
        '600MHz'   : 2,
        '300MHz'   : 3,
        '150MHz'   : 4,
        '75MHz'    : 5,
        '37.5MHz'  : 6,
        '18.75MHz' : 7,
        '9.38MHz'  : 8,
        '4.69MHz'  : 9,
        '2.34MHz'  : 10,
        '1.17MHz'  : 11,
        '585.94kHz': 12,
        '292.97kHz': 13
    })

    @property
    def samplerate(self):
        node = self.nodepaths['samplerate']
        val = self.daq.getInt(node)
        return self._sampleInt.inverse[val]
    @samplerate.setter
    def samplerate(self,val):
        node = self.nodepaths['samplerate']
        if val not in self._sampleInt.keys():
            print('Invalid samplerate, acceptable values are: {}'.format(list(self._sampleInt)))
        else:
            self.daq.setInt(node,self._sampleInt[val])


    #method to load arb waveform to device
    def load_waveform(self, index, wave1, wave2=None, markers=None):
        wave_awg = ziUtils.convert_awg_waveform(wave1, wave2, markers)
        node = self.nodepaths['waves'].replace('index',index)
        #path='/{:s}/awgs/0/waveform/waves/{:d}'.format(self.device, index)
        self.daq.setVector(node, wave_awg)

    #method to read waveform from device, requires it to be part of a 'playwave' command in the sequencer code
    def read_waveform(self,index, channels=1, markers_present=False):
        node = self.nodepaths['waves'].replace('index',index)
        #path='/{:s}/awgs/0/waveform/waves/{:d}'.format(self.device, index) #Path to wave data
        #prepare data to be read out, unclear what this does exactly but in combination with the poll command this will read arbitrarily long data
        checktime=0.1 #time to record data for (s)
        timeout=5 #timeout for poll command (ms)
        self.daq.getAsEvent(node) 
        pollReturn=self.daq.poll(checktime,timeout,flat=True)
        #poll command returns dictionary so we need to access the actual array data
        wave_awg=pollReturn[node][0]['vector']
        [ch1,ch2,markers]=ziUtils.parse_awg_waveform(wave_awg,channels,markers_present)
        return [ch1,ch2,markers]

    #load program to AWG, essentially copied from the example section for HDAWG
    def load_program(self, awg_program):
        awgModule = self.daq.awgModule()
        awgModule.set('device',self.device)
        awgModule.set('index',self.index)
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

    #Turn on AWG and run through program once
    def run(self,AWGcore=0):
        node = self.nodepaths['single']
        #self.daq.setInt('/{}/awgs/{}/single'.format(self.device,AWGcore),1)
        self.daq.setInt(node,1)

    #Run AWG continuously
    def run_loop(self,AWGcore=0):
        node = self.nodepaths['enable']
        #self.daq.setInt('/{}/awgs/{}/enable'.format(self.device,AWGcore),1)
        self.daq.setInt(node,1)

class HDAWGtrigger():
    def __init__(self, daq, triggerID, device='dev8163', AWGcore=0, slope='rising', channel='Tin1'):
        self.device  = device
        self.daq     = daq
        self.ID      = triggerID
        self.AWGcore = AWGcore
        self.nodepaths = self.fill_paths(device, AWGcore, triggerID)
        self.slope   = slope
        self.channel = channel
    
    def fill_paths(self, device, AWGcore, triggerID):
        nodes = {}
        nodes['slope']='/{}/awgs/{}/auxtriggers/{}/slope'.format(device, AWGcore, triggerID)
        nodes['channel']='/{}/awgs/{}/auxtriggers/{}/channel'.format(device, AWGcore, triggerID)
        return nodes

    def configureTrigger(self,slope, channel):
        self.slope=slope
        self.channel=channel

    def getSettings(self):
        settings={}
        for key in self.nodepaths.keys():
            settings[key]=getattr(self,key)
        return settings

    def setSettings(self,settings):
        for key in settings.keys():
            setattr(self,key,settings[key])

    _slopeInt = bidict({
        'level':0,
        'rising':1,
        'falling':2,
        'both':3
    })
    _channelInt = bidict({
        'Trigger in 1':0,
        'Trigger in 2':1,
        'Trigger in 3':2,
        'Trigger in 4':3,
        'Trigger out 1':4,
        'Trigger out 2':5,
        'Trigger out 3':6,
        'Trigger out 4':7,
    })    

    @property
    def slope(self):
        node = self.nodepaths['slope']
        val = self.daq.getInt(node)
        return self._slopeInt.inverse[val]
    @slope.setter
    def slope(self, val):
        node = self.nodepaths['slope']
        if val not in self._slopeInt.keys():
            print('Invalid trigger slope, acceptable values are: {}'.format(list(self._slopeInt)))
        else:
            self.daq.setInt(node,self._slopeInt[val])

    @property
    def channel(self):
        node = self.nodepaths['channel']
        val = self.daq.getInt(node)
        return self._channelInt.inverse[val]
    @channel.setter
    def channel(self, val):
        node = self.nodepaths['channel']
        if val not in self._channelInt.keys():
            print('Invalid trigger channel, acceptable values are: {}'.format(list(self._channelInt)))
        else:
            self.daq.setInt(node,self._channelInt[val])

class HDAWGchannel():
    def __init__(self, daq, channelID, device='dev8163', amp = 1.0, fullscale = 1.0, AWGamp = 1.0, offset = 0.0, delay = 0.0, markers = 'On'):
        self.device    = device
        self.daq       = daq
        self.ID        = channelID
        self.AWGcore   = np.floor(channelID/2)
        self.nodepaths = self.fill_paths(device, channelID)
        self.status    = 'Off'
        self.fullscale = fullscale
        self.AWGamp    = AWGamp
        self.amp       = amp
        self.offset    = offset
        self.delay     = delay
        self.markers   = markers

    def fill_paths(self, device, channelID):
        nodes = {}
        nodes['fullscale'] = '/{}/sigouts/{}/range'.format(device,channelID)
        nodes['offset'] = '/{}/sigouts/{}/offset'.format(device,channelID)
        nodes['delay'] = '/{}/sigouts/{}/delay'.format(device,channelID)
        nodes['AWGamp'] = '/{}/awgs/{}/outputs/{}/amplitude'.format(device, self.AWGcore, channelID)
        nodes['markers'] = '/{}/triggers/out/{}/source'.format(device, channelID)
        nodes['status'] = '/{}/sigouts/{}/on'.format(device, channelID)
        return nodes

    def configureChannel(self, amp=1.0, fullscale=1.0, AWGamp=1.0, delay=0.0, offset=0.0, markers = 'Off'):
        if amp < 1.0 or amp > 1.0:
            self.amp = amp
        else:
            self.fullscale = fullscale
            self.AWGamp    = AWGamp
            print('Manual range configuration')

        self.delay   = delay
        self.offset  = offset
        self.markers = markers
        self.status  = 'On'

    def getSettings(self):
        settings={}
        for key in self.nodepaths.keys():
            settings[key] = getattr(self,key)
        return settings

    def setSettings(self, settings):
        for key in settings:
            setattr(self,key,settings[key])

    @property
    def status(self):
        node = self.nodepaths['status']
        val = self.daq.getInt(node)
        if val == 0:
            return 'Off'
        else:
            return 'On'
    @status.setter
    def status(self,val):
        node = self.nodepaths['status']
        if val == 'On':
            self.daq.setInt(node,1)
        elif val == 'Off':
            self.daq.setInt(node,0)
        else:
            raise ValueError('Status must be "On" or "Off"')

    @property
    def amp(self):
        return self.fullscale*self.AWGamp
    @amp.setter
    def amp(self,val):
        print('Automatically adjusting fullscale and AWGamp values')
        ranges = [0.2,0.4,0.6,0.8,1.0,2.0,3.0,4.0,5.0]
        for x in ranges:
            if val > x:
                #print('{} is greater than {}'.format(val,x))
                continue
            else:
                self.fullscale = x
                self.AWGamp    = val/x
                break

    @property
    def markers(self):
        node      = self.nodepaths['markers']
        markindex = self.ID+4
        status    = self.daq.getDouble(node)
        if status == markindex:
            return 'On'
        else:
            return 'Off'
    @markers.setter
    def markers(self,val):
        node      = self.nodepaths['markers']
        markindex = self.ID+4
        if val == 'On':
            print('Enabling markers')
            self.daq.setDouble(node,markindex)
        else:
            print('Disabling markers')
            self.daq.setDouble(node,self.ID)
     
    @property
    def fullscale(self):
        node = self.nodepaths['fullscale']
        val  = self.daq.getDouble(node)
        return val
    @fullscale.setter
    def fullscale(self,val):
        ranges = [0.2,0.4,0.6,0.8,1.0,2.0,3.0,4.0,5.0]
        node   = self.nodepaths['fullscale']
        if val not in ranges:
            raise ValueError('Fullscale values must one of {}'.format(ranges))
        else:
            self.daq.setDouble(node,val)
    
    @property
    def AWGamp(self):
        node = self.nodepaths['AWGamp']
        val  = self.daq.getDouble(node)
        return val
    @AWGamp.setter
    def AWGamp(self,val):
        node = self.nodepaths['AWGamp']
        if(val<=0.0 or val>1.0):
            raise ValueError('AWG amp out of range (0.,1]')
        else:
            self.daq.setDouble(node, val)

    @property
    def offset(self):
        node = self.nodepaths['offset']
        val  = self.daq.getDouble(node)
        return val
    @offset.setter
    def offset(self,val):
        node = self.nodepaths['offset']
        self.daq.setDouble(node,val)

    @property
    def delay(self):
        node = self.nodepaths['delay']
        val  = self.daq.getDouble(node)
        return val
    @delay.setter
    def delay(self,val):
        node = self.nodepaths['delay']
        self.daq.setDouble(node,val)
    
    #Set amplitude of AWG
    # def set_AWGamp(self,amps, channels):
    #     ranges = [0.2,0.4,0.6,0.8,1.0,2.0,3.0,4.0,5.0]
    #     for amp, channel in zip(amps, channels):
    #         for maxval in ranges:
    #             if amp <= maxval:
    #                 print('Channel: {}, range: {}'.format(channel, maxval))
    #                 self.daq.setDouble('/{}/sigouts/{}/range'.format(self.device, channel), maxval) #set the max range of output channel
    #                 scaled_amp = amp/maxval
    #                 AWGmod = int(np.floor(channel/2)) #get the AWG module number (for 8163, there are two built in modules)
    #                 AWGchannel = np.mod(channel, 2) #get channel number within module
    #                 self.daq.setDouble('/{}/awgs/{}/outputs/{}/amplitude'.format(self.device, AWGmod, AWGchannel), scaled_amp) #scale amplitude of AWG data to match desired voltage
    #                 while(self.daq.getInt('/{}/sigouts/{}/busy'.format(self.device, channel))>0): #wait until channel settings are completed 
    #                     time.sleep(0.1)
    #                 break
    #             else:
    #                 if amp > max(ranges):
    #                     raise Exception("Amplitude is too large, max amplitude is {}".format(max(ranges)))
    #                 continue

    # #Functions to do settings stuff, there's got to be a better way but this is a start
    # def enable_channels(self, channels):
    #     for channel in channels:
    #         self.daq.setDouble('/{}/sigouts/{}/on'.format(self.device,channel),1)
    #         while(self.daq.getInt('/{}/sigouts/{}/busy'.format(self.device, channel))>0): #wait until channel settings are completed 
    #             time.sleep(0.1)

    # #Disable/ turn off channels 
    # def disable_channels(self, channels={}):
    #     for channel in channels:
    #         self.daq.setDouble('/{}/sigouts/{}/on'.format(self.device,channel),0)
    #     if not channels:
    #         for channel in range(4): 
    #             self.daq.setDouble('/{}/sigouts/{}/on'.format(self.device,channel),0)
        
    # #Enable Marker channels (by default AWG will use Triggers sequencer commands but they are less precise)
    # def enable_markers(self,markers={}):
    #     for marker in markers:
    #         self.daq.setDouble('/{}/triggers/out/{}/source'.format(self.device,marker),marker+4)
    #     if not markers:
    #         for marker in range(4):
    #             self.daq.setDouble('/{}/triggers/out/{}/source'.format(self.device,marker),marker+4)

    # #Use Triggers instead of markers for waveforms
    # #TRIGGERS don't seem to work at the moment (can't get them in GUI...)
    # def enable_triggers(self,markers={}):
    #     for marker in markers:
    #         self.daq.setdouble('/{}/triggers/out/{}/source'.format(self.device,marker),marker)
    #     if not markers:
    #         for marker in range(4):
    #             self.daq.setDouble('/{}/triggers/out/{}/source'.format(self.device,marker),marker)

