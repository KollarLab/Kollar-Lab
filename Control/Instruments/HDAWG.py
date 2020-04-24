import os
import time
import numpy as np

import zhinst.ziPython as ziPython
import zhinst.utils as ziUtils
from . import settingTools as sT
#import settingTools as sT
from bidict import bidict

class HDAWG():
    '''
    Class representing an HDAWG instrument
    Attributes:
        channelgrouping (str) : string representation of channelgrouping (2x2 or 1x4)
        referenceclock (str) : string representation of clock used to lock HDAWG (internal or external)
        Channels (HDAWGchannel) : list of channels on the instrument
        AWGs (HDAWGawg) : list of awgs accessible in the instrument
    '''
    apilevel = 6 # Sets the 'modernity' of the api used (higher is more modern, 6 is current max version)
    defaults_path  = os.path.dirname(os.path.realpath(__file__))
    defaults_filename='defaults.xml'

    def __init__(self, device):
        '''
        Constructor for HDAWG class
        Arguments:
            device (str): device identifier (example 'dev8163')
        Returns:
            Initialized HDAWG object with 4 channels and 2 AWGs with all settings set to defaults stored in default file
        '''
        (self.daq,self.device,self.props) = ziUtils.create_api_session(device, self.apilevel) #connect to device specified by string
        ziUtils.disable_everything(self.daq,self.device) #disable all outputs of device
        self.load_default()
        self.nodepaths = self.fill_paths()
        self.Channels = []
        self.AWGs     = []
        for i in range(4):
            self.Channels.append(HDAWGchannel(self.daq,i))
        for i in range(2):
            self.AWGs.append(HDAWGawg(self.daq,i))

    def fill_paths(self):
        '''
        Function to fill in dictionary of nodes on HDAWG that we need to access to set settings
        Arguments:
            None
        Returns:
            Dictionary of nodes on HDAWG identified by their setting name
        '''
        nodes={}
        nodes['channelgrouping']='/{}/system/awg/channelgrouping'.format(self.device)
        nodes['referenceClock']='/{}/system/clocks/referenceclock/source'.format(self.device)
        return nodes
    
    ###################################
    # Methods
    ###################################

    def done(self):
        '''
        Closes connection to API
        '''
        self.daq.disconnectDevice(self.device)

    def save_settings(self, filename):
        '''
        Calls settingTools function to turn settings dictionary into file and save it
        Arguments:
            filename (str): filename of settings file
        '''
        sT.save_settings(self.getSettings(),filename)
    
    def show_keysettings(self):
        '''
        Calls settingTools function to display all key settings
        '''
        sT.print_settings(self.getSettings())

    def load_settings(self,filename):
        '''
        Load settings from file and apply to HDAWG
        Arguments:
            filename (str): full name of file to be loaded (provide full path if not in directory)
        '''
        settings = sT.load_settings(filename)
        self.setSettings(settings)

    ###################################
    # Settings Methods
    ###################################

    def getSettings(self):
        '''
        Function that returns key settings for HDAWG and subclasses (Channels, AWGs)
        Calls getSettings functions for AWGs and Channels respectively
        Arguments:
            None
        Returns:
            Nested dictionary with all settings. Top level contains main keywords: 'System','AWGs','Channels' 
            that contain all the relevant info for their respective categories. Channels will only be stored if 
            the user has configured them or they have been set from a settings file/ dictionary 
        '''
        self.daq.sync()
        settings = {}
        settings['System']={}
        for key in self.nodepaths.keys():
            settings['System'][key] = getattr(self,key)
        settings['AWGS']={}
        if self.channelgrouping == '1x4':
            settings['AWGS']['AWG0'] = self.AWGs[0].getSettings()
        else:
            for i in range(2):
                settings['AWGS']['AWG{}'.format(i)] = self.AWGs[i].getSettings()
        settings['Channels']={}
        for i in range(4):
            if self.Channels[i].configured:
                settings['Channels']['Channel{}'.format(i)] = self.Channels[i].getSettings()
        return settings

    def setSettings(self, settings):
        '''
        Function to set settings from a settings dictionary
        Arguments:
            settings (dict): dictionary structured as follows:
            'System':{All main HDAWG settings}
            'Channels':{Channel0:{Main channel settings}, Channel1:{Main channel settings},...}
            'AWGs':{AWGs0:{Main AWGs settings}, AWGs1:{Main AWGs settings}}
        Returns:
            Nothing
        '''
        for key in settings.keys():
            if key == 'AWGs':
                for i in range(2):
                    AWGname = 'AWG{}'.format(i)
                    if AWGname in settings['AWGs'].keys():
                        self.AWGs[i].setSettings(settings['AWGS'][AWGname])
            elif key == 'Channels':
                for i in range(4):
                    Channelname = 'Channel{}'.format(i)
                    if Channelname in settings['Channels'].keys():
                        self.Channels[i].setSettings(settings['Channels'][Channelname])
            else:
                setattr(self,key,settings['System'][key])
        self.daq.sync()

    
    ###################################
    # XML load and save methods
    ###################################

    def load_default(self):
        '''
        Loads default settings from xml file to device, path is hardcoded in class definition
        '''
        file = self.defaults_path + os.sep + self.defaults_filename
        ziUtils.load_settings(self.daq, self.device, file)

    def load_settings_xml(self,path,filename):
        '''
        Load settings from an xml file
        Arguments:
            path (str) : path to file
            filename (str) : name of file
        '''
        file = path + os.sep + filename
        ziUtils.load_settings(self.daq,self.device, file)

    def save_settings_xml(self, path, filename):
        '''
        Save all HDAWG settings to xml file
        Arguments:
            path (str): path to file
            filename (str): name of file
        '''
        file = path + os.sep + filename
        ziUtils.save_settings(self.daq, self.device, file)

    ##########################
    # Properties
    ##########################
    
    #Structure used to convert between int and string rep of channelgrouping property
    _groupInt=bidict({
        '2x2':0,
        '1x4':1
    })

    #Structure used to convert between int and string rep of reference clock property
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
            if val == '1x4':
                print('AWG settings and programs will be loaded from AWG0 (settings for AWG1 have not been overwriting but will not apply)')

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


#################################################################################################
# HDAWGawg subclass
#################################################################################################

class HDAWGawg:
    '''
    Helper class used to access the awg cores on HDAWG. Contains two trigger objects that can be
    used to configure the triggering of the awg core
    Attributes:
        index (int) : AWGcore identifier (0 or 1)
        channelgrouping (str) : string rep of AWGcore configuration (2x2 with 2 independent cores or 1x4)
        Triggers ([HDAWGtrigger]) : list of HDAWGtrigger objects
        samplerate (str) : string rep of samplerate used by awg core, must be 2.4GHz/2^n 
        (examples '2.4GHz', '1.2MHz', etc.)
    '''
    def __init__(self, daq, AWGcore, samplerate='2.4GHz',device='dev8163', channelgrouping='1x4'):
        '''
        Constructor for HDAWGawg class, initializes attributes and populates Triggers list
        Arguments:
            daq (daq) : handle to daq object from HDAWG, used to do interfacing with instrument
            AWGcore (int) : identified of AWGcore being programmed (0 or 1)
            Optional:
                samplerate (str) : string rep of samplerate (default '2.4GHz')
                device (str) : device serial number (default 'dev8163')
                channelgrouping (str) : string rep of channelgrouping (default '1x4')
        '''
        self.device          = device
        self.daq             = daq
        self.index           = AWGcore
        self.nodepaths       = self.fill_paths()
        self.channelgrouping = channelgrouping
        self.Triggers        = []
        for i in range (2):
            self.Triggers.append(HDAWGtrigger(self.daq,i,AWGcore=self.index, device = self.device))
        self.samplerate      = samplerate

    def fill_paths(self):
        '''
        Function to fill in dictionary of nodes on AWG that we need to access to set settings
        Arguments:
            None
        Returns:
            Dictionary of nodes on HDAWG identified by their setting name
        '''
        nodes = {}
        nodes['waves']='/{}/awgs/{}/waveform/waves/index'.format(self.device, self.index)
        nodes['single']='/{}/awgs/{}/single'.format(self.device,self.index)
        nodes['enable']='/{}/awgs/{}/enable'.format(self.device,self.index)
        nodes['samplerate']='/{}/awgs/{}/time'.format(self.device,self.index)
        return nodes
    
    ###################################
    # Methods
    ###################################

    def load_waveform(self, index, wave1, wave2=None, markers=None):
        '''
        Load arbitrary waveform to AWGcore and replace waveform{index}. 
        Requires waveform to be part of a playWave command
        Arguments:
            index (int) : index of waveform on AWG
            wave1 (array) : array of values to be played on channel 1
            wave2 (array) : array of values to be played on channel 2 (Optional)
            markers (array) : array of int values for marker configuration (0 both off, 1 marker 1 on, 2 marker 2 on, 3 both on)
        '''
        wave_awg = ziUtils.convert_awg_waveform(wave1, wave2, markers)
        node = self.nodepaths['waves'].replace('index',index)
        #path='/{:s}/awgs/0/waveform/waves/{:d}'.format(self.device, index)
        self.daq.setVector(node, wave_awg)

    def read_waveform(self,index, channels=1, markers_present=False):
        '''
        Read waveform generated by HDAWG and convert back into regular arrays
        Requires waveform to be part of playWave command 
        Arguments:
            index (int) : index of waveform on AWG
            Optional:
                channels (int) : number of channels wave is being played on (default 1)
                markers_present (bool) : boolean indicating whether markers were used (default False)
        Returns:
            [ch1,ch2,markers]: list of 3 arrays containing data from channel 1, 2 and markers.
            Arrays will be empty if nothing was played or set on them
        '''
        self.daq.sync()

        node = self.nodepaths['waves'].replace('index',index)
        #prepare data to be read out, unclear what this does exactly but in combination with the poll command this will read arbitrarily long data
        checktime = 0.1 #time to record data for (s)
        timeout   = 5 #timeout for poll command (ms)
        self.daq.getAsEvent(node) 

        #poll command returns dictionary so we need to access the actual array data
        pollReturn = self.daq.poll(checktime,timeout,flat=True)
        wave_awg   = pollReturn[node][0]['vector']

        #Use built in parse_awg_waveform to recover the channel data
        [ch1,ch2,markers] = ziUtils.parse_awg_waveform(wave_awg,channels,markers_present)

        return [ch1,ch2,markers]

    def load_program(self, awg_program):
        '''
        Load program to awg, copied from example code
        Arguments:
            awg_program (raw string): raw string containing program (use triple quotes to make sure formatting stays the same)
        '''
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

    def run(self):
        '''
        Turn on AWG and run code once
        '''
        self.daq.sync()
        node = self.nodepaths['single']
        self.daq.setInt(node,1)

    def run_loop(self):
        '''
        Run AWG continuously
        '''
        self.daq.sync()
        node = self.nodepaths['enable']
        self.daq.setInt(node,1)
    
    def stop(self):
        '''
        Stop AWG program
        '''
        node = self.nodepaths['single']
        self.daq.setInt(node,0)
        node = self.nodepaths['enable']
        self.daq.setInt(node,0)

    ###################################
    # Settings methods
    ###################################

    def getSettings(self):
        '''
        Function to build dictionary from all key settings of AWG object
        Returns:
            Nested dictionary:
                key AWG settings
                'Triggers':{Trigger0:{Key trigger0 settings}, etc.}
        '''
        self.daq.sync()
        settings = {}
        settings['samplerate'] = self.samplerate
        settings['Triggers'] = {}
        for i in range(2):
            trigger = 'Trigger{}'.format(i)
            if self.Triggers[i].configured:
                settings['Triggers'][trigger]=self.Triggers[i].getSettings()
        return settings
    
    def setSettings(self, settings):
        '''
        Function to set settings on AWG from settings dictionary
        Arguments:
            settings (dict): structured as follows: key AWG settings,'Triggers':{Trigger0:{Key trigger0 settings}, etc.}
        '''
        self.samplerate = settings['samplerate']
        if 'Triggers' in settings.keys():
            for i in range(2):
                trigger = 'Trigger{}'.format(i)
                if trigger in settings['Triggers'].keys():
                    self.Triggers[i].setSettings(settings['Triggers'][trigger])
                else:
                    continue
        self.daq.sync()
        
    ###################################
    # Properties
    ###################################

    #Structure to convert between int and string rep of sample rate
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

#################################################################################################
# HDAWGtrigger subclass
#################################################################################################
class HDAWGtrigger():
    '''
    Helper class to access trigger settings for AWG cores
    Attributes:
        ID (int) : trigger identifier (0 or 1)
        AWGcore (int) : AWGcore that trigger is applying to
        slope (str) : string rep of trigger slope ('rising',etc.)
        channel (str) : string rep for channel trigger is configured to listen to ('Trigger in 1',etc.)
    '''
    def __init__(self, daq, triggerID, device, AWGcore, slope='rising', channel='Trigger in 1'):
        '''
        Constructor for Trigger object
        Arguments:
            daq (daq) : handle to daq object used to interface with instrument
            triggerID (int) : trigger identifier (0 or 1)
            device (str) : device serial number
            AWGcore (int) : AWGcore this trigger is applying to
            Optional:
                slope (str) : trigger slope (default: 'rising')
                channel (str) : trigger source (default: 'Trigger in 1')
        '''
        self.device     = device
        self.daq        = daq
        self.ID         = triggerID
        self.AWGcore    = AWGcore
        self.nodepaths  = self.fill_paths(device, AWGcore, triggerID)
        self.slope      = slope
        self.channel    = channel
        self.configured = False
    
    def fill_paths(self, device, AWGcore, triggerID):
        nodes = {}
        nodes['slope']='/{}/awgs/{}/auxtriggers/{}/slope'.format(device, AWGcore, triggerID)
        nodes['channel']='/{}/awgs/{}/auxtriggers/{}/channel'.format(device, AWGcore, triggerID)
        return nodes

    ###################################
    # Methods
    ###################################

    def configureTrigger(self,slope, channel):
        '''
        Configure trigger to provided slope and channel
        Arguments:
            slope (str) : string rep of trigger slope
            channel (str) : string rep of trigger source
        '''
        self.slope      = slope
        self.channel    = channel
        self.configured = True
        self.daq.sync()

    ###################################
    # Settings Methods
    ###################################

    def getSettings(self):
        '''
        Create dictionary of key trigger settings
        Returns:
            settings (dict) : dictionary of trigger settinsg
        '''
        self.daq.sync()
        settings={}
        for key in self.nodepaths.keys():
            settings[key]=getattr(self,key)
        return settings

    def setSettings(self,settings):
        '''
        Set key trigger attributes from dictionary
        Arguments:
            settings (str) : dictionary of trigger settings 
        '''
        for key in settings.keys():
            setattr(self,key,settings[key])
        self.configured = True
        self.daq.sync()

    ###################################
    # Properties
    ###################################

    #Structure to convert between int and string rep of slope property
    _slopeInt = bidict({
        'level':0,
        'rising':1,
        'falling':2,
        'both':3
    })

    #Structure to convert between int and string rep of channel property
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
            self.configured = True

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
            self.configured = True

#################################################################################################
# HDAWGchannel subclass
#################################################################################################
class HDAWGchannel():
    '''
    Helper class used to access all the output settings of the 4 HDAWG channels
    Attributes:
        ID (int) : channel number
        AWGcore (int) : AWGcore controlling output
        AWGout (int) : channel number representation for AWGcore (i.e. channel 2 is on AWG1/channel0 etc.)
        status (str) : string rep of channel state ('On'/'Off')
        fullscale (double) : output voltage range on channel (V)
        AWGamp (double) : scaling factor relative to fullscale for output voltage 
        amp (double) : output voltage amplitude of channel (fullscale*AWGamp). If set, will automatically set fullscale and AWGamp to optimal settings
        offset (double) : output voltage offset (V)
        delay (double) : delay applied to channel (s)
        markers (str) : determines whether this channel has markers associated with it
    '''

    def __init__(self, daq, channelID, device='dev8163', amp = 1.0, fullscale = 1.0, AWGamp = 1.0, offset = 0.0, delay = 0.0, marker_out = 'Trigger', hold = 'False'):
        '''
        Constructor for channel class.
        Arguments:
            daq (daq) : handle to daq used to interface with hardware
            channelID : channel number
            Optional:
                device: device serial number (default: 'dev8163')
                amp: output amplitude in V (default: 1)
                fullscale: channel fullscale setting in V (default: 1)
                AWGamp: awg scaling factor applied on fullscale (default: 1)
                delay: delay applied on channel to align waves (default: 0)
                offset: offset applied to output of channel (default: 0)
                markers: determines whether channel has markers associated with it
        '''
        self.device     = device
        self.daq        = daq
        self.ID         = channelID
        self.AWGcore    = int(np.floor(channelID/2))
        self.AWGout     = channelID%2
        self.nodepaths  = self.fill_paths()
        self.status     = 'Off'
        self._notinit    = False
        self.fullscale  = fullscale
        self.AWGamp     = AWGamp
        self.amp        = amp
        self.offset     = offset
        self.delay      = delay
        self.marker     = marker_out
        self.hold       = hold
        self._notinit   = True
        self.configured = False

    def fill_paths(self):
        '''
        Helper function to fill paths to nodes in hardware containing settings of interest
        '''
        nodes = {}
        nodes['fullscale'] = '/{}/sigouts/{}/range'.format(self.device,self.ID)
        nodes['offset'] = '/{}/sigouts/{}/offset'.format(self.device,self.ID)
        nodes['delay'] = '/{}/sigouts/{}/delay'.format(self.device,self.ID)
        nodes['AWGamp'] = '/{}/awgs/{}/outputs/{}/amplitude'.format(self.device, self.AWGcore, self.AWGout)
        nodes['marker'] = '/{}/triggers/out/{}/source'.format(self.device,self.ID)
        nodes['status'] = '/{}/sigouts/{}/on'.format(self.device,self.ID)
        nodes['hold'] = '/{}/awgs/{}/outputs/{}/hold'.format(self.device,self.AWGcore,self.AWGout)
        return nodes

    ###################################
    # Methods
    ###################################

    def configureChannel(self, amp=1.0, fullscale=1.0, AWGamp=1.0, delay=0.0, offset=0.0, marker_out = 'Trigger', hold = 'False'):
        '''
        Configure channel output settings. If amp is provided, ignores fullscale and AWGamp parameters
        Arguments:
            Optional:
                device: device serial number (default: 'dev8163')
                amp: output amplitude in V (default: 1)
                fullscale: channel fullscale setting in V (default: 1)
                AWGamp: awg scaling factor applied on fullscale (default: 1)
                delay: delay applied on channel to align waves (default: 0)
                offset: offset applied to output of channel (default: 0)
                markers: determines whether channel has markers associated with it
        '''
            
        if amp < 1.0 or amp > 1.0:
            self.amp = amp
        else:
            self.fullscale = fullscale
            self.AWGamp    = AWGamp
            print('Manual range configuration')

        self.delay      = delay
        self.offset     = offset
        self.marker     = marker_out
        self.hold       = hold
        self.configured = True
        self.status     = 'On'
        self.daq.sync()

    ###################################
    # Settings Methods
    ###################################

    def getSettings(self):
        '''
        Build dictionary from all key channel settings
        Returns:
            settings (dict) : dictionary where keys are attributes of channel object
        '''
        self.daq.sync()
        settings={}
        settings['amp'] = getattr(self,'amp')
        for key in self.nodepaths.keys():
            settings[key] = getattr(self,key)
        return settings

    def setSettings(self, settings):
        '''
        From settings dictionary, configure channel
        Arguments:
            settings (dict) : dictionary containing key channel settings
            If settings are missing, Channel will keep default values
        '''
        for key in settings:
            setattr(self,key,settings[key])
        self.configured = True
        self.status     = 'On'
        self.daq.sync()

    ###################################
    # Properties
    ###################################

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
        if self._notinit:
            print('Automatically adjusting fullscale and AWGamp values')
        ranges = [0.2,0.4,0.6,0.8,1.0,2.0,3.0,4.0,5.0]
        for x in ranges:
            if val > x:
                continue
            else:
                self.fullscale = x
                self.AWGamp    = val/x
                break
        self.configured = True

    _markerInt = bidict({
        'Trigger':0,
        'Marker':4
    })

    @property
    def marker(self):
        node   = self.nodepaths['marker']
        status = self.daq.getInt(node)
        val    = status-self.ID
        output = self._markerInt.inverse[val]   
        print('Using {} function on marker channel'.format(output))
        return output
    @marker.setter
    def marker(self,val):
        node      = self.nodepaths['marker']
        if val not in self._markerInt.keys():
            print('Error, acceptable inputs are {}'.format(list(self._markerInt)))
        else:
            markindex = self.ID+self._markerInt[val]
            print('Setting marker output to {}'.format(val))
            self.daq.setInt(node,markindex)
            self.configured = True

    @property
    def hold(self):
        node = self.nodepaths['hold']
        val  = self.daq.getDouble(node)
        if val:
            return 'True'
        else:
            return 'False'
    @hold.setter
    def hold(self,val):
        node = self.nodepaths['hold']
        if val not in ['True','False']:
            print('Hold setting is either "True" or "False", invalid input')
        else:
            if val:
                self.daq.setDouble(node,1)
            else:
                self.daq.setDouble(node,0)
    
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
        self.configured = True
    
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
        self.configured = True

    @property
    def offset(self):
        node = self.nodepaths['offset']
        val  = self.daq.getDouble(node)
        return val
    @offset.setter
    def offset(self,val):
        node = self.nodepaths['offset']
        self.daq.setDouble(node,val)
        self.configured = True

    @property
    def delay(self):
        node = self.nodepaths['delay']
        val  = self.daq.getDouble(node)
        return val
    @delay.setter
    def delay(self,val):
        node = self.nodepaths['delay']
        self.daq.setDouble(node,val)
        self.configured = True