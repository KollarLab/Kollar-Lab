'''
Author: Martin Ritter
Date: 10/19/2020
VNA class
'''
from time import sleep, time

import numpy
import pyvisa

import matplotlib.pyplot as plt

from Instruments.SCPIinst import SCPIinst

class VNA(SCPIinst):
    '''
    Class representing RS VNA instrument
    Attributes:
        output (str or int): output state of instrument (ON/1 or OFF/0)
        power (str or float): output power of instrument ('10 dBm' or 10 work)
        error: read only value, returns error code and string
    Methods:
        __init__(address): initialize connection and resets the instrument
        get_errors(): return all errors form instrument and clears them
        spec_default_settings(): return default settings for a spec measurement
        meas_spec(settings): perform spec measurement using settings dict
        trans_default_settings(): return default settings for transmission meas
        meas_trans(settings): perform transmission measurement using settings
        wait_complete(channel, averages): block execution until acquisition is done
        autoscale(window, ref_trace): autoscale window using ref_trace as reference
        configure_averages(channel, averages): setup VNA to meas and average x times
        configure_frequency(channel, start, stop, sweep_points): configure freq axis
        configure_trace(channel, name, meastype, measformat): set up trace on a
            channel to measure 'meastype' (e.g.'S21') in 'measformat' (e.g.'MLOG')
        get_channel_traces(channel): get a list of traces on the channel
        get_channel_axis(channel): return numpy array of the x-axis of channel
        get_trace(channel, name): return numpy array of trace data on channel
        clear_all_traces(): remove all traces defined on all channels
        clear_channel_traces(channel): remove traces defined in channel
        display_trace(trace, tracenum, window): display trace in window with number
            tracenum
        close(): close VISA connection with instrument
    '''
    errcmds           = {}
    errcmds['error']  = 'SYST:ERR?'
    errcmds['serror'] = 'SYST:SERR?'
    
    commandlist = {}
    commandlist['core']   = {}
    commandlist['Ref']    = {}
    
    core = {}
    core['Output'] = 'OUTPut:STATe'
    core['Power']  = 'SOURce:POWer'
    core['SweepType'] = 'SENS:SWE:TYPE'
    core['Freq']   = 'SOURce:FREQuency:CW'
    core['ifBW']   = 'SENS:BAND'
    Ref = {}
    Ref['Source']    = 'ROSCillator'
    Ref['Frequency'] = 'ROSCillator:EXTernal:FREQuency'
    
    commandlist['core']   = core
    commandlist['Ref']    = Ref


    def __init__(self, address, reset = True):
        '''
        Initialize connection with instrument, reset it and clear all errors
        '''
        super().__init__(address, self.commandlist, self.errcmds, reset) 
        #rm = pyvisa.ResourceManager()
        #self.inst = rm.open_resource(address)
        self.ext_ref = 'EXT'
        self.Output = 'Off'
        #if reset:
        #    self.reset()

    @property
    def error(self):
        '''Return error code and string'''
        return (self.inst.query('SYST:ERR?')).strip().split(',')

    @property
    def output(self):
        '''Return output state of instrument ('ON'/'OFF')'''
        state = int(self.inst.query('OUTP?'))
        if state:
            return 'ON'
        return 'OFF'
    @output.setter
    def output(self, state):
        '''
        Set output state of instrument
        Arguments:
            state(str/int): desired state, can be a string ('ON') or an int (1)
        '''
        self.inst.write('OUTP {}'.format(state))
    
    @property
    def freq(self):
        return float(self.inst.query('FREQ:CW?'))
    @freq.setter
    def freq(self, freq):
        self.inst.write('FREQ:CW {}'.format(freq))
        
    @property
    def power(self):
        '''Return output RF power'''
        return float(self.inst.query('SOUR:POW?'))
    @power.setter
    def power(self, power):
        '''
        Set output RF power
        Arguments:
            power(str/float): desired output power, can be string ('-10 dBm') or float (-11.2')
        '''
        self.inst.write('SOUR:POW {}'.format(power))

    @property
    def ifBW(self):
        '''Return IF bandwidth'''
        return float(self.inst.query('SENS:BAND?'))
    @ifBW.setter
    def ifBW(self, ifBW):
        '''
        Set IF bandwidth
        Argument:
            ifBW(str/float): IF BW to use, can be a string ('1 kHz') of a float (1e3)
            the VNA will auto round to 1,1.5,2,3,5,7 in the nearest decade
        '''
        self.inst.write('SENS:BAND {}'.format(ifBW))
        
    @property
    def ext_ref(self):
        '''Return reference source'''
        return self.inst.query('ROSC?').strip()
    @ext_ref.setter
    def ext_ref(self, source):
        '''
        Set reference source
        Argument:
            source(str): 'EXT' or 'INT' for external or internal freq reference
        '''
        self.inst.write('ROSC {}'.format(source))
        
    def get_errors(self, verbose = True):
        '''Print out all errors and clear queue'''
        code, *string = self.error
        while int(code) != 0:
            print('{},{}'.format(code, string))
            code, *string = self.error
        if verbose:
            print('{},{}'.format(code, string))
        self.inst.write('*CLS')

    def reset(self):
        self.inst.write('*RST; *CLS')
        self.inst.write('OUTP OFF')
        self.ext_ref = 'EXT'
    
    def close(self):
        '''Close VISA connection'''
        self.inst.close()
    
    ### User functions (spec, trans, powersweep etc.)     
    def spec_default_settings(self):
        '''Return default settings for spec measurement'''
        settings = {}

        settings['channel']      = 1
        settings['avg_time']     = 10
        settings['measurement']  = 'S21'
        settings['start_freq']   = '3 GHz'
        settings['stop_freq']    = '6 GHz'
        settings['freq_points'] = 501
        settings['RFpower']      = -10
        settings['RFport']       = 3
        settings['Mport']        = 2
        settings['CAVport']      = 1
        settings['CAVpower']     = -5
        settings['CAVfreq']      = '7 GHz'
        settings['ifBW']         = 1e3
        
        return settings

    def spec_meas(self, settings):
        '''
        Perform spec measurement using settings dict.
        Measurement consists of:
            -setting up an cavity port (fixed frequency/power), CAVport
            -setting up a probe port (swept freq), RFport
            -configuring the measurement port (Mport) to measure
                S(Mport-CAVport) as function of RFfreq
        Currently just plots data on VNA
        '''
        channel      = settings['channel']
        time         = settings['avg_time']
        measurement  = settings['measurement']
        start        = settings['start_freq']
        stop         = settings['stop_freq']
        sweep_points = settings['freq_points']
        rf_power     = settings['RFpower']
        rf_port      = settings['RFport']
        m_port       = settings['Mport']
        cav_port     = settings['CAVport']
        cav_power    = settings['CAVpower']
        cav_freq     = settings['CAVfreq']
        ifBW         = settings['ifBW']

        #Turn off output and switch to single sweep mode instead of continuous sweeps
        self.output = 'OFF'
        self.inst.write('INIT:CONT OFF')

        #Configure averaging
        self.configure_averages(channel, 1e4)

        #Clear old traces on channels and define a new trace to measure 'S21'
        self.configure_measurement(channel, measurement)
  
        #Configure frequency sweep and RF power
        self.configure_frequency(channel, start, stop, sweep_points)
        self.power = rf_power
        self.ifBW = ifBW

        #Configure ports
        #RF
        #configures rf_port to always be on during measurements
        self.inst.write('SOUR{}:POW{}:PERM ON'.format(channel, rf_port))
        #Measure
        #configure to only receive
#        self.inst.write('SOUR{}:POW{}:STATE OFF'.format(channel, m_port))
        #only look at what is happening at cav_freq
        self.inst.write('SOUR{}:FREQ{}:CONV:ARB:IFR 1, 1, {}, FIX'.format(channel, m_port, cav_freq))
        #LO
        #always keep cav_port on during measurements
#        self.inst.write('SOUR{}:POW{}:PERM ON'.format(channel, cav_port))
        #only output at 60MHz
        self.inst.write('SOUR{}:FREQ{}:CONV:ARB:IFR 1, 1, {}, FIX'.format(channel, cav_port, cav_freq))
        #fixed power of cav_power (ignores the power level of Ch1)
        self.inst.write('SOUR{}:POW{}:OFFS {}, ONLY'.format(channel, cav_port, cav_power))
        self.inst.query('*OPC?')

        self.get_errors(verbose = False)
        self.avg_time(channel, time)
        self.autoscale(window=1)

        data = self.get_channel_data(channel)
#        self.output = 'Off'
        
        return data

    def trans_default_settings(self):
        '''Return default settings for trans measurement'''
        settings = {}

        settings['channel']      = 1
        settings['avg_time']     = 10
        settings['measurement']  = 'S21'
        settings['start_freq']        = '3 GHz'
        settings['stop_freq']         = '9 GHz'
        settings['freq_points'] = 501
        settings['RFpower']      = -20
        settings['ifBW']         = 1e3

        return settings

    def trans_meas(self, settings):
        '''
        Perform transmission measurement between user specified ports
        Uses settings dictionary to configure instrument
        Returns:
            mag: array containing magnitude info for measurement
            phase: array containing phase info
            freqs: x axis for channel (frequency here)
        '''
        channel      = settings['channel']
        time         = settings['avg_time']
        measurement  = settings['measurement']
        start        = settings['start_freq']
        stop         = settings['stop_freq']
        sweep_points = settings['freq_points']
        rf_power     = settings['RFpower']
        ifBW         = settings['ifBW']

        #Turn off output and switch to single sweep mode instead of continuous sweeps
        self.output = 'OFF'
        self.inst.write('INIT:CONT OFF')
        
#        self.inst.write('SOUR{}:POW{}:PERM ON'.format(channel, 1))
#        self.inst.write('SOUR{}:POW{}:STATE OFF'.format(channel, 2))
        #Configure averaging
        self.configure_averages(channel, 1e4) #high number so that VNA keeps averaging regardless of user averaging time
        self.configure_measurement(channel, measurement)

        #Configure frequency sweep and RF power
        self.configure_frequency(channel, start, stop, sweep_points)
        self.power = rf_power
        self.ifBW = ifBW

        self.get_errors(verbose = False)

        self.avg_time(channel, time)

        data = self.get_channel_data(channel)
#        self.output = 'Off'
        
        return data

    def power_sweep(self, settings, powers):

        channel = settings['channel']
        ifBW = settings['ifBW']
        sweep_points = settings['sweep_points']
        averages = settings['averages']
        freq = settings['frequency']
        span = settings['span']
        
        self.ifBW = ifBW
        
        self.configure_trace(channel, 'Trc1', 'S21', 'MLOG')
        self.configure_trace(channel, 'Trc2', 'S21', 'UPH')
        self.display_trace('Trc1',1,1)
        self.display_trace('Trc2',2,1)
        
        self.configure_frequency(channel, center=freq, span=span, sweep_points=sweep_points)
        
        mags = []
        phases = []
        axes = []
        
        for power in powers:
            self.output = 'OFF'
            self.inst.write('INIT:CONT OFF')
            self.power = power
            self.configure_averages(channel, averages)
            self.wait_complete(channel, averages)
            mag = self.get_trace(channel, 'Trc1')
            phase = self.get_trace(channel, 'Trc2')
            freq = self.get_channel_axis(channel)
            mags.append(mag)
            phases.append(phase)
            axes.append(freq)
        
        return mags, phases, axes
    
    def meas_lines(self, freqs, spans, settings):
        
        power = settings['RFpower']
        channel = settings['channel']
        ifBW = settings['ifBW']
        sweep_points = settings['sweep_points']
        time = settings['avg_time']
        
        self.power = power
        self.ifBW = ifBW
        
#        self.configure_trace(channel, 'Trc1', 'S21', 'MLOG')
#        self.configure_trace(channel, 'Trc2', 'S21', 'UPH')
#        self.display_trace('Trc2',2,1)
#        self.display_trace('Trc1',1,1)
        self.configure_measurement(channel, 'S21')
        
        
        mags = []
        phases = []
        axes = []

        for freq, span in list(zip(freqs, spans)):
            print('{},{}'.format(freq, span))
            self.output = 'OFF'
            self.inst.write('INIT:CONT OFF')
            self.configure_frequency(channel, center=freq, span=span, sweep_points=sweep_points)
            self.configure_averages(channel, 10e3)
#            self.wait_complete(channel, averages)
            self.avg_time(1, time)
            data = self.get_channel_data(channel)
            mag = data['mag']
            phase = data['phase']
            freq = data['xaxis']
            mags.append(mag)
            phases.append(phase)
            axes.append(freq)
        
        return mags, phases, axes
    
    ### Helpers and utility functions
    def autoscale(self, window=1, ref_trace=None):
        '''Autoscale window to match ref_trace range'''
        tracestr = ''
        if ref_trace is not None:
            tracestr = ",'{}'".format(ref_trace)
            self.inst.write("DISP:WIND{}:TRAC:Y:AUTO ONCE{}".format(window, tracestr))
        if ref_trace is None:
            traces = self.get_channel_traces(1)[::2]
            for trace in traces:
                tracestr = ",'{}'".format(trace)
                self.inst.write("DISP:WIND{}:TRAC:Y:AUTO ONCE{}".format(window, tracestr))
                
    def avg_time(self, channel, time):
        '''
        Set the instrument to measure for 'time' seconds. The instrument will
        autoscale after 20% of the time has elapsed so that the traces can be
        examined on the VNA
        Arguments:
            channel (int): channel being used for the measurement
            time (int): time in seconds that the instrument will measure
        '''
        self.clear_averages(channel)
        self.output = 'On'
        self.inst.write('INIT{}:IMM'.format(channel))
        sleep(time/5)
        self.autoscale()
        sleep(4*time/5)
#        self.inst.write('SENS:CORR:EDEL:AUTO ONCE')
        
    def clear_all_traces(self):
        '''Clear all traces defined on instrument'''
        self.inst.write('CALC:PAR:DEL:ALL')
        self.inst.query('*OPC?')
        
    def clear_averages(self, channel):
        '''Clear averages on channel'''
        self.inst.write('SENS{}:AVER:CLE'.format(channel))

    def clear_channel_traces(self, channel):
        '''Clear all traces defined in channel'''
        self.inst.write('CALC{}:PAR:DEL:CALL'.format(channel))
        self.inst.query('*OPC?')
        
    def configure_averages(self, channel, averages):
        '''Set up channel to measure averages traces'''
        self.inst.write('SENS{}:AVER ON'.format(channel))
        self.inst.write('SENS{}:AVER:COUN {}'.format(channel, averages))
        self.inst.write('SENS{}:SWE:COUN {}'.format(channel, averages))
        self.clear_averages(channel)
        self.inst.query('*OPC?')
        
    def configure_frequency(self, channel, start=None, stop=None, sweep_points=None, center=None, span=None):
        '''
        Set up frequency axis
        Arguments:
            channel(int): channel number, channel must have been created beforehand
            if specified:
                start(float/str): frequency of sweep start ('10 MHz' or 10e6)
                stop(float/str): frequency of sweep end ('100 MHz' or 100e6)
            if specified:
                center(float/str): center of frequency sweep
                span(float/str): span of frequency sweep
            sweep_points(int): number of points in sweep
        '''
        if sweep_points is None:
            sweep_points = 501
        if start is not None and stop is not None:
            self.inst.write('SENS{}:FREQ:STAR {}; STOP {}'.format(channel, start, stop))
            self.inst.write('SENS{}:SWE:POIN {}'.format(channel, sweep_points))
        if center is not None and span is not None:
            self.inst.write('SENS{}:FREQ:CENT {}; SPAN {}'.format(channel, center, span))
            self.inst.write('SENS{}:SWE:POIN {}'.format(channel, sweep_points))
        self.inst.query('*OPC?')
    
    def configure_measurement(self, channel, measurement, window = 1):
        '''
        Set up a basic 'Sij' measurement configuring traces for magnitude and 
        phase. If the traces already exist, the code will ignore the clearing
        commands and return without doing anything
        Arguments:
            channel (int): channel number for measurement
            measurement (string): measurement to be performed (e.g 'S21')
            window (int): window to display the traces in (default 1)
        '''
        traces = self.get_channel_traces(channel)
        reinit = True
        if 'mag' in traces and 'phase' in traces:
            reinit = False
        if reinit:
            self.clear_channel_traces(channel)
            self.configure_trace(channel, 'mag', measurement, 'MLOG')
            self.configure_trace(channel, 'phase', measurement, 'PHAS')
            self.display_trace('mag',tracenum=1, window=window)
            self.display_trace('phase',tracenum=2, window=window)

    def configure_trace(self, channel, name, meastype, measformat):
        '''
        Configure trace to measure a specific parameter in a given format
        Argument:
            channel(int): channel number
            name(str): name to assign to trace
            meastype(str): parameter being measured (e.g. 'S21')
            measformat(str): format to measure parameter (e.g. 'MLOG' for dB scale)
        '''
        self.inst.write("CALC{}:PAR:SDEF '{}', '{}'".format(channel, name, meastype))
        self.inst.write("CALC{}:FORM {}".format(channel, measformat))
        self.inst.query('*OPC?')
    
    def display_trace(self, trace, tracenum, window):
        '''
        Display trace in window assigning it the tracenum number
        Arguments:
            trace (str): name of trace to be displayed
            tracenum (int): number to assign to trace in window
            window (int): window to display trace in (will create a new one if needed)
        '''
        self.inst.write('DISP:WIND{}:STAT ON'.format(window))
        self.inst.query('*OPC?')
        self.inst.write('DISP:WIND{}:TRAC{}:FEED "{}"'.format(window, tracenum, trace))
        self.inst.query('*OPC?')

    def get_channel_axis(self, channel):
        '''Return numpy array of channel x-axis'''
        xaxis = self.inst.query_ascii_values("CALC{}:DATA:STIM?".format(channel))
        return numpy.asarray(xaxis)
    
    def get_channel_data(self, channel):
        '''Return dictionary with all the data from a channel'''
        traces = self.get_channel_traces(channel)[::2]
        data = {}
        for trace in traces:
            data[trace] = self.get_trace(channel, trace)
        data['xaxis'] = self.get_channel_axis(channel)
        return data
    
    def get_channel_traces(self, channel):
        '''Return list of traces defined in channel'''
        instresp  = self.inst.query('CALC{}:PAR:CAT?'.format(channel))
        tracelist = instresp.strip("\n'").split(',')
        return tracelist

    def get_trace(self, channel, name):
        '''
        Return measurement data for trace with name 'name'
        Arguments:
            channel(int): channel number
            name(str): trace name, same as name given at trace creation
        '''
        tracelist = self.get_channel_traces(channel)[::2]
        if name not in tracelist:
            print('Trace not found, possible options are: {}'.format(tracelist))
            return numpy.zeros(1000)
        data = self.inst.query_ascii_values("CALC{}:DATA:TRAC? '{}',FDAT".format(channel, name))
        return numpy.asarray(data)
    
    ####Legacy functions/ unused functions
    def wait_complete(self, channel, averages):
        '''
        Block execution until instrument has completed measurement
        Will estimate the sweep time and update the user periodically on progress
        Arguments:
            channel (int): channel number
            averages (int): number of traces the instrument is set to average
        '''
        sweep_time = float(self.inst.query('SENS{}:SWE:TIME?'.format(channel)))
        total_time = averages*sweep_time

        #Turn on output and wait until sweep is complete to autoscale the trace
        self.output = 'ON'
        self.inst.write('INIT{}:IMM; *OPC'.format(channel))

        opc = 0
        auto_scale = True
        t1 = time()
        while not opc&1:
            t2 = time()
            print('waiting for command to complete, elapsed time:{}, total time est. {}'.format(t2-t1, total_time))
            opc = int(self.inst.query('*ESR?').strip("\n'"))
            if opc&1:
                break
            sleep(min(max(sweep_time*0.2*averages, 5), total_time))
            if auto_scale:
                self.autoscale()
                auto_scale = False


##Basic commands:
#    Query error:
#        inst.query('SYST:ERR?')
#    inst.write('*RST')
##Doing S measurements:
#
#    Clear all traces:
#        inst.write('CALC:PAR:DEL:ALL')
#    List all traces in channel:
#        inst.query('CALC2:PAR:CAT?')
#    Make s param group:
#        inst.write("CALC2:PAR:DEF:SGR 1,2")
#        Display results on screen:
#            inst.write("DISP:WIND:TRAC2:FEED 'Ch2_SG_S11'")
#    Read data from s param group:
#        sdat = inst.query("CALC2:DATA:SGR? FDAT")

##Other useful stuff:
#
#    Bandwidth:
#        inst.write('SENS1:BAND:RES 1')
#        inst.write('SENS1:BAND:RES:SEL NORM')
#
#    Electrical delay:
#        inst.write('SENS2:CORR:EDEL:AUTO ONCE')
#        inst.write('SENS2:CORR:LOSS:AUTO ONCE') (loss and edelay)
#        inst.write('SENS2:CORR:OFFS:COMP ON') (toggle compensation)