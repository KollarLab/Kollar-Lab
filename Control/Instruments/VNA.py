'''
Author: Martin Ritter
Date: 10/19/2020
VNA class
'''
from time import sleep, time

import numpy
import pyvisa

class VNA():
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
    def __init__(self, address):
        '''
        Initialize connection with instrument, reset it and clear all errors
        '''
        rm = pyvisa.ResourceManager()
        self.inst = rm.open_resource(address)
        self.inst.write('*RST; *CLS')
        self.inst.write('OUTP OFF')

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

    def get_errors(self):
        '''Print out all errors and clear queue'''
        code, *string = self.error
        while int(code) != 0:
            print('{},{}'.format(code, string))
            code, *string = self.error
        print('{},{}'.format(code, string))
        self.inst.write('*CLS')

    def spec_default_settings(self):
        '''Return default settings for spec measurement'''
        settings = {}

        settings['channel'] = 1
        settings['averages'] = 100
        settings['measurement'] = 'S23'
        settings['meas_form'] = 'MLOG'
        settings['start'] = '1 MHz'
        settings['stop'] = '40 MHz'
        settings['sweep_points'] = 501
        settings['RFpower'] = -10
        settings['RFport'] = 1
        settings['Mport'] = 2
        settings['CAVport'] = 3
        settings['CAVpower'] = -5
        settings['CAVfreq'] = '60 MHz'

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
        averages     = settings['averages']
        measurement  = settings['measurement']
        meas_format  = settings['meas_form']
        start        = settings['start']
        stop         = settings['stop']
        sweep_points = settings['sweep_points']
        rf_power     = settings['RFpower']
        rf_port      = settings['RFport']
        m_port       = settings['Mport']
        cav_port     = settings['CAVport']
        cav_power    = settings['CAVpower']
        cav_freq     = settings['CAVfreq']

        #Turn off output and switch to single sweep mode instead of continuous sweeps
        self.output = 'OFF'
        self.inst.write('INIT:CONT OFF')

        #Configure averaging
        self.configure_averages(channel, averages)

        #Clear old traces on channels and define a new trace to measure 'S23'
        self.configure_trace(channel, 'spec', measurement, meas_format)

        #Configure frequency sweep and RF power
        self.configure_frequency(channel, start, stop, sweep_points)
        self.power = rf_power

        #Configure ports
        #RF
        #configures rf_port to always be on during measurements
        self.inst.write('SOUR{}:POW{}:PERM ON'.format(channel, rf_port))
        #Measure
        #configure to only receive
        self.inst.write('SOUR{}:POW{}:STATE OFF'.format(channel, m_port))
        #only look at what is happening at cav_freq
        self.inst.write('SOUR{}:FREQ{}:CONV:ARB:IFR 1, 1, {}, FIX'.format(channel, m_port, cav_freq))
        #LO
        #always keep cav_port on during measurements
        self.inst.write('SOUR{}:POW{}:PERM ON'.format(channel, cav_port))
        #only output at 60MHz
        self.inst.write('SOUR{}:FREQ{}:CONV:ARB:IFR 1, 1, {}, FIX'.format(channel, cav_port, cav_freq))
        #fixed power of cav_power (ignores the power level of Ch1)
        self.inst.write('SOUR{}:POW{}:OFFS {}, ONLY'.format(channel, cav_port, cav_power))
        self.inst.query('*OPC?')

        self.get_errors()

        self.wait_complete(channel, averages)

        self.autoscale(window=1, ref_trace='spec')

    def wait_complete(self, channel, averages):
        '''
        Block execution until instrument has completed measurement
        Will estimate the sweep time and update the user periodically on progress
        Arguments:
            channel (int): channel number
            averages (int): number of traces the instrument is set to average
        '''
        sweep_time = float(self.inst.query('SENS{}:SWE:TIME?'.format(channel)))
        total_time = 2.5*averages*sweep_time

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
            sleep(min(max(sweep_time*0.2*averages, 5), total_time))
            if auto_scale:
                self.autoscale()
                auto_scale = False

    def trans_default_settings(self):
        '''Return default settings for trans measurement'''
        settings = {}

        settings['channel'] = 1
        settings['averages'] = 100
        settings['measurement'] = 'S21'
        settings['start'] = '1 MHz'
        settings['stop'] = '40 MHz'
        settings['sweep_points'] = 501
        settings['RFpower'] = -10

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
        averages     = settings['averages']
        measurement  = settings['measurement']
        start        = settings['start']
        stop         = settings['stop']
        sweep_points = settings['sweep_points']
        rf_power      = settings['RFpower']

        #Turn off output and switch to single sweep mode instead of continuous sweeps
        self.output = 'OFF'
        self.inst.write('INIT:CONT OFF')

        #Configure averaging
        self.configure_averages(channel, averages)

        #Clear old traces on channels and define a new trace to measure 'S21'
        self.configure_trace(channel, 'mag', measurement, 'MLOG')
        self.configure_trace(channel, 'phase', measurement, 'UPH')
        self.display_trace(trace='mag', tracenum=1, window=1)
        self.display_trace(trace='phase', tracenum=2, window=1)

        #Configure frequency sweep and RF power
        self.configure_frequency(channel, start, stop, sweep_points)
        self.power = rf_power

        self.get_errors()

        #Query instrument until operation is complete
        self.wait_complete(channel, averages)

        self.autoscale(window=1, ref_trace='mag')
        self.autoscale(window=1, ref_trace='phase')

        mag   = self.get_trace(channel, 'mag')
        phase = self.get_trace(channel, 'phase')
        freqs = self.get_channel_axis(channel)

        return mag, phase, freqs

    def autoscale(self, window=1, ref_trace=None):
        '''Autoscale window to match ref_trace range'''
        tracestr = ''
        if ref_trace is not None:
            tracestr = ",'{}'".format(ref_trace)
        self.inst.write("DISP:WIND{}:TRAC:Y:AUTO ONCE{}".format(window, tracestr))

    def configure_averages(self, channel, averages):
        '''Set up channel to measure averages traces'''
        self.inst.write('SENS{}:AVER ON'.format(channel))
        self.inst.write('SENS{}:AVER:COUN {}'.format(channel, averages))
        self.inst.write('SENS{}:AVER:CLE'.format(channel))
        self.inst.write('SENS{}:SWE:COUN {}'.format(channel, averages))
        self.inst.query('*OPC?')

    def configure_frequency(self, channel, start, stop, sweep_points):
        '''
        Set up frequency axis
        Arguments:
            channel(int): channel number, channel must have been created beforehand
            start(float/str): frequency of sweep start ('10 MHz' or 10e6)
            stop(float/str): frequency of sweep end ('100 MHz' or 100e6)
            sweep_points(int): number of points in sweep
        '''
        self.inst.write('SENS{}:FREQ:STAR {}; STOP {}'.format(channel, start, stop))
        self.inst.write('SENS{}:SWE:POIN {}'.format(channel, sweep_points))
        self.inst.query('*OPC?')

    def clear_all_traces(self):
        '''Clear all traces defined on instrument'''
        self.inst.write('CALC:PAR:DEL:ALL')
        self.inst.query('*OPC?')

    def clear_channel_traces(self, channel):
        '''Clear all traces defined in channel'''
        self.inst.write('CALC{}:PAR:DEL:CALL'.format(channel))
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

    def get_channel_traces(self, channel):
        '''Return list of traces defined in channel'''
        instresp  = self.inst.query('CALC{}:PAR:CAT?'.format(channel))
        tracelist = instresp.strip("\n'").split(',')
        return tracelist

    def get_channel_axis(self, channel):
        '''Return numpy array of channel x-axis'''
        xaxis = self.inst.query_ascii_values("CALC{}:DATA:STIM?".format(channel))
        return numpy.asarray(xaxis)

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

    def close(self):
        '''Close VISA connection'''
        self.inst.close()


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
#    Run single sweep:
#        inst.write('INIT2:CONT OFF')
#        inst.write('INIT2:IMM')
#
##Averaging stuff:
#
#    Enable/disable:
#        inst.write('SENS1:AVER ON')
#    Clear average:
#        inst.write('SENS1:AVER:CLE')
#    Number of averages:
#        inst.write('SENS1:AVER:COUN 15')
#    Average mode:
#        inst.write('SENS1:AVER:MODE AUTO' )
#
##Frequency span control:
#        inst.write('FREQ:STAR 1 MHz')
#        inst.write('FREQ:STOP 40 MHz')
#        inst.write('SWE:POIN 501')
#        inst.query('SWE:TIME?')

## Physical port control:
#        inst.write('SOUR<Channel>:FREQ<Port>:CONV:ARB:IFR num,denom, offset, FIX')
#        inst.write('SOUR<Channel>:POW<Port>:STATE ON') turn RF on/off
#        inst.write('SOUR<Channel>:POW<Port>:PERM ON') turn 'gen' on/off
#        inst.write('SOUR<Channel>:POW<Port>:OFFS offset, ONLY')

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