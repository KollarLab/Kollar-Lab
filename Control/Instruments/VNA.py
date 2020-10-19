import numpy
from time import sleep, time

import pyvisa

class VNA():
    '''
    '''
    def __init__(self, address):
        '''
        '''
        rm = pyvisa.ResourceManager()
        self.inst = rm.open_resource(address)
        self.inst.write('*RST; *CLS')
        self.inst.write('OUTP OFF')

    @property
    def error(self):
        return (self.inst.query('SYST:ERR?')).strip().split(',')
    
    @property
    def output(self):
        state = int(self.inst.query('OUTP?'))
        if state:
            return 'ON'
        return 'OFF'
    @output.setter
    def output(self, state):
        self.inst.write('OUTP {}'.format(state))
    
    @property
    def power(self):
        return float(self.inst.query('SOUR:POW?'))
    @power.setter
    def power(self, power):
        self.inst.write('SOUR:POW {}'.format(power))
        
    def getErrors(self):
        '''
        '''
        code, *string = self.error
        while int(code) != 0:
            print('{},{}'.format(code, string))
            code, *string = self.error
        print('{},{}'.format(code, string))
        self.inst.write('*CLS')
            
    def SpecDefaultSettings(self):
        '''
        '''
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
    
    def SpecMeas(self, settings):
        '''
        '''
        channel      = settings['channel']
        averages     = settings['averages']
        measurement  = settings['measurement']
        meas_format  = settings['meas_form']
        start        = settings['start']
        stop         = settings['stop']
        sweep_points = settings['sweep_points']
        RFpower      = settings['RFpower']
        RFport       = settings['RFport']
        Mport        = settings['Mport']
        CAVport      = settings['CAVport']
        CAVpower     = settings['CAVpower']
        CAVfreq      = settings['CAVfreq']
        
        #Turn off output and switch to single sweep mode instead of continuous sweeps
        self.output = 'OFF'
        self.inst.write('INIT:CONT OFF')
        
        #Configure averaging
        self.configureAverages(channel, averages)
        
        #Clear old traces on channels and define a new trace to measure 'S23'
        self.configureTrace(channel, 'spec', measurement, meas_format)
        
        #Configure frequency sweep and RF power 
        self.configureFrequency(channel, start, stop, sweep_points)
        self.power = RFpower
        
        #Configure ports
        #RF
        self.inst.write('SOUR{}:POW{}:PERM ON'.format(channel, RFport)) #configures Ch1 to always be on during measurements
        #Measure
        self.inst.write('SOUR{}:POW{}:STATE OFF'.format(channel, Mport)) #configure to only receive 
        self.inst.write('SOUR{}:FREQ{}:CONV:ARB:IFR 1, 1, {}, FIX'.format(channel, Mport, CAVfreq)) #only look at what is happening at 60MHz
        #LO
        self.inst.write('SOUR{}:POW{}:PERM ON'.format(channel, CAVport)) #always keep ch3 on during measurements
        self.inst.write('SOUR{}:FREQ{}:CONV:ARB:IFR 1, 1, {}, FIX'.format(channel, CAVport, CAVfreq)) #only output at 60MHz
        self.inst.write('SOUR{}:POW{}:OFFS {}, ONLY'.format(channel, CAVport, CAVpower)) #fixed power of -5dBm (ignores the power level of Ch1)    
        self.inst.query('*OPC?')
        
        self.getErrors()
        
        self.wait_complete(channel, averages)
                
        self.autoscale(window=1, ref_trace='spec')
    
    def wait_complete(self, channel, averages):
        '''
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

    def TransDefaultSettings(self):
        '''
        '''
        settings = {}

        settings['channel'] = 1
        settings['averages'] = 100
        settings['measurement'] = 'S21'
        settings['start'] = '1 MHz'
        settings['stop'] = '40 MHz'
        settings['sweep_points'] = 501
        settings['RFpower'] = -10

        return settings
    
    def TransMeas(self, settings):
        '''
        '''
        channel      = settings['channel']
        averages     = settings['averages']
        measurement  = settings['measurement']
        start        = settings['start']
        stop         = settings['stop']
        sweep_points = settings['sweep_points']
        RFpower      = settings['RFpower']

        #Turn off output and switch to single sweep mode instead of continuous sweeps
        self.output = 'OFF'
        self.inst.write('INIT:CONT OFF')
        
        #Configure averaging
        self.configureAverages(channel, averages)
        
        #Clear old traces on channels and define a new trace to measure 'S21'
        self.configureTrace(channel, 'mag', measurement, 'MLOG')
        self.configureTrace(channel, 'phase', measurement, 'UPH')
        self.displayTrace(trace= 'mag', tracenum= 1, window= 1)
        self.displayTrace(trace= 'phase', tracenum= 2, window= 1)
        
        #Configure frequency sweep and RF power 
        self.configureFrequency(channel, start, stop, sweep_points)
        self.power = RFpower
        
        self.getErrors()
        
        #Query instrument until operation is complete
        self.wait_complete(channel, averages)
        
        self.autoscale(window= 1, ref_trace= 'mag')
        self.autoscale(window= 1, ref_trace= 'phase')
        
        mag   = self.getTrace(channel, 'mag')
        phase = self.getTrace(channel, 'phase')
        freqs = self.getChannelAxis(channel)
        
        return mag, phase, freqs
    
    def autoscale(self, window= 1, ref_trace= None):
        '''
        '''
        tracestr = ''
        if ref_trace is not None:
            tracestr = ",'{}'".format(ref_trace)
        self.inst.write("DISP:WIND{}:TRAC:Y:AUTO ONCE{}".format(window, tracestr))
        
    def configureAverages(self, channel, averages):
        '''
        '''
        self.inst.write('SENS{}:AVER ON'.format(channel))
        self.inst.write('SENS{}:AVER:COUN {}'.format(channel, averages)) 
        self.inst.write('SENS{}:AVER:CLE'.format(channel)) 
        self.inst.write('SENS{}:SWE:COUN {}'.format(channel, averages)) 
        self.inst.query('*OPC?')
      
    def configureFrequency(self, channel, start, stop, sweep_points):
        '''
        '''
        self.inst.write('SENS{}:FREQ:STAR {}; STOP {}'.format(channel, start, stop))
        self.inst.write('SENS{}:SWE:POIN {}'.format(channel, sweep_points))
        self.inst.query('*OPC?')
        
    def clearAllTraces(self):
        '''
        '''
        self.inst.write('CALC:PAR:DEL:ALL')
        self.inst.query('*OPC?')
    
    def clearChannelTraces(self, channel):
        '''
        '''
        self.inst.write('CALC{}:PAR:DEL:CALL'.format(channel))
        self.inst.query('*OPC?')
        
    def displayTrace(self, trace, tracenum, window):
        '''
        '''
        self.inst.write('DISP:WIND{}:STAT ON'.format(window))
        self.inst.query('*OPC?')
        self.inst.write('DISP:WIND{}:TRAC{}:FEED "{}"'.format(window, tracenum, trace))
        self.inst.query('*OPC?')
        
    def configureTrace(self, channel, name, meastype, measformat):
        '''
        '''
        self.inst.write("CALC{}:PAR:SDEF '{}', '{}'".format(channel, name, meastype))
        self.inst.write("CALC{}:FORM {}".format(channel, measformat))
        self.inst.query('*OPC?')
    
    def getChannelTraces(self, channel):
        '''
        '''
        instresp  = self.inst.query('CALC{}:PAR:CAT?'.format(channel))
        tracelist = instresp.strip("\n'").split(',')
        return tracelist
    
    def getChannelAxis(self, channel):
        '''
        '''
        xaxis = self.inst.query_ascii_values("CALC{}:DATA:STIM?".format(channel))
        return numpy.asarray(xaxis)
    
    def getTrace(self, channel, name):
        '''
        '''
        tracelist = self.getChannelTraces(channel)[::2]
        if name not in tracelist:
            print('Trace not found, possible options are: {}'.format(tracelist))
            return numpy.zeros(1000)
        data = self.inst.query_ascii_values("CALC{}:DATA:TRAC? '{}',FDAT".format(channel, name))
        return numpy.asarray(data)
        
    def close(self):
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