import pyvisa
import pylab
from time import sleep, time


class VNA():
    def __init__(self, address):
        rm = pyvisa.ResourceManager()
        self.inst = rm.open_resource(address)
        self.inst.write('*RST')
        self.inst.write('OUTP OFF')
    
    @property
    def error(self):
        return self.inst.query('SYST:ERR?')
    
    def getErrors(self):
        (code,string) = eval(self.error)
        while code != 0:
            print('{},{}'.format(code,string))
            (code,string) = eval(self.error)
        print('{},{}'.format(code,string))
            
    def configureSmeasurement(self):
        self.inst.write('INIT:CONT OFF')
        self.inst.write('CALC:PAR:DEL:ALL')
        
        self.inst.write("CALC:PAR:DEF:SGR 1,2")
        self.inst.write("CALC:FORM MLOG")
        self.inst.write("DISP:WIND1:TRAC12:FEED 'Ch1_SG_S11'")
        self.inst.write("DISP:WIND1:TRAC13:FEED 'Ch1_SG_S12'")
        self.inst.write("DISP:WIND1:TRAC14:FEED 'Ch1_SG_S21'")
        self.inst.write("DISP:WIND1:TRAC15:FEED 'Ch1_SG_S22'")
        
        self.inst.write("CALC2:PAR:DEF:SGR 1,2")
        self.inst.write("CALC2:FORM UPH")
        self.inst.write("DISP:WIND2:STAT ON")
        self.inst.write("DISP:WIND2:TRAC22:FEED 'Ch2_SG_S11'")
        self.inst.write("DISP:WIND2:TRAC23:FEED 'Ch2_SG_S12'")
        self.inst.write("DISP:WIND2:TRAC24:FEED 'Ch2_SG_S21'")
        self.inst.write("DISP:WIND2:TRAC25:FEED 'Ch2_SG_S22'")
        
        
    def getAvgTrace(self, numAvg):
        self.inst.write('SENS:AVER ON')
        self.inst.write('SENS:AVER:CLE')
        self.inst.write('SENS:AVER:COUN {}'.format(numAvg))
        self.inst.write('SENS:AVER:MODE AUTO')
        self.inst.write('SWE:COUN {}'.format(numAvg)) 
        print(self.error)
        
        self.inst.write('INIT:IMM; *OPC')
        opc = 0
        t1 = time()
        while not (opc&1):
            t2 = time()
            print('waiting for command to complete, elapsed time:{}'.format(t2-t1))
            opc = eval(self.inst.query('*ESR?').rstrip())
            sleep(1)
            
        sdat = self.inst.query_ascii_values("CALC:DATA:SGR? FDAT")
        print(self.error)
        pylab.figure()
        for i in range(4):
            l = int(len(sdat)/4)
            pylab.plot(sdat[i*l:(i+1)*l])
            
    def SpecDefaultSettings(self):
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
        self.inst.write('OUTP OFF')
        self.inst.write('INIT:CONT OFF')
        
        #Configure averaging
        self.inst.write('SENS{}:AVER ON'.format(channel))
        self.inst.write('SENS{}:AVER:COUN {}'.format(channel, averages)) #set the number of averages
        self.inst.write('SENS{}:AVER:CLE'.format(channel)) #clear old average data
        self.inst.write('SENS{}:SWE:COUN {}'.format(channel, averages)) #set the number of sweeps
        self.inst.query('*OPC?')
        
        #Clear old traces on channels and define a new trace to measure 'S23'
        self.inst.write('CALC{}:PAR:DEL:ALL'.format(channel))
        self.inst.write('CALC{}:PAR:SDEF "Ch1Tr1", "{}"'.format(channel, measurement)) #make a new trace to measure 'S23'
        self.inst.write('CALC{}:FORM {}'.format(channel, meas_format))
        self.inst.write('DISP:WIND1:TRAC1:FEED "Ch1Tr1"') #display trace in window
        self.inst.query('*OPC?')
        
        #Configure frequency sweep and RF power 
        self.inst.write('SENS{}:FREQ:STAR {}; STOP {}'.format(channel, start, stop)) #set the frequency range of the sweep
        self.inst.write('SENS{}:SWE:POIN {}'.format(channel, sweep_points)) #set the number of points in sweep
        self.inst.write('SOUR{}:POW {}'.format(channel, RFpower)) #set base output power (for RF)
        self.inst.query('*OPC?')
        
         #    Configure ports
         #    RF
        self.inst.write('SOUR{}:POW{}:PERM ON'.format(channel, RFport)) #configures Ch1 to always be on during measurements
        #    Measure
        self.inst.write('SOUR{}:POW{}:STATE OFF'.format(channel, Mport)) #configure to only receive 
        self.inst.write('SOUR{}:FREQ{}:CONV:ARB:IFR 1, 1, {}, FIX'.format(channel, Mport, CAVfreq)) #only look at what is happening at 60MHz
        #    LO
        self.inst.write('SOUR{}:POW{}:PERM ON'.format(channel, CAVport)) #always keep ch3 on during measurements
        self.inst.write('SOUR{}:FREQ{}:CONV:ARB:IFR 1, 1, {}, FIX'.format(channel, CAVport, CAVfreq)) #only output at 60MHz
        self.inst.write('SOUR{}:POW{}:OFFS {}, ONLY'.format(channel, CAVport, CAVpower)) #fixed power of -5dBm (ignores the power level of Ch1)    
        self.inst.query('*OPC?')
        
        self.getErrors()
        
        sweep_time = eval(self.inst.query('SENS{}:SWE:TIME?'.format(channel)))
        
        #    Turn on output and wait until sweep is complete to autoscale the trace
        self.inst.write('OUTP ON')
        self.inst.write('INIT{}:IMM; *OPC'.format(channel))
        
        opc = 0
        auto_scale = True
        t1 = time()
        while not (opc&1):
            t2 = time()
            print('waiting for command to complete, elapsed time:{}, total time est. {}'.format(t2-t1, 2.5*averages*sweep_time))
            opc = eval(self.inst.query('*ESR?').rstrip())
            sleep(min(max(sweep_time*0.2*averages,5),2.5*averages*sweep_time))
            if auto_scale:
                self.inst.write('DISP:WIND1:TRAC1:Y:AUTO ONCE')
                auto_scale = False
                
        self.inst.write('DISP:WIND1:TRAC1:Y:AUTO ONCE')
    
    def wait_complete(self, channel, averages):
        sweep_time = eval(self.inst.query('SENS{}:SWE:TIME?'.format(channel)))
        
        #    Turn on output and wait until sweep is complete to autoscale the trace
        self.inst.write('OUTP ON')
        self.inst.write('INIT{}:IMM; *OPC'.format(channel))
        
        opc = 0
        auto_scale = True
        t1 = time()
        while not (opc&1):
            t2 = time()
            print('waiting for command to complete, elapsed time:{}, total time est. {}'.format(t2-t1, 2.5*averages*sweep_time))
            opc = eval(self.inst.query('*ESR?').rstrip())
            sleep(min(max(sweep_time*0.2*averages,5),2.5*averages*sweep_time))
            if auto_scale:
                self.inst.write('DISP:WIND1:TRAC1:Y:AUTO ONCE')
                auto_scale = False

    def TransDefaultSettings(self):

        settings = {}

        settings['channel'] = 1
        settings['averages'] = 100
        settings['measurment'] = 'S21'
        settings['start'] = '1 MHz'
        settings['stop'] = '40 MHz'
        settings['sweep_points'] = 501
        settings['RFpower'] = -10

        return settings
    
    def TransMeas(self, settings):
        channel      = settings['channel']
        averages     = settings['averages']
        measurement  = settings['measurement']
        start        = settings['start']
        stop         = settings['stop']
        sweep_points = settings['sweep_points']
        RFpower      = settings['RFpower']

        #Turn off output and switch to single sweep mode instead of continuous sweeps
        self.inst.write('OUTP OFF')
        self.inst.write('INIT:CONT OFF')
        
        #Configure averaging
        self.inst.write('SENS{}:AVER ON'.format(channel))
        self.inst.write('SENS{}:AVER:COUN {}'.format(channel, averages)) #set the number of averages
        self.inst.write('SENS{}:AVER:CLE'.format(channel)) #clear old average data
        self.inst.write('SENS{}:SWE:COUN {}'.format(channel, averages)) #set the number of sweeps
        self.inst.query('*OPC?')
        
        #Clear old traces on channels and define a new trace to measure 'S21'
        self.inst.write('CALC{}:PAR:DEL:CALL'.format(channel))
        self.inst.write('CALC{}:PAR:SDEF "mag", "{}"'.format(channel, measurement)) #make a new trace to measure 'S21'
        self.inst.write('CALC{}:FORM {}'.format(channel, 'MLOG'))
        self.inst.write('CALC{}:PAR:SDEF "phase", "{}"'.format(channel, measurement)) #make a new trace to measure 'S21'
        self.inst.write('CALC{}:FORM {}'.format(channel, 'UPH'))
        self.inst.write('DISP:WIND1:STAT ON')
        self.inst.write('DISP:WIND2:STAT ON')
        self.inst.write('DISP:WIND1:TRAC1:FEED "mag"') #display trace in window
        self.inst.write('DISP:WIND2:TRAC2:FEED "phase"') #display trace in window
        self.inst.query('*OPC?')
        
        #Configure frequency sweep and RF power 
        self.inst.write('SENS{}:FREQ:STAR {}; STOP {}'.format(channel, start, stop)) #set the frequency range of the sweep
        self.inst.write('SENS{}:SWE:POIN {}'.format(channel, sweep_points)) #set the number of points in sweep
        self.inst.write('SOUR{}:POW {}'.format(channel, RFpower)) #set base output power (for RF)
        self.inst.query('*OPC?')

        self.wait_complete(channel, averages)
        mag = self.inst.query_ascii_values("CALC{}:DATA:TRAC 'mag',FDAT".format(channel))
        phase = self.inst.query_ascii_values("CALC{}:DATA:TRAC 'phase',FDAT".format(channel))

        return mag, phase

    
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
##        axis = inst.query("SENS2:X?") incorrect syntax
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
###Other useful stuff:
##    
##    Bandwidth:
##        inst.write('SENS1:BAND:RES 1')
##        inst.write('SENS1:BAND:RES:SEL NORM')
##    
##    Electrical delay:
##        inst.write('SENS2:CORR:EDEL:AUTO ONCE')
##        inst.write('SENS2:CORR:LOSS:AUTO ONCE') (loss and edelay)
##        inst.write('SENS2:CORR:OFFS:COMP ON') (toggle compensation)