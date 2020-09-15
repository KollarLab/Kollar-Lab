import pyvisa
import pylab
from time import sleep, time


class VNA():
    def __init__(self, address):
        rm = pyvisa.ResourceManager()
        self.inst = rm.open_resource(address)
        self.inst.write('*RST')
    
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