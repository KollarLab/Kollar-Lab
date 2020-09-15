#import visa
import pyvisa

class RFgen():

    def __init__(self, address):
#        rm=visa.ResourceManager()
        rm=pyvisa.ResourceManager()
        self.inst=rm.open_resource(address)
        self.inst.write('*RST')
        
    def set_amp_v(self,amp):
        self.inst.write('SOURce:POWer {} V'.format(amp))
        
    def set_amp(self,amp):
        """Set power in dB

        Arguments:
            amp {float} -- amplitude in dB
        """
        #set power in dB
        self.inst.write('SOURce:POWer {}'.format(amp))
        
    def set_freq(self,freq):
        """Set the frequency in GHz

        Arguments:
            freq {float} -- frequency of the CW in GHz 
        """
        self.inst.write('SOURce:FREQuency:CW {} GHz'.format(freq))
        
        
    def set_external_reference(self):
        self.inst.write(':SOURce:ROSCillator:SOURce EXTernal')
        
    def set_internal_reference(self):
        self.inst.write(':SOURce:ROSCillator:SOURce Internal')
        
        
    def set_phase(self, phase):
        #degrees
        self.inst.write(':SOURce:PHASe {}'.format(phase))
        
        
    def power_on(self):
        self.inst.write(':OUTPut:STATe ON')
        
    def power_off(self):
        self.inst.write(':OUTPut:STATe OFF')
        
        
    def mod_on(self):
        self.inst.write(':SOURce:IQ:STATe ON')
        
    def mod_off(self):
        self.inst.write(':SOURce:IQ:STATe OFF')
        
        
    def PEP(self, pep):
        self.inst.write(':SOURce:POWer:PEP?')
           
    def iqi_imparement(self):
        self.inst.write('SOURce:IQ:IMPairment:LEAKage:I')
               
    def iqq_imparement(self):
        self.inst.write('SOURce:IQ:IMPairment:LEAKage:Q')

    def iq_angle(self):
        self.inst.write('SOURce:IQ:IMPairment:QUADrature:ANGLe')
               
    def iq_ratio(self):
        self.inst.write('SOURce:IQ:IMPairment:IQRatio:MAGNitude')
               
        
    def send_cmd(self,cmd):
        self.inst.write(cmd)

    def query(self,cmd):
        print(self.inst.query(cmd))
        
        

#    def run(self):
#        self.inst.write('OUTPut:STATe ON')
