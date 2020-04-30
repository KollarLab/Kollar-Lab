#import visa
import pyvisa

class RFgen():

    def __init__(self, address):
#        rm=visa.ResourceManager()
        rm=pyvisa.ResourceManager()
        self.inst=rm.open_resource(address)
        self.inst.write('*RST')

    def set_Amp_V(self,amp):
        self.inst.write('SOURce:POWer {} V'.format(amp))
        
    def set_Amp(self,amp):
        #set power in dB
        self.inst.write('SOURce:POWer {}'.format(amp))
        
        
        

    def set_Freq(self,freq):
        self.inst.write('SOURce:FREQuency:CW {} GHz'.format(freq))
        
        
    def set_External_Reference(self):
        self.inst.write(':SOURce:ROSCillator:SOURce EXTernal')
        
    def set_Internal_Reference(self):
        self.inst.write(':SOURce:ROSCillator:SOURce Internal')
        
        
    def set_Phase(self, phase):
        #degrees
        self.inst.write(':SOURce:PHASe {}'.format(phase))
        
        
    def power_On(self):
        self.inst.write(':OUTPut:STATe ON')
        
    def power_Off(self):
        self.inst.write(':OUTPut:STATe OFF')
        
        
    def mod_On(self):
        self.inst.write(':SOURce:IQ:STATe ON')
        
    def mod_Off(self):
        self.inst.write(':SOURce:IQ:STATe OFF')
        
        

    def send_cmd(self,cmd):
        self.inst.write(cmd)

    def query(self,cmd):
        print(self.inst.query(cmd))

#    def run(self):
#        self.inst.write('OUTPut:STATe ON')