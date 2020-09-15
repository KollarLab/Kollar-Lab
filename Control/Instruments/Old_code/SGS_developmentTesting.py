#import visa
import pyvisa

class SGSgen():

    def __init__(self, address):
        ''' version of an RFgen object for SGS generator
        
        syntax:
        rfgen = SGSgen(hardwareaddress)
        rfgen.Power_On()
        
        or
        rfgen = SGS.RFgen(hardwareaddress)
        rfgen.PowerOn()
        
        different make and model would look like
        rfgen = Holzworth(hardwareaddress)
        
        or
        rfgen = Holzworth.RFgen(hardwareaddress)
        rfgen.PowerOn()
        
        
        
        #properties
        ###########
        %settings
        %hardware_adress
        %freq
        %freq_GHz
        %power
        %?amp_V
        %phase 
        
        %reference_source (internal v external)
        %reference_in_freq
        reference_out_freq
        send_reference (boolean)
        
        send_LO (boolean) (send_reference, and send_LO cannot both be true)
        receive_LO (boolean)
        
        %mod_enabled (boolean)
        
        #IQsettings
        I_offset
        Q_offset
        IQ_imparement
        IQ_angle
        IQ_ratio
        
        would not be crazy to have
        rfgen.IQ.ratio
        v 
        rfgen.IQ_ratio
        
        #methods
        ##########
        %power_on
        %power_off
        %mod_on, or maybe enable_mod
        %mod_off, or maybe disable_mod
        
        1/2%configure_IQ
        
        1/2%configure_mod
        
        configure_ref_out
        %configure_ref_in
        
        configure_LO_out
        configure_LO_in
        
        %send_cmd
        %querry
        
        
        #housekeeping
        ############
        load_config #turn a set of settings into a functioning object
        _loadDefaultConfig #store a default internally
        print (a bit redundant with settings)
                
        
        '''
#        rm=visa.ResourceManager()
        rm=pyvisa.ResourceManager()
        self.inst=rm.open_resource(address)
        self.inst.write('*RST')

    #property decorators
    ####################
    @property
    def settings(self):
        fullsettings = {}
        for setting in dir(self):
            if setting == 'inst' or setting == 'settings':
                continue
            if not setting.startswith('__') and not callable(getattr(self,setting)):
                fullsettings[setting] = getattr(self,setting)
        return fullsettings
    @settings.setter
    def settings(self, fullsettings):
        for key in fullsettings.keys():
            try:
                setattr(self, key, fullsettings[key])
            except:
                print('Unknown attribute {}, ignoring'.format(key))
                
                
    @property
    def frequency(self):
        freq = self.get_Freq()
        return freq
    @frequency.setter
    def frequency(self, freq):
        self.set_Freq(freq)
        
    @property
    def power(self):
        '''Doc string '''
        power = self.get_Power()
        return power
    @power.setter
    def power(self, poww):
        self.set_Power(poww)    
        
                
    #advanced configuration functions
    #############            
                
                
    #low-level get, set functions
    ###################            

    def set_Amp_V(self,amp):
        self.inst.write('SOURce:POWer {} V'.format(amp))
        
    def set_Power(self,amp):
        #set power in dB
        self.inst.write('SOURce:POWer {}'.format(amp))
    def get_Power(self):
        #set power in dB
        temp = self.inst.query_ascii_values('SOURce:POWer?')
        poww = temp[0]
        return poww

    def set_Freq(self,freq):
        freq_GHz = freq/1e9
        self.inst.write('SOURce:FREQuency:CW {} GHz'.format(freq_GHz))      
    def get_Freq(self):    
        temp = self.inst.query_ascii_values('SOURce:FREQuency:CW?')
        freq = temp[0]
        return freq
        
    def set_Offset(self, offset):
        self.inst.write('SOURce:FREQuency:OFFSet {} MHz'.format(offset))
        
    def set_External_Reference(self, freq=10):
        self.inst.write(':SOURce:ROSCillator:SOURce EXTernal')
        self.inst.write(':SOURce:ROSCillator:EXTernal:FREQuency {} MHz'.format(freq))
        
    def set_Internal_Reference(self):
        self.inst.write(':SOURce:ROSCillator:SOURce Internal')
        
    def set_External_LO(self):
        self.inst.write(':SOURce:LOSCillator:SOURce EXTernal')
    
    def set_Internal_LO(self):
        self.inst.write(':SOURce:LOSCillator:SOURce INTernal')
    
    def set_RefLO_output(self, output = 'Ref', freq = 10):
        self.inst.write(':CONNector:REFLo:OUTPut {}'.format(output))
        if output == 'Ref':
            self.inst.write(':SOURce:ROSCillator:OUTPut:FREQuency {} MHz'.format(freq))
        
    def set_Phase(self, phase):
        #degrees
        self.inst.write(':SOURce:PHASe {}'.format(phase))
        
    
    #low-level configuration functions
    ############################
    
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
        
        
if __name__ == '__main__':

    hardwareAddress = 'TCPIP0::rssgs100a110738::inst0::INSTR'

    testgen = SGSgen(hardwareAddress)

    
        
        
        
        
        
        