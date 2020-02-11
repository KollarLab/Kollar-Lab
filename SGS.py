import visa

class RFgen():

    def __init__(self, address):
        rm=visa.ResourceManager()
        self.inst=rm.open_resource(address)
        self.inst.write('*RST')

    def set_Amp(self,amp):
        self.inst.write('SOURce:POWer {} V'.format(amp))

    def set_Freq(self,freq):
        self.inst.write('SOURce:FREQuency:CW {} GHz'.format(freq))

    def send_cmd(self,cmd):
        self.inst.write(cmd)

    def query(self,cmd):
        print(self.inst.query(cmd))

    def run(self):
        self.inst.write('OUTPut:STATe ON')