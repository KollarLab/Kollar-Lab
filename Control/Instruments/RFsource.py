import pyvisa
from userfuncs import freeze

class RFsource(object):

    def __init__(self, address, commands):
        rm = pyvisa.ResourceManager()
        self.inst = rm.open_resource(address)
        self.inst.write('*RST')
        self.commandset = commands
    
    def __getattribute__(self, name):
        if name in self.commandset.keys():
            return self.inst.query(self.commandset[name])
        else:
            print('Setting not defined: {}'.format(name))
    
    def __setattr__(self, name, value):
        if name in self.commandset.keys():
            self.inst.write(self.commandset[name])