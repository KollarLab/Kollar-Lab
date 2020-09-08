import pyvisa
from userfuncs import freeze

class RFsource(object):

    def __init__(self, address, commands):
        rm = pyvisa.ResourceManager()
        self.inst = rm.open_resource(address)
        self.inst.write('*RST')
        self.commandset = commands
    
    def __getattr__(self, name):
        cmds = self.__dict__['commandset']
        if name in cmds.keys():
            return self.inst.query(cmds[name])
        else:
            print('Setting not defined: {}'.format(name))
    
    def __setattr__(self, name, value):
        cmds = self.__dict__['commandset']
        if name in cmds.keys():
            self.inst.write(cmds['self'])
        else:
            super(RFsource, self).__setattr__(name, value)