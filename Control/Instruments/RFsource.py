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
            return self.inst.query(cmds[name]+'?')
        else:
            super().__getattribute__(self, name)
    
    def __setattr__(self, name, value):
        print('Called setattribute1')
        cmds = {}
        try:
            cmds = self.__dict__['commandset']
            configured = True
        except:
            print('Doin stupid')
            super().__setattr__(name, value)
        
        if name in cmds.keys() and configured:
            cmd = '{} {}'.format(cmds[name], value)
            print(cmd)
            self.inst.write(cmd)
        else:
            print('Doin stupid')
            super().__setattr__(name, value)
            
#class testsub(object):
#
#    def __init__(self,subname):
#        self.subname = subname
#    def __getattr__(self,name):
#        print('Called getattribute2')
#    def __setattr__(self, name, value):
#        print('Called setattribute2')
#        super().__setattr__(name, value)