import pyvisa

class SCPIinst(object):

    def __init__(self, address, commands):
        rm = pyvisa.ResourceManager()
        self.inst = rm.open_resource(address)
        self.inst.write('*RST')
        self.modules = []
        for key in commands.keys():
            if key == 'core':
                self.commandset = commands['core']
            elif isinstance(commands[key],dict):
                setattr(self,key,Module(self.inst, commands[key]))
                self.modules.append(key)
    
    def __getattr__(self, name):
        cmds = self.__dict__['commandset']
        if name in cmds.keys():
            val = self.inst.query(cmds[name]+'?').rstrip()
            try:
                return eval(val)
            except:
                return val
        else:
            super().__getattribute__(self, name)
    
    def __setattr__(self, name, value):
#        print('Called setattribute1')
        cmds = {}
        try:
            cmds = self.__dict__['commandset']
        except:
            pass 
        
        if name in cmds.keys():
            cmd = '{} {}'.format(cmds[name], value)
#            print(cmd)
            self.inst.write(cmd)
            self.inst.write('*OPC?')
        else:
#            print('Doin stupid')
            super().__setattr__(name, value)
    @property
    def settings(self):
        fullsettings = {}
        for setting in self.commandset:
            fullsettings[setting] = self.__getattr__(setting)
        for module in self.modules:
            fullsettings[module] = {}
            mod = getattr(self, module)
            fullsettings[module] = mod.settings
        return fullsettings
    @settings.setter
    def settings(self, fullsettings):
        for setting in fullsettings.keys():
            if not isinstance(fullsettings[setting], dict):
                self.__setattr__(setting, fullsettings[setting])
            else:
                try:
                    mod = getattr(self, setting)
                    mod.settings = fullsettings[setting]
                except:
                    print('invalid')
            
class Module(object):

    def __init__(self,inst, commands):
        self.inst = inst
        self.commandset = commands
        
    def __getattr__(self, name):
        cmds = self.__dict__['commandset']
        if name in cmds.keys():
            val = self.inst.query(cmds[name]+'?').rstrip()
            try:
                return eval(val)
            except:
                return val
        else:
            super().__getattribute__(self, name)
    
    def __setattr__(self, name, value):
#        print('Called setattribute1')
        cmds = {}
        try:
            cmds = self.__dict__['commandset']
        except:
            pass 
        
        if name in cmds.keys():
            cmd = '{} {}'.format(cmds[name], value)
#            print(cmd)
            self.inst.write(cmd)
            self.inst.write('*OPC?')
        else:
#            print('Doin stupid')
            super().__setattr__(name, value)
    @property
    def settings(self):
        fullsettings = {}
        for setting in self.commandset:
            fullsettings[setting] = self.__getattr__(setting)
#        for module in modules:
#            fullsettings[module] = {}
#            mod = getattr(self, module)
#            fullsettings[module] = mod.settings
        return fullsettings
    @settings.setter
    def settings(self, fullsettings):
        for setting in fullsettings.keys():
            self.__setattr__(setting, fullsettings[setting])