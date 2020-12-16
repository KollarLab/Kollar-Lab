import pyvisa

class SCPIinst(object):

    def __init__(self, address, commands, errcmd, reset = True):
        rm = pyvisa.ResourceManager()
        self.inst = rm.open_resource(address)
        if reset:
            self.inst.write('*RST; *CLS')
        self.errcmd = errcmd
        self.modules = []
        for key in commands.keys():
            if key == 'core':
                self.commandset = commands['core']
            elif isinstance(commands[key],dict):
                setattr(self,key,Module(self.inst, commands[key], errcmd))
                self.modules.append(key)
    
    def __getattr__(self, name):
        cmds = self.__dict__['commandset']
        if name in cmds.keys():
            if isinstance(cmds[name], list):
                val = self.inst.query(cmds[name][0]+'?').rstrip()
                try:
                    val = eval(val)
                    strings = cmds[name][1]
                    return strings[val]
                except:
                    return val
            else:
                val = self.inst.query(cmds[name]+'?').rstrip()
                try:
                    return eval(val)
                except:
                    return val
        else:
            super().__getattribute__(self, name)
    
    def __setattr__(self, name, value):
        cmds = {}
        try:
            cmds = self.__dict__['commandset']
        except:
            pass 
        
        if name in cmds.keys():
            setting = cmds[name]
            if isinstance(setting, list):
                command = setting[0]
                try:
                    value = setting[1].inverse[value]
                except:
                    print('Invalid input, valid options are: {}'.format(list(setting[1].inverse.keys())))
            else:
                command = setting
            cmd = '{} {}'.format(command, value)
            self.inst.write(cmd)
            self.inst.write('*OPC')
            self.error
        else:
            super().__setattr__(name, value)
    
    def close(self):
        self.inst.close()
        
    @property
    def error(self):
        errcmd = self.errcmd
        if isinstance(errcmd, dict):
            for errtype in errcmd.keys():
                code = 0
                string = {}
                error_struct = errcmd[errtype]
                if isinstance(error_struct, list):
                    errorstring = self.inst.query(errcmd[errtype][0])
                    errstr = errcmd[errtype][1]
                else:
                    errorstring = self.inst.query(errcmd[errtype])
                try:
                    (code, string) = eval(errorstring)
                except:
                    code = eval(errorstring)
                    string = errstr[code]
                if code != 0:
                    print('{}:{},{}'.format(errtype, code, string))
        else:
            (code, string) = eval(self.inst.query(errcmd))
            if code != 0:
                print('{},{}'.format(code, string))
        
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

    def __init__(self,inst, commands, errcmd):
        self.inst = inst
        self.commandset = commands
        self.errcmd = errcmd
        
    def __getattr__(self, name):
        cmds = self.__dict__['commandset']
        if name in cmds.keys():
            if isinstance(cmds[name], list):
                val = self.inst.query(cmds[name][0]+'?').rstrip()
                try:
                    val = eval(val)
                    strings = cmds[name][1]
                    return strings[val]
                except:
                    return val
            else:
                val = self.inst.query(cmds[name]+'?').rstrip()
                try:
                    return eval(val)
                except:
                    return val
        else:
            super().__getattribute__(self, name)
    
    def __setattr__(self, name, value):
        cmds = {}
        try:
            cmds = self.__dict__['commandset']
        except:
            pass 
        
        if name in cmds.keys():
            setting = cmds[name]
            if isinstance(setting, list):
                command = setting[0]
                try:
                    value = setting[1].inverse[value]
                except:
                    print('Invalid input, valid options are: {}'.format(list(setting[1].inverse.keys())))
            else:
                command = setting
            cmd = '{} {}'.format(command, value)
            self.inst.write(cmd)
            self.inst.write('*OPC')
            self.error
        else:
            super().__setattr__(name, value)
            
    @property
    def error(self):
        errcmd = self.errcmd
        if isinstance(errcmd, dict):
            for errtype in errcmd.keys():
                code = 0
                string = {}
                error_struct = errcmd[errtype]
                if isinstance(error_struct, list):
                    errorstring = self.inst.query(errcmd[errtype][0])
                    errstr = errcmd[errtype][1]
                else:
                    errorstring = self.inst.query(errcmd[errtype])
                try:
                    (code, string) = eval(errorstring)
                except:
                    code = eval(errorstring)
                    string = errstr[code]
                if code != 0:
                    print('{}:{},{}'.format(errtype, code, string))
        else:
            (code, string) = eval(self.inst.query(errcmd))
            if code != 0:
                print('{},{}'.format(code, string))
                
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