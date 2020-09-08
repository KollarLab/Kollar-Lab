class test(object):
    def __init__(self, name, cmds):
        self.prop = name
        self.subobj = testsub('subObj')
        self.commandset = cmds
        print('Initialized')
        for key in cmds.keys():
            if key == 'core':
                continue
            else:
                setattr(self, key, testsub(key))
    def __getattr__(self,name):
        print('Called getattribute1')
        cmds = self.__dict__['commandset']['core']
        if name in cmds.keys():
            print(name, cmds[name])
        else:
            print('No attribute of that name exists')
    def __setattr__(self, name, value):
        print('Called setattribute1')
        try:
            cmds = self.__dict__['commandset']['core']
            if name in cmds.keys():
                print('Name: {}, value {}'.format(name, value))
            else:
                super().__setattr__(name, value)
        except:
            super().__setattr__(name, value)

class testsub(object):

    def __init__(self,subname):
        self.subname = subname
    def __getattr__(self,name):
        print('Called getattribute2')
    def __setattr__(self, name, value):
        print('Called setattribute2')
        super().__setattr__(name, value)