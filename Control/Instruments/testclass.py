class test(object):
    names = {}
    names['prop'] = 'wer'
    names['number'] = 5
    names['list'] = ['asd','qwe','zxc']
    names['subobj'] = 'asd'

    def __init__(self, prop):
        self.prop = prop
        self.subobj = testsub('def')
        print('Initialized')
    def __getattr__(self,name):
        print('Called getattribute1')
        if name in self.names.keys():
            print(name, object.__getattribute__(self,name))
        else:
            print('No attribute of that name exists')
    def __setattr__(self, name, value):
        print('Called setattribute1')
        if name in self.names.keys():
            object.__setattr__(self, name, value)
        else:
            return

class testsub(object):
    names ={}
    names['subname'] = 'sdf'

    def __init__(self,subname):
        self.subname = subname
    def __getattr__(self,name):
        print('Called getattribute2')
        if name in self.names.keys():
            print(name, object.__getattribute__(self,name))
        else:
            print('No attribute of that name exists')
    def __setattr__(self, name, value):
        print('Called setattribute2')
        if name in self.names.keys():
            object.__setattr__(self, name, value)
        else:
            return