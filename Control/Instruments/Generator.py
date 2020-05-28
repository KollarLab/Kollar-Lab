import pyvisa
from userfuncs import freeze

@freeze
class Generator():

    def __init__(self, address):
        rm=pyvisa.ResourceManager()
        self.inst=rm.open_resource(address)
        self.inst.write('*RST')

    @property
    
