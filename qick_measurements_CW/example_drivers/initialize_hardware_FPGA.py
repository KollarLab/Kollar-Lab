#pip install Pyro4
import Pyro4

# from qick import QickSoc
from qick.qick_asm import QickConfig
#import pyvisa
#from kollar_instruments.SGS import SGS100A

#rm = pyvisa.ResourceManager()
#rm.close()

#logen = SGS100A('TCPIP0::rssgs100a110738::inst0::INSTR')

#logen.Ref.Source = 'INT'
#logen.Ref.Frequency = 10e6

#extgen = SGS100A('TCPIP0::rssgs100a110425::inst0::INSTR')
#
#extgen.Ref.Source = 'INT'
#extgen.Ref.Frequency = 10e6


Pyro4.config.SERIALIZER="pickle"
Pyro4.config.PICKLE_PROTOCOL_VERSION=4

#ns = Pyro4.locateNS(host="192.168.10.132", port=8888)
ns = Pyro4.locateNS(host="192.168.10.170", port=8888)
soc = Pyro4.Proxy(ns.lookup("myqick"))


soccfg = QickConfig(soc.get_cfg())


print(soccfg)

# we have a loopback from DAC 0 to ADC 1