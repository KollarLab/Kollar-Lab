from Instruments.HDAWG import HDAWG
#from Instruments.Acqiris import Acqiris
from Instruments.SGS import RFgen

from calibration import cwhomodynestability
from calibration import PulsedHomodyneStability
#card  = Acqiris('PXI23::0::0::INSTR')
#hdawg = HDAWG('dev8163')
#logen = RFgen('TCPIP0::rssgs100a110738::inst0::INSTR')
#rfgen = RFgen('TCPIP0::rssgs100a110739::inst0::INSTR')
#
instruments = {}
#instruments['AWG'] = hdawg
#instruments['Digitizer'] = card
#instruments['LO'] = logen
#instruments['RFgen'] = rfgen

defaults = PulsedHomodyneStability.GetDefaultSettings()
settings = defaults
print(settings)
#CWHomodyneStabilityTest(instruments, settings)
