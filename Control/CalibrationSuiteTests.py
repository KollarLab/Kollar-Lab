from Instruments.HDAWG import HDAWG
#from Instruments.Acqiris import Acqiris
from Instruments.SGS import RFgen

from Calibration import cwhomodynestability as cwhomodyne
from Calibration import mixerIQcal as mixer
from Calibration import PulsedHomodyneStability as pulsedhomodyne

import userfuncs as uf

#card  = Acqiris('PXI23::0::0::INSTR')
#hdawg = HDAWG('dev8163')
#logen = RFgen('TCPIP0::rssgs100a110738::inst0::INSTR')
#rfgen = RFgen('TCPIP0::rssgs100a110739::inst0::INSTR')
#
instruments = {}
instruments['AWG'] = hdawg
instruments['Digitizer'] = card
instruments['LO'] = logen
instruments['RFgen'] = rfgen

defaults_pulsed = pulsedhomodyne.GetDefaultSettings()
settings = defaults_pulsed
settings['tau_max'] = 100e-6
settings['num_points'] = 10
settings['save'] = True

## Example usage:
#mixer.calibrate_mixer_IQ(instruments, defaults_mixer)
#cwhomodyne.CWHomodyneStabilityTest(instruments, settings)
pulsedhomodyne.pulsedhomodynestability(instruments, settings)
