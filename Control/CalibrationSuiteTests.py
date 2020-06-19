from Instruments.HDAWG import HDAWG
#from Instruments.Acqiris import Acqiris
from Instruments.SGS import RFgen

from Calibration import cwhomodynestability as cwhomodyne
from Calibration import mixerIQcal as mixer
from Calibration import pulsedhomodynestability as pulsedhomodyne
from Calibration import cwheterodyne

import userfuncs as uf

instruments = {}
instruments['AWG'] = hdawg
instruments['Digitizer'] = card
instruments['LO'] = logen
instruments['RFgen'] = rfgen


#defaults_mixer = mixer.GetDefaultSettings()
#defaults_mixer['showFig'] = True

#settings = cwhomodyne.GetDefaultSettings()
#settings['measure_time'] = 36000

#defaults_pulsed = pulsedhomodyne.GetDefaultSettings()
#settings = defaults_pulsed
#settings['tau_max'] = 1000e-6
#settings['num_points'] = 99
#settings['save'] = False

settings = cwheterodyne.GetDefaultSettings()
settings['frequency_IF'] = 1e6
settings['measure_time'] = 1000
settings['one_shot_time'] = 10e-6
## Example usage:
#mixer.calibrate_mixer_IQ(instruments, defaults_mixer)
#cwhomodyne.CWHomodyneStabilityTest(instruments, settings)
#pulsedhomodyne.pulsedhomodynestability(instruments, settings)
cwheterodyne.CWHeterodyneStabilityTest(instruments, settings)
