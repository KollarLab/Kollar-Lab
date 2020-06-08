from Instruments.HDAWG import HDAWG
#from Instruments.Acqiris import Acqiris
from Instruments.SGS import RFgen

from calibration import cwhomodynestability as cwhomodyne
from calibration import mixerIQcal as mixer
from calibration import PulsedHomodyneStability as pulsedhomodyne

import userfuncs as uf

#card  = Acqiris('PXI23::0::0::INSTR')
#hdawg = HDAWG('dev8163')
#logen = RFgen('TCPIP0::rssgs100a110738::inst0::INSTR')
#rfgen = RFgen('TCPIP0::rssgs100a110739::inst0::INSTR')
#
#instruments = {}
#instruments['AWG'] = hdawg
#instruments['Digitizer'] = card
#instruments['LO'] = logen
#instruments['RFgen'] = rfgen

defaults_mixer  = mixer.GetDefaultSettings()
defaults_cw     = cwhomodyne.GetDefaultSettings()
defaults_pulsed = pulsedhomodyne.GetDefaultSettings()

print('Mixer defaults: {}'.format(defaults_mixer))
print('CW Homodyne defaults: {}'.format(defaults_cw))
print('Pulsed Homodyne defaults: {}'.format(defaults_pulsed))

## Example usage:
#mixer.calibrate_mixer_IQ(instruments, defaults_mixer)
#cwhomodyne.CWHomodyneStabilityTest(instruments, defaults_cw)
#pulsedhomodyne.pulsedhomodynestability(instruments, defaults_pulsed)
data = []
figures = []
settings = {}
settings['Mixer'] = defaults_mixer
settings['CW'] = defaults_cw
settings['Pulsed'] = defaults_pulsed
testpath = r'C:\Users\Martin\Kollar-Lab\Control'
uf.SaveData(testpath, 'test.pkl', data, settings, figures)