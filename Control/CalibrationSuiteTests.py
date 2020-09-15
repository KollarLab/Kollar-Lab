import numpy
import pylab 

from Calibration import cwhomodynestability as cwhomodyne
from Calibration import mixerIQcal as mixer
from Calibration import pulsedhomodynestability as pulsedhomodyne
from Calibration import cwheterodyne
from Calibration import pulsedheterodynestability as pulsedheterodyne

import userfuncs as uf

instruments = {}
instruments['AWG'] = hdawg
instruments['Digitizer'] = card
instruments['LO'] = logen
instruments['RFgen'] = rfgen


defaults_mixer = mixer.GetDefaultSettings()
defaults_mixer['showFig'] = True
mixer.calibrate_mixer_IQ(instruments, defaults_mixer)

settings = cwhomodyne.GetDefaultSettings()
settings['measure_time'] = 35
settings['channelRange'] = 2.5
settings['save'] = False
cwhomodyne.CWHomodyneStabilityTest(instruments, settings)

defaults_pulsed = pulsedhomodyne.GetDefaultSettings()
settings = defaults_pulsed
settings['tau_max'] = 100e-6
settings['num_points'] = 30
settings['save'] = False
pulsedhomodyne.PulsedHomodyneStability(instruments, settings)

settings = cwheterodyne.GetDefaultSettings()
settings['frequency_IF'] = 1e6
settings['measure_time'] = 30
settings['one_shot_time'] = 10e-6
settings['save'] = False
cwheterodyne.CWHeterodyneStabilityTest(instruments, settings)

settings = pulsedheterodyne.GetDefaultSettings()
settings['IF_freq'] = 1.1e6
settings['pulse_width'] = 1./settings['IF_freq']
settings['tau_max'] = 15e-6
settings['num_points']  = 17
settings['tau_min'] = 11e-6
settings['segments'] = 20
settings['save'] = False
pulsedheterodyne.PulsedHeterodyneStability(instruments, settings)

#pplots = numpy.linspace(1.,2.,10, endpoint = False)
#
#fig, axs = pylab.subplots(2,5, figsize=(15, 6), facecolor='w', edgecolor='k')
#fig.subplots_adjust(hspace = .5, wspace=1)
#axs = axs.ravel()
#
#for pind in range(len(pplots)):
#    print('Period multiple: {}'.format(pplots[pind]))
#    settings = pulsedheterodyne.GetDefaultSettings()
#    settings['IF_freq'] = 1.1e6
#    settings['pulse_width'] = pplots[pind]/settings['IF_freq']
#    settings['tau_max'] = 15e-6
#    settings['num_points']  = 201
#    settings['tau_min'] = 11e-6
#    settings['segments'] = 20
#    settings['save'] = False
#    
#    phaseDiff,ts = pulsedheterodyne.PulsedHeterodyneStability(instruments, settings)
#    axs[pind].scatter(ts*1e3, phaseDiff, c = 'deepskyblue', s = 7)
#
#pylab.tight_layout()
#pylab.show()