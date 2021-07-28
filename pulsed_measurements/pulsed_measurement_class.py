from userfuncs import saveDir
from utility.measurement_helpers import configure_card, read_and_process, estimate_time, check_inputs

class pulsed_measurement():
    def __init__(self, instruments, settings):

        self.exp_settings = settings['exp_settings']
        self.exp_globals  = settings['exp_globals']

        self.qubitgen  = instruments['qubitgen']
        self.cavitygen = instruments['cavitygen']
        self.card      = instruments['card']
        self.hdawg     = instruments['hdawg']
        self.LO        = instruments['LO']

        self.saveDir = saveDir(settings)

        self.hdawg_path = r"C:\Users\Kollarlab\Desktop\Kollar-Lab\pulsed_measurements\HDAWG_sequencer_codes"

    def configure_hdawg(self):
        self.hdawg.AWGs[0].samplerate = '2.4GHz'
        self.hdawg.channelgrouping = '1x4'
        self.Channels[0].configureChannel(amp=1.0, marker_out='Marker', hold='False')
        self.Channels[1].configureChannel(amp=1.0, marker_out='Marker', hold='False')
        self.AWGs[0].Triggers[0].configureTrigger(slope='rising', channel='Trigger in 1')



