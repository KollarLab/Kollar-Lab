from Instruments.HDAWG import HDAWG

hdawg = HDAWG('dev8163')
hdawg.AWGs[0].samplerate = '2.4GHz'
hdawg.channelgrouping = '1x4'
hdawg.Channels[0].configureChannel(amp=1.0,marker_out='Marker', hold='True')
hdawg.Channels[1].configureChannel(marker_out='Trigger', hold='True')
hdawg.AWGs[0].Triggers[0].configureTrigger(slope='rising',channel='Trigger in 1')
hdawg.OSCs[1].freq = 10e6
hdawg.Channels[2].analog_outs = [0.5,0]
hdawg.Channels[3].analog_outs = [0,1.0]
hdawg.Channels[2].configureChannel(amp=1.0)
hdawg.Channels[3].configureChannel(amp=2.0)

## Experimental parameters
pulseamp = 0.5
pulselength = 200e-9
maxtime = 100e-6
tau = 1e-6

## Read in the sequencer program we want to use and configure it
progFile = open("HDAWG_sequencer_codes/HomodynePulse.cpp",'r')
rawprog  = progFile.read()
loadprog = rawprog
progFile.close
loadprog = loadprog.replace('_Amp_', str(pulseamp))
loadprog = loadprog.replace('_Time_', str(pulselength))
loadprog = loadprog.replace('_max_time_', str(maxtime))
loadprog = loadprog.replace('_tau_', str(tau))

hdawg.AWGs[0].load_program(loadprog)