import time


saveDir = r'Z:\Data\HouckTaTransmon\20201111'
name = 'PulseSpec'
stamp = userfuncs.timestamp()
scanname = name + '_' + stamp
    
settings = {}

settings['CAVpower']        = -15
settings['CAV_freq']        = 7.173145e9

freqs = numpy.linspace(4.9, 5.2, 101)*1e9
#powers = numpy.linspace(-10, 5, 5)
powers = [20]

CAV_freq  = settings['CAV_freq']
CAVpower  = settings['CAVpower'] 

## Generator settings
cavitygen.Freq   = CAV_freq
cavitygen.Power  = CAVpower
cavitygen.IQ.Mod = 'On'

qbitgen.Freq   = 4e9
qbitgen.Power  = -20
qbitgen.IQ.Mod = 'On'

cavitygen.Output = 'On'
qbitgen.Output = 'On'

#Card settings
settings['segments']         = 1
settings['reads']            = 1
settings['averages']         = 1e3
settings['activeChannels']   = [1,2]
settings['sampleRate']       = 2e9/8
settings['trigger_buffer']   = 0e-6

card.averages       = settings['averages']
card.segments       = settings['segments']
card.sampleRate     = settings['sampleRate']
card.activeChannels = settings['activeChannels']
card.triggerDelay   = settings['trigger_buffer']
card.timeout        = 30
card.samples        = 400e2
card.channelRange   = 1

xaxis = (numpy.array(range(card.samples))/card.sampleRate)*1e6

powerdat = []

t1 = time.time()

for powerind in range(len(powers)):
    qbitgen.Power = powers[powerind]
    time.sleep(0.2)
    
    amps = []
    phases = []
    
    Is = []
    Qs = []
    
    print('Current power:{}, max:{}'.format(powers[powerind], powers[-1]))
    for freq in freqs:
        qbitgen.Freq = freq
        time.sleep(0.2)
        
        card.ArmAndWait()
        
        I,Q = card.ReadAllData()
        
        DC_I = numpy.mean(I[0][-50:])
        DC_Q = numpy.mean(Q[0][-50:])
        
        Idat = I[0]-DC_I
        Qdat = Q[0]-DC_Q
        
        amp = numpy.sqrt(Idat**2+Qdat**2)
        phase = numpy.arctan2(Qdat, Idat)*180/numpy.pi
        
        amps.append(amp)
        phases.append(phase)
        Is.append(Idat)
        Qs.append(Qdat)
    
    ampsarray = numpy.asarray(amps)
    powerslice = numpy.mean(ampsarray[:,0:int(len(amp)/2)], axis=1)/(10**(power/20))
    powerdat.append(powerslice)
    
    userfuncs.SaveFull(saveDir, scanname, ['powers','freqs', 'powerdat','xaxis'], locals(), expsettings=settings)

t2 = time.time()



###plot the full power dependent data 
#fig = plt.figure(1)
#ax = plt.subplot(1,1,1)
#base_power_plot(fig, ax, freqs, powerdat, powers, 'Freq sweep', 'mag', HWattenuation = -30)  
#plt.savefig(os.path.join(saveDir, scanname+'_fullColorPlot.png'), dpi = 150)

#plot the raw data at (the last?) power
fig = plt.figure(2)
ax = plt.subplot(1,1,1)
base_power_plot_spec(fig, ax, xaxis, amps, freqs, 'Raw traces','',0)
#plt.savefig(os.path.join(saveDir, scanname+'_singleRawData.png'), dpi = 150)


#plotthe last row as a trace.
fig = plt.figure(3)
ax = plt.subplot(1,1,1)
plt.plot(freqs/1e9, powerslice)
plt.xlabel('freqs (GHz)')
plt.ylabel('spec voltage')
#plt.savefig(os.path.join(saveDir, scanname+'rowTrace.png'), dpi = 150)


















