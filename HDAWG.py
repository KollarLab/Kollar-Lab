#Set frequency for sine oscillator

import zhinst.ziPython as zi
import time

device = 'DEV8163' #unique identifier for HDAWG

discovery = zi.ziDiscovery()
discovery.find(device)
device_props = discovery.get(device)
daq = zi.ziDAQServer(device_props['serveraddress'], device_props['serverport'], device_props['apilevel'])

h = daq.awgModule()
h.set('awgModule/device', 'dev8163')
h.set('awgModule/index', 0)
h.execute()

h = daq.awgModule()
h.set('awgModule/device', 'dev8163')
h.set('awgModule/index', 1)
h.execute()

daq.setInt('/dev8163/sigouts/3/on', 0)
daq.setInt('/dev8163/sigouts/2/on', 0)
daq.setInt('/dev8163/sigouts/1/on', 0)
daq.setInt('/dev8163/sigouts/0/on', 0)

#initial setting for oscillators
daq.setDouble('/dev8163/oscs/0/freq', 500000)
daq.setDouble('/dev8163/oscs/1/freq', 1000000)

daq.setDouble('/dev8163/sines/0/phaseshift', 45)
daq.setDouble('/dev8163/sines/1/phaseshift', 0)

#controls for oscillator amplitude for channel 0
daq.setInt('/dev8163/sines/0/enables/0', 1)
daq.setDouble('/dev8163/sines/0/amplitudes/0', 0.8)
daq.setInt('/dev8163/sines/1/enables/0', 1)
daq.setDouble('/dev8163/sines/1/amplitudes/0', 0.2)
#turn on channel 0
daq.setDouble('/dev8163/sigouts/0/range', 2)
daq.setInt('/dev8163/sigouts/0/on', 1)

#controls for oscillator amplitude for channel 1
daq.setInt('/dev8163/sines/0/enables/1', 1)
daq.setDouble('/dev8163/sines/0/amplitudes/1', 0.8)
daq.setInt('/dev8163/sines/1/enables/1', 1)
daq.setDouble('/dev8163/sines/1/amplitudes/1', 0.8)
#turn on channel 1
daq.setDouble('/dev8163/sigouts/1/range', 2)
daq.setInt('/dev8163/sigouts/1/on', 1)

#controls for oscillator amplitude for channel 2
daq.setInt('/dev8163/sines/2/enables/0', 1)
daq.setDouble('/dev8163/sines/2/amplitudes/0', 0.8)
daq.setInt('/dev8163/sines/3/enables/0', 0)
daq.setDouble('/dev8163/sines/3/amplitudes/1', 0.8)
#turn on channel 2
daq.setDouble('/dev8163/sigouts/2/range', 2)
daq.setInt('/dev8163/sigouts/2/on', 1)

#controls for oscillator amplitude for channel 3
daq.setInt('/dev8163/sines/2/enables/1', 1)
daq.setDouble('/dev8163/sines/2/amplitudes/0', 0.8)
daq.setInt('/dev8163/sines/3/enables/1', 0)
daq.setDouble('/dev8163/sines/3/amplitudes/1', 0.8)
#turn on channel 3
daq.setDouble('/dev8163/sigouts/3/range', 2)
daq.setInt('/dev8163/sigouts/3/on', 1)

## Starting module awgModule on 2020/01/23 18:09:28
#h = daq.awgModule()
#h.set('awgModule/device', 'dev8163')
#h.set('awgModule/index', 0)
#h.execute()
#h.set('awgModule/compiler/sourcestring', 'wave w_gauss = 0.5*gauss(8000,4000,1000);\
#wave w_drag = 0.5*drag(8000,4000,1000);\
#while(true){\
#  setTrigger(1);\
#  playWave(w_gauss,0.5*w_gauss,-0.5*w_gauss,-w_gauss);\
#  playWave(2,w_gauss);\
#  playWave(3,w_gauss);\
#  playWave(4,w_gauss);\
#  waitWave();\
#  setTrigger(0);\
#  wait(10000);\
#}')
#daq.setInt('/dev8163/awgs/0/enable', 1)
#daq.setInt('/dev8163/sigouts/0/on', 1)
#daq.setInt('/dev8163/sigouts/1/on', 1)
#daq.setInt('/dev8163/sigouts/2/on', 1)
#daq.setInt('/dev8163/sigouts/3/on', 1)