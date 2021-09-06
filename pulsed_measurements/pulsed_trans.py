'''
8-25-21 AK modifying to normalize the amplitudes to the drive power.

8-25-21 AK modifying to return the raw data.

'''


import os
import time
import numpy as np 
import matplotlib.pyplot as plt

import userfuncs
from utility.plotting_tools import simplescan_plot
from utility.measurement_helpers import configure_card, configure_hdawg, estimate_time, read_and_process

def get_default_settings():
    settings = {}
    
    settings['scanname'] = 'initial_power_scan_q4'
    settings['meas_type'] = 'PulsedTrans'
    
    #Sweep parameters
    settings['start_freq']   = 4.15*1e9  
    settings['stop_freq']    = 4.25*1e9 
    settings['freq_points']  = 50

    settings['start_power']  = -20
    settings['stop_power']   = 10
    settings['power_points'] = 31

    #Card settings
    settings['segments'] = 1
    settings['reads']    = 1
    settings['averages'] = 5e3
    
    return settings

def pulsed_trans(instruments, settings):
    ##Instruments used
    cavitygen = instruments['cavitygen']
    card      = instruments['card']
    hdawg     = instruments['AWG']
    LO        = instruments['LO']

    exp_globals  = settings['exp_globals']
    exp_settings = settings['exp_settings'] 

    stamp    = userfuncs.timestamp()
    saveDir  = userfuncs.saveDir(settings)
    filename = exp_settings['scanname'] + '_' + stamp

    ##Sweep settings
    start_freq  = exp_settings['start_freq']
    stop_freq   = exp_settings['stop_freq']
    freq_points = exp_settings['freq_points']
    freqs  = np.round(np.linspace(start_freq,stop_freq,freq_points),-3)
    
    CAV_Attenuation = exp_globals['CAV_Attenuation']
    start_power  = exp_settings['start_power'] + CAV_Attenuation
    stop_power   = exp_settings['stop_power'] + CAV_Attenuation
    power_points = exp_settings['power_points']
    powers = np.round(np.linspace(start_power,stop_power,power_points),2)
    
    ## Generator settings
    cavitygen.Freq   = freqs[0]
    cavitygen.Power  = powers[0]
    cavitygen.IQ.Mod = 'On'

    cavitygen.Output = 'On'
    
    LO.Power  = 12
    LO.Freq   = freqs[0] - exp_globals['IF']
    LO.Output = 'On'
    
    ##Card settings
    configure_card(card, settings)

    ##HDAWG settings
    configure_hdawg(hdawg, settings)
    
    progFile = open(r"C:\Users\Kollarlab\Desktop\Kollar-Lab\pulsed_measurements\HDAWG_sequencer_codes\pulsedtrans.cpp",'r')
    rawprog  = progFile.read()
    loadprog = rawprog
    progFile.close()

    m_pulse  = exp_globals['measurement_pulse'] 
    loadprog = loadprog.replace('_max_time_', str(m_pulse['meas_pos']))
    loadprog = loadprog.replace('_meas_window_', str(m_pulse['meas_window']))
    
    hdawg.AWGs[0].load_program(loadprog)
    hdawg.AWGs[0].run_loop()
    time.sleep(0.1)

    ## Start main measurement loop 
    powerdat = np.zeros((len(powers), len(freqs)))
    phasedat = np.zeros((len(powers), len(freqs)))
    
    tstart = time.time()
    first_it = True

    drive_powers_lin = 10**(powers/10)
    drive_amps_lin = np.sqrt(drive_powers_lin)

    for powerind in range(len(powers)):
        cavitygen.Power = powers[powerind]
        time.sleep(0.2)

        total_samples = card.samples 
        amps   = np.zeros((len(freqs), total_samples))
        phases = np.zeros((len(freqs), total_samples))
        
        print('Current power:{}, max:{}'.format(powers[powerind] - CAV_Attenuation, powers[-1] - CAV_Attenuation))
    
        for find in range(0, len(freqs)):
            freq = freqs[find]

            cavitygen.freq = freq

            LO.freq   = freq-exp_globals['IF']
            LO.output = 'On'

            time.sleep(0.2)

            amp, phase, amp_full, phase_full, xaxis = read_and_process(card, settings, plot=first_it)
            
            if first_it:
                tstop = time.time()
                estimate_time(tstart, tstop, len(powers)*len(freqs))
                first_it = False
            
            amps[find,:]   = amp_full
            phases[find,:] = phase_full
            powerdat[powerind, find] = np.mean(amp)
            phasedat[powerind, find] = np.mean(phase)   ####!!!! this does not work with heterodyne. Because the phase is wrapping.
            
        full_data = {}
        full_data['xaxis']  = freqs/1e9
        full_data['mags']   = powerdat[0:powerind+1]
        full_data['phases'] = phasedat[0:powerind+1]
        
        #rescale to fractional amplitude for the plots.
        plot_data = {}
        plot_data['xaxis']  = freqs/1e9
        plot_data['mags']   = np.transpose(   np.transpose(powerdat[0:powerind+1])/drive_amps_lin[0:powerind+1])
        plot_data['phases'] = phasedat[0:powerind+1]

        single_data = {}
        single_data['xaxis'] = freqs/1e9
        single_data['mag']   = powerdat[powerind]
        single_data['phase'] = phasedat[powerind]

        yaxis  = powers[0:powerind+1] - CAV_Attenuation
        labels = ['Freq (GHz)', 'Power (dBm)']
#        simplescan_plot(full_data, single_data, yaxis, filename, labels, identifier='', fig_num=1) #unnormalized plot
        simplescan_plot(plot_data, single_data, yaxis, filename, labels, identifier='', fig_num=1) #normalized to drive level
        plt.savefig(os.path.join(saveDir, filename+'_fullColorPlot.png'), dpi = 150)

        full_time = {}
        full_time['xaxis']  = xaxis*1e6
        full_time['mags']   = amps
        full_time['phases'] = phases

        single_time = {}
        single_time['xaxis'] = xaxis*1e6
        single_time['mag']   = amp_full
        single_time['phase'] = phase_full

        time_labels = ['Time (us)', 'Freq (GHz)']
        identifier = 'Power: {}dBm'.format(powers[powerind]-CAV_Attenuation)
        simplescan_plot(full_time, single_time, freqs/1e9, 'Raw_time_traces\n'+filename, time_labels, identifier, fig_num=2)
        plt.savefig(os.path.join(saveDir, filename+'_RawTimeTraces.png'), dpi = 150)

        userfuncs.SaveFull(saveDir, filename, ['powers','freqs', 'powerdat', 'phasedat','xaxis','full_data', 'single_data', 'full_time', 'single_time'],
                             locals(), expsettings=settings, instruments=instruments)
    t2 = time.time()
    
    print('elapsed time = ' + str(t2-tstart))
    
    return full_data