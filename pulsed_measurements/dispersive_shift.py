
import os
import time
import numpy as np 
import matplotlib.pyplot as plt

import userfuncs
from utility.plotting_tools import simplescan_plot
from utility.measurement_helpers import configure_card, configure_hdawg, estimate_time, read_and_process

import scipy.signal as signal

def get_default_settings():
    settings = {}
    
    settings['scanname']  = 'test'
    settings['meas_type'] = 'Dispersive_shift'
    
    settings['CAV_Power'] = -30
    settings['Q_Power']   = -30
    
    #Sweep parameters
    settings['start_freq']   = 4.15*1e9  
    settings['stop_freq']    = 4.25*1e9 
    settings['freq_points']  = 50

    settings['start_freq_q']  = 4e9
    settings['stop_freq_q']   = 5e9
    settings['q_freq_points'] = 31

    #Card settings
    settings['segments'] = 1
    settings['reads']    = 1
    settings['averages'] = 5e3
    
    return settings

def dispersive_shift(instruments, settings):
    ##Instruments used
    cavitygen = instruments['cavitygen']
    qubitgen  = instruments['qubitgen']
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
    
    start_freq_q  = exp_settings['start_freq_q']
    stop_freq_q   = exp_settings['stop_freq_q']
    q_freq_points = exp_settings['q_freq_points']
    qfreqs  = np.round(np.linspace(start_freq_q,stop_freq_q,q_freq_points),-3)
    
    ## Generator settings
    CAV_Attenuation  = exp_globals['CAV_Attenuation']
    Qbit_Attenuation = exp_globals['Qbit_Attenuation']
    cavitygen.Freq   = freqs[0]
    cavitygen.Power  = exp_settings['CAV_Power'] + CAV_Attenuation
    qubitgen.Power   = exp_settings['Q_Power'] + Qbit_Attenuation
    
    if cavitygen.instrument_type == 'SGS':
        cavitygen.IQ.Mod = 'On'
    else:
        cavitygen.Mod = 'On'
    cavitygen.Output = 'On'

    LO.power  = 12
    LO.freq = '{} GHz'.format((cavitygen.Freq-exp_globals['IF'])/1e9)
    LO.output = 'On'
    
    ##Card settings
    configure_card(card, settings)

    ##HDAWG settings
    configure_hdawg(hdawg, settings)
    
    progFile = open(r"C:\Users\Kollarlab\Desktop\Kollar-Lab\pulsed_measurements\HDAWG_sequencer_codes\T1.cpp",'r')
    rawprog  = progFile.read()
    loadprog = rawprog
    progFile.close()
    
    q_pulse = exp_globals['qubit_pulse']
    m_pulse = exp_globals['measurement_pulse']
    loadprog = loadprog.replace('_max_time_', str(m_pulse['meas_pos']))
    loadprog = loadprog.replace('_meas_window_', str(m_pulse['meas_window']))
    loadprog = loadprog.replace('_tau_',str(q_pulse['delay']))
    loadprog = loadprog.replace('_qsigma_',str(q_pulse['sigma']))
    loadprog = loadprog.replace('_num_sigma_',str(q_pulse['num_sigma']))
    loadprog = loadprog.replace('_piAmp_',str(q_pulse['piAmp']))
    
    hdawg.AWGs[0].load_program(loadprog)
    hdawg.AWGs[0].run_loop()
    time.sleep(0.1)
    
    ##create the digital down conversion filter if needed.
    if exp_globals['IF'] != 0:
        #create Chebychev type II digital filter
        filter_N = exp_globals['ddc_config']['order']
        filter_rs = exp_globals['ddc_config']['stop_atten']
        filter_cutoff = np.abs(exp_globals['ddc_config']['cutoff'])
        LPF = signal.cheby2(filter_N, filter_rs, filter_cutoff, btype='low', analog=False, output='sos', fs=card.sampleRate)
        
        xaxis = np.arange(0, card.samples, 1) * 1/card.sampleRate
        digLO_sin = np.sin(2*np.pi*exp_globals['IF']*xaxis)
        digLO_cos = np.cos(2*np.pi*exp_globals['IF']*xaxis)
        
        #store in settings so that the processing functions can get to them
        settings['digLO_sin'] = digLO_sin 
        settings['digLO_cos'] = digLO_cos
        settings['LPF'] = LPF

    ## Start main measurement loop 
    freqdat = np.zeros((len(qfreqs), len(freqs)))
    phasedat = np.zeros((len(qfreqs), len(freqs)))
    
    tstart = time.time()
    first_it = True

    for qind in range(len(qfreqs)):

        qubitgen.freq = qfreqs[qind]
        
        time.sleep(0.2)

        total_samples = card.samples 
        
        Is  = np.zeros((len(freqs), total_samples)) #the very rawest data is thrown away for heterodyne! Warning
        Qs  = np.zeros((len(freqs), total_samples))
        
        print('Current freq:{}, max:{}'.format(qfreqs[qind], qfreqs[-1]))
    
        for find in range(0, len(freqs)):
            freq = freqs[find]

            cavitygen.freq = freq

            #LO.freq   = freq-exp_globals['IF']
            LO.freq = '{} GHz'.format((cavitygen.Freq-exp_globals['IF'])/1e9)
            LO.output = 'On'

            time.sleep(0.2)

#            amp, phase, amp_full, phase_full, xaxis = read_and_process(card, settings, plot=first_it)
            I_window, Q_window, I_full, Q_full, xaxis = read_and_process(card, settings, 
                                                                         plot=first_it, 
                                                                         IQstorage = True)
            
            I_final = np.mean(I_window) #compute <I> in the data window
            Q_final = np.mean(Q_window) #compute <Q> in the data window
            
            if first_it:
                tstop = time.time()
                estimate_time(tstart, tstop, len(qfreqs)*len(freqs))
                first_it = False
            
            Is[find,:]   = I_full
            Qs[find,:] = Q_full
            freqdat[qind, find] = np.sqrt(I_final**2 + Q_final**2)
            phasedat[qind, find] = np.arctan2(Q_final, I_final)*180/np.pi
            
        full_data = {}
        full_data['xaxis']  = freqs/1e9
        full_data['mags']   = freqdat[0:qind+1]
        full_data['phases'] = phasedat[0:qind+1]

        single_data = {}
        single_data['xaxis'] = freqs/1e9
        single_data['mag']   = freqdat[qind]
        single_data['phase'] = phasedat[qind]

        yaxis  = qfreqs[0:qind+1]
        labels = ['Freq (GHz)', 'Drive (GHz)']
        simplescan_plot(full_data, single_data, yaxis, filename, labels, identifier='', fig_num=1) #unnormalized plot
        plt.savefig(os.path.join(saveDir, filename+'_fullColorPlot.png'), dpi = 150)

        full_time = {}
        full_time['xaxis']  = xaxis*1e6
        full_time['Is']   = Is
        full_time['Qs'] = Qs

        single_time = {}
        single_time['xaxis'] = xaxis*1e6
        single_time['I']   = I_full
        single_time['Q'] = Q_full

        time_labels = ['Time (us)', 'Freq (GHz)']
        identifier = 'Power: {}dBm'.format(qfreqs[qind])
        simplescan_plot(full_time, single_time, freqs/1e9, 
                        'Raw_time_traces\n'+filename, 
                        time_labels, 
                        identifier, 
                        fig_num=2,
                        IQdata = True)
        plt.savefig(os.path.join(saveDir, filename+'_RawTimeTraces.png'), dpi = 150)

        userfuncs.SaveFull(saveDir, filename, ['qfreqs','freqs', 'freqdat', 'phasedat','xaxis','full_data', 'single_data', 'full_time', 'single_time'],
                             locals(), expsettings=settings, instruments=instruments)
    t2 = time.time()
    
    print('elapsed time = ' + str(t2-tstart))
    
    return full_data