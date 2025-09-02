'''
8-25-21 AK modifying to normalize the amplitudes to the drive power. Undid that.

9-2-21 AK made it return the data

'''


import os
import time
import numpy as np
import matplotlib.pyplot as plt
import mplcursors
mplcursors.cursor(multiple=True)

import userfuncs
from utility.plotting_tools import simplescan_plot
from utility.measurement_helpers import configure_card, configure_hdawg, estimate_time, read_and_process, read_and_process_singleshot
from utility.scheduler import scheduler
from pdh_measurements.scheduler_pdh import scheduler_pdh
from utility.userfits_v2 import fit_model

import scipy.signal as signal

def get_default_settings():
    settings = {}
    
    settings['scanname'] = 'scanname'
    settings['meas_type'] = 'PulsedTrans'
    
    #Cavity parameters
    settings['CAV_Power']        = -60
    settings['CAV_freq']        = 7e9
    
    #Qubit parameters
    settings['start_freq']      = 4.15*1e9  
    settings['stop_freq']       = 4.25*1e9 
    settings['freq_points']     = 50

    settings['start_power']     = -30
    settings['stop_power']      = -20
    settings['power_points']    = 5

    #Card settings
    settings['segments']         = 1
    settings['reads']            = 1
    settings['averages']         = 5e3
    
    #Measurement settings
    settings['Quasi_CW']    = False
    settings['reverse']   = False
    settings['num_save'] = 1
    
    #background_subtraction (by taking reference trace with no qubit drive power)
    settings['subtract_background'] = False
    
    return settings


def pulsed_trans_GEF(instruments, settings):
    
        ##Instruments used
    Qbitgen  = instruments['Qbitgen']
    cavitygen = instruments['cavitygen']
    card      = instruments['card']
    hdawg     = instruments['AWG']
    LO        = instruments['LO']
    
    exp_globals  = settings['exp_globals']
    exp_settings = settings['exp_settings']

    ##Data saving and naming
    stamp    = userfuncs.timestamp()
    saveDir  = userfuncs.saveDir(settings)
    filename = exp_settings['scanname'] + '_' + stamp
    
    ##Cavity settings
    CAV_Attenuation = exp_globals['CAV_Attenuation']
    # cav_freq  = exp_settings['CAV_freq']
    cav_start_freq  = exp_settings['cav_start_freq']
    cav_stop_freq   = exp_settings['cav_stop_freq']
    cav_freq_points = exp_settings['cav_freq_points']
    cav_freqs  = np.round(np.linspace(cav_start_freq,cav_stop_freq,cav_freq_points),3)
    cav_start_power = exp_settings['cav_start_power'] + CAV_Attenuation
    cav_stop_power = exp_settings['cav_stop_power'] + CAV_Attenuation
    cav_power_points = exp_settings['cav_power_points']
    cav_powers = np.round(np.linspace(cav_start_power,cav_stop_power,cav_power_points),2)
    
    Qbit_Attenuation = exp_globals['Qbit_Attenuation']
    Qbit_E_freq = exp_settings['Qbit_E_freq']
    Qbit_E_power = exp_settings['Qbit_E_power'] + Qbit_Attenuation
    Qbit_F_freq = exp_settings['Qbit_F_freq']
    Qbit_F_power = exp_settings['Qbit_F_power'] + Qbit_Attenuation
    
    ## Generator settings
    cavitygen.Freq   = cav_start_freq 
    cavitygen.Power  = cav_powers[0]
    cavitygen.Output = 'On'

    Qbitgen.freq = Qbit_E_freq
    Qbitgen.power = Qbit_E_power
    Qbitgen.Output  = 'On'
    
    LO.power = 12
    LO.freq = cav_start_freq - exp_globals['IF']    
    LO.output = 'On'
    
    cavitygen.disable_IQ()
    cavitygen.enable_pulse()
    cavitygen.output = 'On'
    
    if exp_settings['Quasi_CW']:
        Qbitgen.disable_pulse()
        Qbitgen.disable_IQ()
    else:
        Qbitgen.enable_pulse()
        Qbitgen.enable_IQ()
    
    ## Card config
    configure_card(card, settings)
    
    ## Sequencer program
    progFile = open(r"C:\Users\Kollarlab\Desktop\Kollar-Lab\pulsed_measurements\HDAWG_sequencer_codes\hdawg_placeholder.cpp",'r')
    rawprog  = progFile.read()
    loadprog = rawprog
    progFile.close()
    
    m_pulse = exp_globals['measurement_pulse']
    q_pulse = exp_globals['qubit_pulse']
    
    start_time  = m_pulse['meas_pos']
    window_time = m_pulse['meas_window']
    
    awg_sched = scheduler(total_time=start_time+2*window_time, sample_rate=2.4e9)

    awg_sched.add_analog_channel(1, name='Qbit_I')
    awg_sched.add_analog_channel(2, name='Qbit_Q')
    
    awg_sched.add_digital_channel(1, name='Qbit_enable', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
    awg_sched.add_digital_channel(2, name='Cavity_enable', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
    
    
    Qbit_I       = awg_sched.analog_channels['Qbit_I']
    Qbit_marker  = awg_sched.digital_channels['Qbit_enable']
    cavity_marker = awg_sched.digital_channels['Cavity_enable']
    
    delay = q_pulse['delay']
    sigma = q_pulse['sigma']
    num_sigma = q_pulse['num_sigma']
    # init_buffer = m_pulse['init_buffer']

    position = start_time-delay-num_sigma*sigma-q_pulse['hold_time']
    Qbit_I.add_pulse('gaussian_square', position=position, 
                              amplitude=q_pulse['piAmp'], length = q_pulse['hold_time'], 
                              ramp_sigma=q_pulse['sigma'], num_sigma=q_pulse['num_sigma'])
    
    Qbit_marker.add_window(position-160e-9, position+2*160e-9+q_pulse['hold_time'])
    cavity_marker.add_window(start_time, start_time+window_time)
    
    awg_sched.plot_waveforms()
    
    [ch1, ch2, marker] = awg_sched.compile_schedule('HDAWG', ['Qbit_I', 'Qbit_Q'], ['Qbit_enable', 'Cavity_enable'])
    
    loadprog = loadprog.replace('_samples_', str(awg_sched.samples))
    hdawg.AWGs[0].load_program(loadprog)
    hdawg.AWGs[0].load_waveform('0', ch1, ch2, marker)
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
    
    first_it = True
    segments = int(card.segments)
    total_samples = int(card.samples)

    # Data saving arrays
    I_all_powers_g = np.zeros((len(cav_powers), len(cav_freqs), segments))
    Q_all_powers_g = np.zeros((len(cav_powers), len(cav_freqs), segments))

    I_all_powers_e = np.zeros((len(cav_powers), len(cav_freqs), segments))
    Q_all_powers_e = np.zeros((len(cav_powers), len(cav_freqs), segments))

    I_all_powers_f = np.zeros((len(cav_powers), len(cav_freqs), segments))
    Q_all_powers_f = np.zeros((len(cav_powers), len(cav_freqs), segments))

    mags_all_powers_g = np.zeros((len(cav_powers), len(cav_freqs), segments))
    mags_all_powers_e = np.zeros((len(cav_powers), len(cav_freqs), segments))
    mags_all_powers_f = np.zeros((len(cav_powers), len(cav_freqs), segments))

    phases_all_powers_g = np.zeros((len(cav_powers), len(cav_freqs), segments))
    phases_all_powers_e = np.zeros((len(cav_powers), len(cav_freqs), segments))
    phases_all_powers_f = np.zeros((len(cav_powers), len(cav_freqs), segments))

    Is_full_all_powers, Qs_full_all_powers = np.zeros((len(cav_powers), len(cav_freqs), total_samples)), np.zeros((len(cav_powers), len(cav_freqs), total_samples))
    Is_back_all_powers, Qs_back_all_powers = np.zeros((len(cav_powers), len(cav_freqs), total_samples)), np.zeros((len(cav_powers), len(cav_freqs), total_samples))

    for powerind in range(len(cav_powers)):
        cavitygen.Power = cav_powers[powerind]
        
        print('Current power:{}, max:{}'.format(cav_powers[powerind]-CAV_Attenuation, cav_powers[-1]-CAV_Attenuation))
    
        for find in range(0, len(cav_freqs)):
            
            if first_it:
                tstart = time.time()
            
            # Acquire  g state signal
            cav_freq = cav_freqs[find]
            cavitygen.freq = cav_freq
            cavitygen.phase = 0
            cavitygen.output = 'On'
            LO.freq = cav_freq - exp_globals['IF'] 
            LO.phase = 0
            LO.output = 'On'
            Qbitgen.output='Off'  
            Qbitgen.phase = 0

            I_window_g, Q_window_g, I_full_g, Q_full_g, xaxis = read_and_process_singleshot(card, settings, 
                                                                         plot=first_it, 
                                                                         IQstorage = True)
            # Acquire e state signal
            Qbitgen.freq = Qbit_E_freq
            Qbitgen.power = Qbit_E_power
            Qbitgen.output = 'On'

            I_window_e, Q_window_e, I_full_e, Q_full_e, xaxis = read_and_process_singleshot(card, settings, 
                                                                         plot=first_it, 
                                                                         IQstorage = True)
            # Acquire f state signal
            Qbitgen.freq = Qbit_F_freq
            Qbitgen.power = Qbit_F_power
            Qbitgen.output = 'On'

            I_window_f, Q_window_f, I_full_f, Q_full_f, xaxis = read_and_process_singleshot(card, settings, 
                                                                         plot=first_it, 
                                                                         IQstorage = True)
            
            if exp_settings['subtract_background']:
                #Acquire background trace
                Qbitgen.output='Off'
#                time.sleep(0.1)
                I_window_b, Q_window_b, I_full_b, Q_full_b, xaxis_b = read_and_process_singleshot(card, settings, 
                                                                 plot=first_it, 
                                                                 IQstorage = True)
            else:
                I_window_b, Q_window_b, I_full_b, Q_full_b = np.zeros(I_window_g.shape), np.zeros(Q_window_g.shape), np.zeros(I_full_g.shape), np.zeros(Q_full_g.shape)
            
            if first_it:
                first_it=False

            ##Useful handles for variables
            I_final_g = np.mean(I_window_g, axis=1, keepdims=True) - np.mean(I_window_b, axis=1, keepdims=True)
            Q_final_g = np.mean(Q_window_g, axis=1, keepdims=True) - np.mean(Q_window_b, axis=1, keepdims=True)

            I_final_e = np.mean(I_window_e, axis=1, keepdims=True) - np.mean(I_window_b, axis=1, keepdims=True)
            Q_final_e = np.mean(Q_window_e, axis=1, keepdims=True) - np.mean(Q_window_b, axis=1, keepdims=True)

            I_final_f = np.mean(I_window_f, axis=1, keepdims=True) - np.mean(I_window_b, axis=1, keepdims=True)
            Q_final_f = np.mean(Q_window_f, axis=1, keepdims=True) - np.mean(Q_window_b, axis=1, keepdims=True)

            I_final_g, Q_final_g = I_final_g.flatten(), Q_final_g.flatten()
            I_final_e, Q_final_e = I_final_e.flatten(), Q_final_e.flatten()
            I_final_f, Q_final_f = I_final_f.flatten(), Q_final_f.flatten()

            I_all_powers_g[powerind, find, :], Q_all_powers_g[powerind, find, :] = I_final_g, Q_final_g
            mags_all_powers_g[powerind, find, :] = np.sqrt(I_final_g**2 + Q_final_g**2)
            phases_all_powers_g[powerind, find, :] = np.arctan2(Q_final_g, I_final_g)*180/np.pi

            I_all_powers_e[powerind, find, :], Q_all_powers_e[powerind, find, :] = I_final_e, Q_final_e
            mags_all_powers_e[powerind, find, :] = np.sqrt(I_final_e**2 + Q_final_e**2)
            phases_all_powers_e[powerind, find, :] = np.arctan2(Q_final_e, I_final_e)*180/np.pi

            I_all_powers_f[powerind, find, :], Q_all_powers_f[powerind, find, :] = I_final_f, Q_final_f
            mags_all_powers_f[powerind, find, :] = np.sqrt(I_final_f**2 + Q_final_f**2)
            phases_all_powers_f[powerind, find, :] = np.arctan2(Q_final_f, I_final_f)*180/np.pi

            Is_full_all_powers[powerind, find, :], Qs_full_all_powers[powerind, find, :] = I_full_e[0], Q_full_e[0]
            Is_back_all_powers[powerind, find, :], Qs_back_all_powers[powerind, find, :] = I_full_b[0], Q_full_b[0]

    ## Plot IQ data points

    ##########################
    if cav_freq_points > 1:
        fig = plt.figure(11131, figsize=(13,8))
        plt.clf()
        ax1 = plt.subplot(2,2,1)
        xaxis = cav_freqs/1e9
        yaxis_g = mags_all_powers_g[0,:,0]
        params_cav_g = fit_model(xaxis*1e9, yaxis_g, 'lorenz', plot=False, ax=ax1)
        linewidth_g = round(2*params_cav_g['sigma']/1e6,2)
        res_freq_g = round(params_cav_g['center']/1e9,5)
        plt.plot(xaxis, yaxis_g, label = 'G state')
        plt.xlabel('Freq (GHz)')
        plt.title(f'Res Freq:{res_freq_g}GHz | Linewidth:{linewidth_g}MHz')
        plt.legend()
       
        ax2 = plt.subplot(2,2,2)
        xaxis = cav_freqs/1e9
        yaxis_e = mags_all_powers_e[0,:,0]
        params_cav_e = fit_model(xaxis*1e9, yaxis_e, 'lorenz', plot=False, ax=ax2)
        linewidth_e = round(2*params_cav_e['sigma']/1e6,2)
        res_freq_e = round(params_cav_e['center']/1e9,5)
        plt.plot(xaxis, yaxis_e, label = 'E state')
        plt.xlabel('Freq (GHz)')
        plt.title(f'Res Freq:{res_freq_e}GHz | Linewidth:{linewidth_e}MHz')
        plt.legend()

        ax3 = plt.subplot(2,2,3)
        xaxis = cav_freqs/1e9
        yaxis_f = mags_all_powers_f[0,:,0]
        params_cav_f = fit_model(xaxis*1e9, yaxis_f, 'lorenz', plot=False, ax=ax3)
        linewidth_f = round(2*params_cav_f['sigma']/1e6,2)
        res_freq_f = round(params_cav_f['center']/1e9,5)
        plt.plot(xaxis, yaxis_f, label = 'F state')
        plt.xlabel('Freq (GHz)')
        plt.title(f'Res Freq:{res_freq_f}GHz | Linewidth:{linewidth_f}MHz')
        plt.legend()

        ax4 = plt.subplot(2,2,4)
        xaxis = cav_freqs/1e9
        plt.plot(xaxis, yaxis_g, label = 'G state')
        plt.plot(xaxis, yaxis_e, label = 'E state')
        plt.plot(xaxis, yaxis_f, label = 'F state')
        plt.xlabel('Freq (GHz)')
        plt.title(f'$\chi$:{round((res_freq_g - res_freq_e)*1e3/2,3)}MHz | $\kappa$g:{linewidth_g}MHz')
        plt.legend()

        plt.tight_layout()
        mplcursors.cursor(multiple=True)
        plt.savefig(os.path.join(saveDir, filename+'AllSweeps.png'), dpi=150)
        plt.show()
    else:
    #################3########
        fig = plt.figure(1110, figsize=(16,10))
        plt.clf()
        ax = plt.subplot(2, 2, 1)
        ax.scatter(I_all_powers_g[0, 0, :], Q_all_powers_g[0, 0, :], s=4, label = 'G')
        plt.xlabel('I')
        plt.ylabel('Q')
        plt.title('G state IQ')
        plt.legend()
        plt.axhline(y = 0, color ="black", linestyle ="-")
        plt.axvline(x = 0, color ="black", linestyle ="-") 
        y_max = np.max(np.abs(ax.get_ylim()))
        x_max = np.max(np.abs(ax.get_xlim()))
        list = np.array([x_max, y_max])
        max = np.max(list)
        ax.set_ylim(ymin=-max, ymax=max)
        ax.set_xlim(xmin=-max, xmax=max)
        ax.set_aspect('equal') 

        ax = plt.subplot(2, 2, 2)
        ax.scatter(I_all_powers_e[0, 0, :], Q_all_powers_e[0, 0, :], s=4, label = 'E', color = 'orange')
        plt.xlabel('I')
        plt.ylabel('Q')
        plt.title('E state IQ')
        plt.legend()
        plt.axhline(y = 0, color ="black", linestyle ="-")
        plt.axvline(x = 0, color ="black", linestyle ="-") 
        ax.set_ylim(ymin=-max, ymax=max)
        ax.set_xlim(xmin=-max, xmax=max)
        ax.set_aspect('equal') 

        ax = plt.subplot(2, 2, 3)
        ax.scatter(I_all_powers_f[0, 0, :], Q_all_powers_f[0, 0, :], s=4, label = 'F', color='green')
        plt.xlabel('I')
        plt.ylabel('Q')
        plt.title('F state IQ')
        plt.legend()
        plt.axhline(y = 0, color ="black", linestyle ="-")
        plt.axvline(x = 0, color ="black", linestyle ="-") 
        ax.set_ylim(ymin=-max, ymax=max)
        ax.set_xlim(xmin=-max, xmax=max)
        ax.set_aspect('equal') 

        ax = plt.subplot(2, 2, 4)
        ax.scatter(I_all_powers_g[0, 0, :], Q_all_powers_g[0, 0, :], s=4, label = 'G')
        ax.scatter(I_all_powers_e[0, 0, :], Q_all_powers_e[0, 0, :], s=4, label = 'E', color = 'orange')
        ax.scatter(I_all_powers_f[0, 0, :], Q_all_powers_f[0, 0, :], s=4, label = 'F', color = 'green')
        plt.xlabel('I')
        plt.ylabel('Q')
        plt.title('GEF state IQ')
        plt.legend()
        plt.axhline(y = 0, color ="black", linestyle ="-")
        plt.axvline(x = 0, color ="black", linestyle ="-") 
        ax.set_ylim(ymin=-max, ymax=max)
        ax.set_xlim(xmin=-max, xmax=max)
        ax.set_aspect('equal') 
    
        fig.suptitle(f'{segments} points | Readout:{round(cav_start_freq/1e9,3)}GHz, {cav_start_power - CAV_Attenuation}dBm | Fq_e:{round(Qbit_E_freq/1e9,3)}GHz | Fq_f:{round(Qbit_F_freq/1e9,3)}GHz')
        fig.tight_layout()
        plt.savefig(os.path.join(saveDir, filename+'_IQPlot.png'), dpi = 150)
        plt.show()

        # plot above with imshow, cmap = hot
        bins = 100
        fig = plt.figure(1112, figsize=(16,10))
        plt.clf()
        ax = plt.subplot(2, 2, 1)
        I_array_g, Q_array_g = I_all_powers_g[0, 0, :], Q_all_powers_g[0, 0, :]
        heatmap, xedges, yedges = np.histogram2d(I_array_g, Q_array_g, bins=bins)
        max_count = np.max(heatmap)
        heatmap_scaled = heatmap/max_count
        plt.imshow(heatmap_scaled.T, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
                origin='lower', cmap='hot', vmin=0.05, vmax=0.5)
        plt.colorbar(label=fr'Counts $\times$ {int(max_count)}', fraction=0.04, pad=0.04)        
        plt.xlabel('I')
        plt.ylabel('Q')
        plt.title('G state IQ')
        plt.axhline(y = 0, color ="black", linestyle ="-")
        plt.axvline(x = 0, color ="black", linestyle ="-") 
        y_max = np.max(np.abs(ax.get_ylim()))
        x_max = np.max(np.abs(ax.get_xlim()))
        list = np.array([x_max, y_max])
        max = np.max(list)
        ax.set_ylim(ymin=-max, ymax=max)
        ax.set_xlim(xmin=-max, xmax=max)
        plt.gca().set_facecolor('black')

        ax = plt.subplot(2, 2, 2)
        I_array_e, Q_array_e = I_all_powers_e[0, 0, :], Q_all_powers_e[0, 0, :]
        heatmap, xedges, yedges = np.histogram2d(I_array_e, Q_array_e, bins=bins)
        max_count = np.max(heatmap)
        heatmap_scaled = heatmap/max_count
        plt.imshow(heatmap_scaled.T, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
                origin='lower', cmap='hot', vmin=0.05, vmax=0.5)
        plt.colorbar(label=fr'Counts $\times$ {int(max_count)}', fraction=0.04, pad=0.04)        
        plt.xlabel('I')
        plt.ylabel('Q')
        plt.title('E state IQ')
        plt.axhline(y = 0, color ="black", linestyle ="-")
        plt.axvline(x = 0, color ="black", linestyle ="-") 
        ax.set_ylim(ymin=-max, ymax=max)
        ax.set_xlim(xmin=-max, xmax=max)
        plt.gca().set_facecolor('black')

        ax = plt.subplot(2, 2, 3)
        I_array_f, Q_array_f = I_all_powers_f[0, 0, :], Q_all_powers_f[0, 0, :]
        heatmap, xedges, yedges = np.histogram2d(I_array_f, Q_array_f, bins=bins)
        max_count = np.max(heatmap)
        heatmap_scaled = heatmap/max_count
        plt.imshow(heatmap_scaled.T, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
                origin='lower', cmap='hot', vmin=0.05, vmax=0.5)
        plt.colorbar(label=fr'Counts $\times$ {int(max_count)}', fraction=0.04, pad=0.04)        
        plt.xlabel('I')
        plt.ylabel('Q')
        plt.title('F state IQ')
        plt.axhline(y = 0, color ="black", linestyle ="-")
        plt.axvline(x = 0, color ="black", linestyle ="-") 
        ax.set_ylim(ymin=-max, ymax=max)
        ax.set_xlim(xmin=-max, xmax=max)
        plt.gca().set_facecolor('black')

        # Fourth subplot (All states IQ combined)
        ax = plt.subplot(2, 2, 4)

        # Combine all the arrays: G, E, F states
        I_combined = np.concatenate((I_array_g, I_array_e, I_array_f))
        Q_combined = np.concatenate((Q_array_g, Q_array_e, Q_array_f))
        heatmap, xedges, yedges = np.histogram2d(I_combined, Q_combined, bins=bins)
        max_count = np.max(heatmap)
        heatmap_scaled = heatmap / max_count
        plt.imshow(heatmap_scaled.T, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
                origin='lower', cmap='hot', vmin=0.05, vmax=0.5)
        plt.colorbar(label=fr'Counts $\times$ {int(max_count)}', fraction=0.04, pad=0.04)
        plt.xlabel('I')
        plt.ylabel('Q')
        plt.title('All States IQ Combined')
        plt.axhline(y=0, color="black", linestyle="-")
        plt.axvline(x=0, color="black", linestyle="-")
        ax.set_ylim(ymin=-max, ymax=max)
        ax.set_xlim(xmin=-max, xmax=max)
        plt.gca().set_facecolor('black')

        fig.suptitle(f'{segments} points | Readout:{round(cav_start_freq/1e9,3)}GHz, {cav_start_power - CAV_Attenuation}dBm | Fq_e:{round(Qbit_E_freq/1e9,3)}GHz | Fq_f:{round(Qbit_F_freq/1e9,3)}GHz')
        fig.tight_layout()
        plt.savefig(os.path.join(saveDir, filename+'_IQHotPlot.png'), dpi = 150)
        plt.show()
                

    # data to save
    full_data = {}
    full_data['I_all_powers_g'] = I_all_powers_g
    full_data['Q_all_powers_g'] = Q_all_powers_g
    full_data['I_all_powers_e'] = I_all_powers_e
    full_data['Q_all_powers_e'] = Q_all_powers_e
    full_data['I_all_powers_f'] = I_all_powers_f
    full_data['Q_all_powers_f'] = Q_all_powers_f
    full_data['Is_full_all_powers'] = Is_full_all_powers
    full_data['Qs_full_all_powers'] = Qs_full_all_powers
    full_data['Is_background'] = Is_back_all_powers
    full_data['Qs_background'] = Qs_back_all_powers

    userfuncs.SaveFull(saveDir, filename, ['cav_powers','cav_freqs', 'xaxis',
                                                       'full_data'], 
                                                     locals(), 
                                                     expsettings=settings, 
                                                     instruments=instruments, saveHWsettings=True)    


    t2 = time.time()
    
    print('elapsed time = ' + str(t2-tstart))
       
    cavitygen.Output = 'Off'
    Qbitgen.Output = 'Off'
    LO.output = 'Off'
    
    return full_data
