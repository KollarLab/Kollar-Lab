# -*- coding: utf-8 -*-
"""
Created on Thu Jul 17 15:39:22 2025

@author: Kollarlab
"""

import numpy as np
import matplotlib.pyplot as plt
from datetime import timedelta, datetime

import scipy.signal as signal


####import the normal functions that we are using without change
from utility.measurement_helpers import plot_data_extraction, extract_data, remove_IQ_ellipse, extract_data_heterodyne, plot_data_extraction_TripleDDC, extract_data_heterodyne_TripleDDC, extract_finalIQ


def read_only_spec(card, settings, plot, IQstorage = True):
    '''basic data grab function '''
    
    card.ArmAndWait()
    I, Q = card.ReadAllData()
    Ip = np.mean(I, 0)
    Qp = np.mean(Q, 0)
    
    return I, Q, Ip, Qp

def read_only_TripleDDC(card, settings, plot, IQstorage = True):
    '''basic data grab function '''
    
    card.ArmAndWait()
    I, Q = card.ReadAllData()
    Ip = I #np.mean(I, 0)
    Qp = Q #np.mean(Q, 0)
    
    return I, Q, Ip, Qp


def read_only_both(card, settings, plot, IQstorage = True):
    '''basic data grab function '''
    
    card.ArmAndWait()
    I, Q = card.ReadAllData()
    
    Ip_triple = I #np.mean(I, 0)
    Qp_triple = Q #np.mean(Q, 0)
    
    Ip_spec = np.mean(I, 0)
    Qp_spec = np.mean(Q, 0)
    
    return I, Q, Ip_triple, Qp_triple, Ip_spec, Qp_spec



def and_process(I, Q, Ip, Qp, card, settings, plot, IQstorage = True):
    '''The original version of this function wanted to convert to
    amplitude(t) and phase(t). This has been found to be problematic.
    To get the old version of operation, set
    IQstorage = False
    IQstorage = True will instead return I(t) and Q(t) for
    later processing.
    
    
    AK 7-17-25 - pulled the read out of here so we can do two processings on the same data.
    
    '''
    # card.ArmAndWait()
    # I, Q = card.ReadAllData()
    # Ip = np.mean(I, 0)
    # Qp = np.mean(Q, 0)
    
    mixer_config = settings['exp_globals']['mixer_config']
    Ipp, Qpp = remove_IQ_ellipse(Ip, Qp, mixer_config)
#    Ipp, Qpp = remove_slow_drift(Ipp, Qpp, time_el)
    xaxis = np.linspace(0, card.samples, card.samples, endpoint=False)/card.sampleRate
#    Idata, Itime, Ifull, time_full = extract_data(Ipp, xaxis, settings)
#    Qdata, Qtime, Qfull, time_full = extract_data(Qpp, xaxis, settings)
    
    if not IQstorage:
        Idata, Itime, Ifull, time_full = extract_data(Ipp, xaxis, settings)
        Qdata, Qtime, Qfull, time_full = extract_data(Qpp, xaxis, settings)
        amp_full = np.sqrt(Ifull**2+Qfull**2)
        
        if settings['exp_globals']['IF'] == 0:
            phase_full = np.arctan2(Qfull, Ifull)*180/np.pi
            
        else:
            raw_angle = np.arctan2(Qfull, Ifull)*180/np.pi
            phase_full = np.mod(raw_angle + 360*settings['exp_globals']['IF']*xaxis, 360)
    #        plt.figure(101)
    #        plt.clf()
    #        ax = plt.subplot(1,2,1)
    #        plt.plot(xaxis, raw_angle)
    #        plt.title('raw phase angle')
    #        ax = plt.subplot(1,2,2)
    #        plt.plot(xaxis, phase_full)
    #        plt.show()
    #        plt.title('(hopefully) corrected phase angle')
        amp = np.sqrt(Idata**2+Qdata**2)
        
        if settings['exp_globals']['IF'] == 0:
            phase = np.arctan2(Qdata, Idata)*180/np.pi
            
        else:
            raw_angle = np.arctan2(Qdata, Idata)*180/np.pi
            phase = np.mod(raw_angle + 360*settings['exp_globals']['IF']*Itime, 360)
            
        if plot:
            plot_data_extraction(amp, Itime, amp_full, time_full, Ifull, Qfull)
            
        return amp, phase, amp_full, phase_full, xaxis
    
    else:
        #do not convert to amp(t) and phase(t)
#        if plot:
#            amp = np.sqrt(Idata**2+Qdata**2) #this is just a guide to theeye for locating the pulse
#            amp_full = np.sqrt(Ifull**2+Qfull**2)
#            plot_data_extraction(amp, Itime, amp_full, time_full, Ifull, Qfull)

        if settings['exp_globals']['IF'] == 0:
#            Idata, Itime, Ifull, time_full = extract_data(Ipp, xaxis, settings)
#            Qdata, Qtime, Qfull, time_full = extract_data(Qpp, xaxis, settings)
#            if plot:
#                amp = np.sqrt(Idata**2+Qdata**2) #this is just a guide to theeye for locating the pulse
#                amp_full = np.sqrt(Ifull**2+Qfull**2)
#                plot_data_extraction(amp, Itime, amp_full, time_full, Ifull, Qfull)
            #matching newer names
            I_window, Itime, I_full, time_full = extract_data(Ipp, xaxis, settings)
            Q_window, Qtime, Q_full, time_full = extract_data(Qpp, xaxis, settings)
            
            if plot:
                amp = np.sqrt(I_window**2+Q_window**2) #this is just a guide to theeye for locating the pulse
                amp_full = np.sqrt(I_full**2+Q_full**2)
                plot_data_extraction(amp, Itime, amp_full, time_full, I_full, Q_full)
                
        else:
            #Quick fix to the problem of full cancellation in the DDC data. Turns out the channels got swapped
            #during a reshuffle and combining the channels incorrectly leads to near perfect cancellation (offset by
            #the mixer impairements). For now I just swapped the definition of I and Q to match the combination phase
            I_cos, I_sin, Itime, I_cos_full, I_sin_full, time_full = extract_data_heterodyne(Ipp, xaxis, settings)
            Q_cos, Q_sin, Qtime, Q_cos_full, Q_sin_full, time_full = extract_data_heterodyne(Qpp, xaxis, settings)
            
            #rotate the mixer Q signal back into I so they can be averaged properly
            theta = -np.pi/2
            Qprime_sin = Q_sin*np.cos(theta) + -Q_cos *np.sin(theta)
            Qprime_cos = Q_sin*np.sin(theta) + Q_cos*np.cos(theta)
            
            Qprime_sin_full = Q_sin_full*np.cos(theta) + -Q_cos_full *np.sin(theta)
            Qprime_cos_full = Q_sin_full*np.sin(theta) + Q_cos_full*np.cos(theta)
            
            Q_window  = (I_sin +  Qprime_sin)/2
            I_window = (I_cos + Qprime_cos)/2
            
            Q_full  = (I_sin_full +  Qprime_sin_full)/2
            I_full = (I_cos_full + Qprime_cos_full)/2
            
            if plot:
                fig0 = plt.figure(97)
                plt.clf()
                ax = plt.subplot(1,2,1)
                plt.plot(xaxis*1e6, Ip, label = 'raw V1')
                plt.plot(xaxis*1e6, Qp, label = 'raw V2')
                plt.xlabel('Time (us)')
                plt.ylabel('Voltage')
                plt.title('Card Voltages')
                ax.legend(loc = 'upper right')
                
                ax = plt.subplot(1,2,2)
                plt.plot(xaxis*1e6, np.sqrt(Ip**2 + Qp**2), label = 'crude amplitude')
                plt.plot(xaxis*1e6, np.sqrt(Ipp**2 + Qpp**2), label = 'mixer corrected')
                plt.xlabel('Time (us)')
                plt.ylabel('Voltage')
                plt.title('Rough Pulse Amplitude')
                ax.legend(loc = 'upper right')
                plt.show()
                fig0.canvas.draw()
                fig0.canvas.flush_events()
                
                amp = np.sqrt(Q_window**2+I_window**2) #this is just the I signal
                amp_full = np.sqrt(I_full**2+Q_full**2)
                plot_data_extraction(amp, Itime, amp_full, time_full, I_full, Q_full)
        return I_window, Q_window, I_full, Q_full, xaxis
    
    
    
    
# def and_process_TripleDDC(I, Q, Ip, Qp, card, settings, plot, IQstorage = True):

    card.ArmAndWait()
    I, Q = card.ReadAllData()
    Ip = I #np.mean(I, 0)
    Qp = Q #np.mean(Q, 0)
    
    
    mixer_config = settings['exp_globals']['mixer_config']
    Ipp, Qpp = remove_IQ_ellipse(Ip, Qp, mixer_config)
    xaxis = np.linspace(0, card.samples, card.samples, endpoint=False)/card.sampleRate
    
    I_filtered_carr, I_filtered_usb, I_filtered_lsb = extract_data_heterodyne_TripleDDC(Ipp, xaxis, settings)
    Q_filtered_carr, Q_filtered_usb, Q_filtered_lsb = extract_data_heterodyne_TripleDDC(Qpp, xaxis, settings)
    
    # Extract final IQ for carrier
    I_window_carr, Q_window_carr, I_time_carr, I_full_carr, Q_full_carr, I_time_full_carr = extract_finalIQ(I_filtered_carr, Q_filtered_carr)
        
    # Extract final IQ for upper sideband
    I_window_usb, Q_window_usb, I_time_usb, I_full_usb, Q_full_usb, I_time_full_usb = extract_finalIQ(I_filtered_usb, Q_filtered_usb)   
    
    # Extract final IQ for lower sideband
    I_window_lsb, Q_window_lsb, I_time_lsb, I_full_lsb, Q_full_lsb, I_time_full_lsb = extract_finalIQ(I_filtered_lsb, Q_filtered_lsb) 
def and_process_TripleDDC(I, Q, Ip, Qp,card, settings, plot, IQstorage = True):

    # card.ArmAndWait()
    # I, Q = card.ReadAllData()
    # Ip = I #np.mean(I, 0)
    # Qp = Q #np.mean(Q, 0)
    
    
    mixer_config = settings['exp_globals']['mixer_config']
    Ipp, Qpp = remove_IQ_ellipse(Ip, Qp, mixer_config)
    xaxis = np.linspace(0, card.samples, card.samples, endpoint=False)/card.sampleRate
    
    I_filtered_carr, I_filtered_usb, I_filtered_lsb = extract_data_heterodyne_TripleDDC(Ipp, xaxis, settings)
    Q_filtered_carr, Q_filtered_usb, Q_filtered_lsb = extract_data_heterodyne_TripleDDC(Qpp, xaxis, settings)
    
    # Extract final IQ for carrier
    I_window_carr, Q_window_carr, I_time_carr, I_full_carr, Q_full_carr, I_time_full_carr = extract_finalIQ(I_filtered_carr, Q_filtered_carr)
        
    # Extract final IQ for upper sideband
    I_window_usb, Q_window_usb, I_time_usb, I_full_usb, Q_full_usb, I_time_full_usb = extract_finalIQ(I_filtered_usb, Q_filtered_usb)   
    
    # Extract final IQ for lower sideband
    I_window_lsb, Q_window_lsb, I_time_lsb, I_full_lsb, Q_full_lsb, I_time_full_lsb = extract_finalIQ(I_filtered_lsb, Q_filtered_lsb) 

    ### For debugging - plot fft of Ip and Qp 
    fft_mag_Ip = np.abs(np.fft.fft(Ip[0]))
    fft_freq_Ip = np.fft.fftfreq(len(Ip[0]),1/card.sampleRate)
    
    fft_mag_Qp = np.abs(np.fft.fft(Qp[0]))
    fft_freq_Qp = np.fft.fftfreq(len(Qp[0]),1/card.sampleRate)
    
    if plot:
        fig0 = plt.figure(971)
        plt.clf()
        ax = plt.subplot(2,2,1)
        plt.plot(xaxis*1e6, Ip[0], label = 'raw V1')
        plt.plot(xaxis*1e6, Qp[0], label = 'raw V2')
        plt.xlabel('Time (us)')
        plt.ylabel('Voltage')
        plt.title('Card Voltages')
        ax.legend(loc = 'upper right')
        
        ax = plt.subplot(2,2,2)
        plt.plot(xaxis*1e6, np.sqrt(Ip[0]**2 + Qp[0]**2), label = 'crude amplitude')
        plt.plot(xaxis*1e6, np.sqrt(Ipp[0]**2 + Qpp[0]**2), label = 'mixer corrected')
        plt.xlabel('Time (us)')
        plt.ylabel('Voltage')
        plt.title('Rough Pulse Amplitude')
        ax.legend(loc = 'upper right')
        
        ax = plt.subplot(2,2,3)
        plt.plot(fft_freq_Ip[:len(fft_freq_Ip)//2]/1e6, np.log10(fft_mag_Ip[:len(fft_mag_Ip)//2]), label='raw_I')
        plt.xlabel("Freq (MHz)")
        plt.xlim(0,50)
        plt.title('FFT for raw V1')
        ax.legend(loc = 'upper right')
        
        ax = plt.subplot(2,2,4)
        plt.plot(fft_freq_Qp[:len(fft_freq_Qp)//2]/1e6, np.log10(fft_mag_Qp[:len(fft_mag_Qp)//2]), label='raw_Q')
        plt.xlabel("Freq (MHz)")
        plt.xlim(0,50)
        plt.title('FFT for raw V2')
        ax.legend(loc = 'upper right')
        
        plt.show()
        fig0.canvas.draw()
        fig0.canvas.flush_events()
        
        carr = {}
        carr['amp'] = np.sqrt(Q_window_carr[0]**2 + I_window_carr[0]**2)
        carr['amp_full'] = np.sqrt(Q_full_carr[0]**2 + I_full_carr[0]**2)
        carr['Ifull'], carr['Qfull'] = I_full_carr[0], Q_full_carr[0] 
        
        usb = {}
        usb['amp'] = np.sqrt(Q_window_usb[0]**2 + I_window_usb[0]**2)
        usb['amp_full'] = np.sqrt(Q_full_usb[0]**2 + I_full_usb[0]**2)
        usb['Ifull'], usb['Qfull'] = I_full_usb[0], Q_full_usb[0] 
        
        lsb = {}
        lsb['amp'] = np.sqrt(Q_window_lsb[0]**2 + I_window_lsb[0]**2)
        lsb['amp_full'] = np.sqrt(Q_full_lsb[0]**2 + I_full_lsb[0]**2)
        lsb['Ifull'], lsb['Qfull'] = I_full_lsb[0], Q_full_lsb[0] 
        # carr = {}
        # carr['amp'] = np.sqrt(Q_window_carr**2 + I_window_carr**2)
        # carr['amp_full'] = np.sqrt(Q_full_carr**2 + I_full_carr**2)
        # carr['Ifull'], carr['Qfull'] = I_full_carr, Q_full_carr 
        
        # usb = {}
        # usb['amp'] = np.sqrt(Q_window_usb**2 + I_window_usb**2)
        # usb['amp_full'] = np.sqrt(Q_full_usb**2 + I_full_usb**2)
        # usb['Ifull'], usb['Qfull'] = I_full_usb, Q_full_usb 
        
        # lsb = {}
        # lsb['amp'] = np.sqrt(Q_window_lsb**2 + I_window_lsb**2)
        # lsb['amp_full'] = np.sqrt(Q_full_lsb**2 + I_full_lsb**2)
        # lsb['Ifull'], lsb['Qfull'] = I_full_lsb, Q_full_lsb 
        
        plot_data_extraction_TripleDDC(carr, usb, lsb, I_time_carr, I_time_full_carr)
        
    # save data for return
    carr_data = {}
    carr_data['I_window'], carr_data['Q_window'], carr_data['I_full'], carr_data['Q_full'] = I_window_carr, Q_window_carr, I_full_carr, Q_full_carr
    
    usb_data = {}
    usb_data['I_window'], usb_data['Q_window'], usb_data['I_full'], usb_data['Q_full'] = I_window_usb, Q_window_usb, I_full_usb, Q_full_usb
    
    lsb_data = {}
    lsb_data['I_window'], lsb_data['Q_window'], lsb_data['I_full'], lsb_data['Q_full'] = I_window_lsb, Q_window_lsb, I_full_lsb, Q_full_lsb
        
    return carr_data, usb_data, lsb_data, I_time_full_carr
    ### For debugging - plot fft of Ip and Qp 
    fft_mag_Ip = np.abs(np.fft.fft(Ip[0]))
    fft_freq_Ip = np.fft.fftfreq(len(Ip[0]),1/card.sampleRate)
    
    fft_mag_Qp = np.abs(np.fft.fft(Qp[0]))
    fft_freq_Qp = np.fft.fftfreq(len(Qp[0]),1/card.sampleRate)
    
    if plot:
        fig0 = plt.figure(971)
        plt.clf()
        ax = plt.subplot(2,2,1)
        plt.plot(xaxis*1e6, Ip[0], label = 'raw V1')
        plt.plot(xaxis*1e6, Qp[0], label = 'raw V2')
        plt.xlabel('Time (us)')
        plt.ylabel('Voltage')
        plt.title('Card Voltages')
        ax.legend(loc = 'upper right')
        
        ax = plt.subplot(2,2,2)
        plt.plot(xaxis*1e6, np.sqrt(Ip[0]**2 + Qp[0]**2), label = 'crude amplitude')
        plt.plot(xaxis*1e6, np.sqrt(Ipp[0]**2 + Qpp[0]**2), label = 'mixer corrected')
        plt.xlabel('Time (us)')
        plt.ylabel('Voltage')
        plt.title('Rough Pulse Amplitude')
        ax.legend(loc = 'upper right')
        
        ax = plt.subplot(2,2,3)
        plt.plot(fft_freq_Ip[:len(fft_freq_Ip)//2]/1e6, np.log10(fft_mag_Ip[:len(fft_mag_Ip)//2]), label='raw_I')
        plt.xlabel("Freq (MHz)")
        plt.xlim(0,50)
        plt.title('FFT for raw V1')
        ax.legend(loc = 'upper right')
        
        ax = plt.subplot(2,2,4)
        plt.plot(fft_freq_Qp[:len(fft_freq_Qp)//2]/1e6, np.log10(fft_mag_Qp[:len(fft_mag_Qp)//2]), label='raw_Q')
        plt.xlabel("Freq (MHz)")
        plt.xlim(0,50)
        plt.title('FFT for raw V2')
        ax.legend(loc = 'upper right')
        
        plt.show()
        fig0.canvas.draw()
        fig0.canvas.flush_events()
        
        carr = {}
        carr['amp'] = np.sqrt(Q_window_carr[0]**2 + I_window_carr[0]**2)
        carr['amp_full'] = np.sqrt(Q_full_carr[0]**2 + I_full_carr[0]**2)
        carr['Ifull'], carr['Qfull'] = I_full_carr[0], Q_full_carr[0] 
        
        usb = {}
        usb['amp'] = np.sqrt(Q_window_usb[0]**2 + I_window_usb[0]**2)
        usb['amp_full'] = np.sqrt(Q_full_usb[0]**2 + I_full_usb[0]**2)
        usb['Ifull'], usb['Qfull'] = I_full_usb[0], Q_full_usb[0] 
        
        lsb = {}
        lsb['amp'] = np.sqrt(Q_window_lsb[0]**2 + I_window_lsb[0]**2)
        lsb['amp_full'] = np.sqrt(Q_full_lsb[0]**2 + I_full_lsb[0]**2)
        lsb['Ifull'], lsb['Qfull'] = I_full_lsb[0], Q_full_lsb[0] 
        # carr = {}
        # carr['amp'] = np.sqrt(Q_window_carr**2 + I_window_carr**2)
        # carr['amp_full'] = np.sqrt(Q_full_carr**2 + I_full_carr**2)
        # carr['Ifull'], carr['Qfull'] = I_full_carr, Q_full_carr 
        
        # usb = {}
        # usb['amp'] = np.sqrt(Q_window_usb**2 + I_window_usb**2)
        # usb['amp_full'] = np.sqrt(Q_full_usb**2 + I_full_usb**2)
        # usb['Ifull'], usb['Qfull'] = I_full_usb, Q_full_usb 
        
        # lsb = {}
        # lsb['amp'] = np.sqrt(Q_window_lsb**2 + I_window_lsb**2)
        # lsb['amp_full'] = np.sqrt(Q_full_lsb**2 + I_full_lsb**2)
        # lsb['Ifull'], lsb['Qfull'] = I_full_lsb, Q_full_lsb 
        
        plot_data_extraction_TripleDDC(carr, usb, lsb, I_time_carr, I_time_full_carr)
        
    # save data for return
    carr_data = {}
    carr_data['I_window'], carr_data['Q_window'], carr_data['I_full'], carr_data['Q_full'] = I_window_carr, Q_window_carr, I_full_carr, Q_full_carr
    
    usb_data = {}
    usb_data['I_window'], usb_data['Q_window'], usb_data['I_full'], usb_data['Q_full'] = I_window_usb, Q_window_usb, I_full_usb, Q_full_usb
    
    lsb_data = {}
    lsb_data['I_window'], lsb_data['Q_window'], lsb_data['I_full'], lsb_data['Q_full'] = I_window_lsb, Q_window_lsb, I_full_lsb, Q_full_lsb
        
    return carr_data, usb_data, lsb_data, I_time_full_carr




