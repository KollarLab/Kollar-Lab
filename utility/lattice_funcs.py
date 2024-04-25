# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 15:48:20 2024

@author: kollarlab
"""

import numpy as np
from scipy.optimize import curve_fit
import userfuncs
from utility.userfits_v2 import fit_model
import os


################################################################################
def full_data_fit(spec_file,hanger = False):
    #spec_file: str, filename of relevant data
    #hanger: bool, whether the device is a hanger
    #Loads spec flux data, finds max/min for trans and spec values to extract cavity/qubit frequencies at each voltage point.
    full_data = userfuncs.LoadFull(spec_file)

    if 'powers' in full_data[0]:
        xaxis = full_data[0]['powers']
    elif 'voltages' in full_data[0]:
        xaxis = full_data[0]['voltages']
    elif 'fluxes' in full_data[0]:
        xaxis = full_data[0]['fluxes']


    trans_freqs_GHz = full_data[0]['transdata']['xaxis']
    trans_mags = full_data[0]['transdata']['phases']
    
    cfreqs = []
    
    for f in range(0,len(xaxis)):
        if hanger:
            fholder = trans_freqs_GHz[np.argmin(trans_mags[f])]
        else:
            fholder = trans_freqs_GHz[np.argmax(trans_mags[f])]
        cfreqs.append(fholder)
    
    spec_freqs_GHz = full_data[0]['specdata']['xaxis']
    spec_mags = full_data[0]['specdata']['phases']
    
    qfreqs = []
    
    for f in range(0,len(xaxis)):
        if hanger:
            qholder = spec_freqs_GHz[np.argmax(spec_mags[f])]
        else:
            qholder = spec_freqs_GHz[np.argmin(spec_mags[f])]
        qfreqs.append(qholder)
    
    
    return xaxis, cfreqs, qfreqs



def full_data_fit_mags(spec_file,hanger = False):
    #spec_file: str, filename of relevant data
    #hanger: bool, whether the device is a hanger
    #Loads spec flux data, finds max/min for trans and spec values to extract cavity/qubit frequencies at each voltage point.
    full_data = userfuncs.LoadFull(spec_file)
    
    if 'powers' in full_data[0]:
        xaxis = full_data[0]['powers']
    elif 'voltages' in full_data[0]:
        xaxis = full_data[0]['voltages']
    elif 'fluxes' in full_data[0]:
        xaxis = full_data[0]['fluxes']
    

    trans_freqs_GHz = full_data[0]['transdata']['xaxis']
    trans_mags = full_data[0]['transdata']['phases']
    
    cfreqs = []
    
    for f in range(0,len(xaxis)):
        if hanger:
            fholder = trans_freqs_GHz[np.argmin(trans_mags[f])]
        else:
            fholder = trans_freqs_GHz[np.argmax(trans_mags[f])]
        cfreqs.append(fholder)
    
    spec_freqs_GHz = full_data[0]['specdata']['xaxis']
    spec_mags = full_data[0]['specdata']['mags']
    
    qfreqs = []
    
    for f in range(0,len(xaxis)):
        if hanger:
            qholder = spec_freqs_GHz[np.argmax(spec_mags[f])]
        else:
            qholder = spec_freqs_GHz[np.argmin(spec_mags[f])]
        qfreqs.append(qholder)
    
    
    return xaxis, cfreqs, qfreqs


def trans_data_fit(spec_file,hanger = True):
    #spec_file: str, filename of relevant data
    #hanger: bool, whether the device is a hanger
    #Loads trans flux data, finds max/min for trans values to extract cavity frequencies at each voltage point.
    full_data = userfuncs.LoadFull(spec_file)

    if 'powers' in full_data[0]:
        xaxis = full_data[0]['powers']
    elif 'voltages' in full_data[0]:
        xaxis = full_data[0]['voltages']
    elif 'fluxes' in full_data[0]:
        xaxis = full_data[0]['fluxes']
    
    trans_freqs_GHz = full_data[0]['full_data']['xaxis']
    trans_mags = full_data[0]['full_data']['phases']
    
    cfreqs = []
    
    for f in range(0,len(xaxis)):
        if hanger:
            fholder = trans_freqs_GHz[np.argmin(trans_mags[f])]
        else:
            fholder = trans_freqs_GHz[np.argmax(trans_mags[f])]
        cfreqs.append(fholder)
    
    return xaxis, cfreqs

################################################################################

def lin_fun(x,m,b):
    #x, m, b: float
    #You get it
    return m*x + b

def kerr_fit(data_trans, data_transoe, data_spec, RF_atten, index=False):
    if index:
        data_trans_mags = data_trans['full_data']['mags'][index]
        data_trans_phases = data_trans['full_data']['phases'][index]
        data_transoe_phases = data_transoe['full_data']['phases'][index]
    else:
        data_trans_mags = data_trans['full_data']['mags'][0]
        data_trans_phases = data_trans['full_data']['phases'][0]
        data_transoe_phases = data_transoe['full_data']['phases'][0]
    
    
    zero_point_ind = np.argmax(data_trans_mags)

    data_spec_phase = data_spec['full_data']['phases']
    data_spec_phase_min = np.zeros(len(data_spec_phase))
    for i in range(len(data_spec_phase)):
        data_spec_phase_min[i] = np.min(data_spec_phase[i])

    lin_ax = 10**((data_spec['powers']-RF_atten)/10)
    
    normalize = data_transoe_phases[zero_point_ind] - data_trans_phases[zero_point_ind]
    normalized_phases = data_transoe_phases - normalize
    normalized_phases = np.unwrap(normalized_phases,period=360)
    
    if normalized_phases[zero_point_ind]>180:
        normalized_phases = normalized_phases - 360
    elif normalized_phases[zero_point_ind]<-180:
        normalized_phases = normalized_phases + 360
    
    monitor_tone = np.zeros(len(data_spec_phase))
    for i in range(len(data_spec_phase_min)):
        a = np.argmin(np.abs(normalized_phases - data_spec_phase_min[i]))
        monitor_tone[i] = data_transoe['full_data']['xaxis'][zero_point_ind] - data_transoe['full_data']['xaxis'][a] 
    
    popt, pcov = curve_fit(lin_fun, lin_ax, monitor_tone, method='lm', maxfev=10000)
    # txtstr = 'fit : ' + str(int(popt[0])) + '[KHz/mW]'
    
    
    output_dict = {'power_mW': lin_ax, 'phase_shift' : data_spec_phase_min, 
                   'freq_shift' : monitor_tone,'slope_GHz_mW' : popt[0], 'full_fit':popt}
    return output_dict

def full_flux_vector(flux_start,flux_stop,flux_pts):
    #flux_start: array of length n
    #flux_start: array of length n
    #flux_pts: integer
    #Function for generating a series of flux vectors to feed to a flux scan
    flux_holder = []
    if len(flux_start) != len(flux_stop):
        raise ValueError('requesting invalid array of flux points')
    for i in range(len(flux_start)):
        flux_range = np.linspace(flux_start[i],flux_stop[i],flux_pts)
        flux_holder.append(flux_range)

    flux_holder = tuple(flux_holder)
    full_fluxes = np.stack(flux_holder,axis=1)

    return full_fluxes    

def flux2v_generator(v2f,v_offsets,full_fluxes):
    #Function for calculating the corresponding voltage vectors for an 
    #array of flux vectors
    f2v = np.linalg.inv(v2f)
    diags = np.diagonal(v2f)
    phase_offsets = v_offsets * (diags)
    
    full_voltages = np.zeros(full_fluxes.shape)
    
    for i in range(len(full_fluxes)):
        desired_phases = full_fluxes[i] + phase_offsets
        full_voltages[i] = f2v@desired_phases
    return full_voltages




def phase_fun(volts,v2f,v_offsets):
    #volts: (n,) array of voltages
    #v2f: (n,n) array, the volt to flux matrix
    #v_offsets: (n,) array, collection of voltage offsets
    #Transforms SRS voltages into the actual phase response of the qubits
    diags = np.diagonal(v2f) # Diagonal elements of the volt to flux matrix
    phase_offsets = v_offsets * (diags)
    phases = v2f@volts - phase_offsets
    return phases


def phase_finder(volt,v2f,v_offsets,SRS_ind):
    #volt: float, voltage point from a 1D SRS sweep you'd like to know the qubit's phase for
    #v2f: (n,n) array, the volt to flux matrix
    #v_offsets: (n,) array, collection of voltage offsets
    #SRS_ind: int, labels the qubit being driven
    #Finds the phase for a qubit when one SRS is at a voltage point.
    diags = np.diagonal(v2f)
    return diags[SRS_ind-1] * (volt - v_offsets[SRS_ind-1])

def volt_finder(phase,v2f,v_offsets,SRS_ind):
    #phase: float, qubit phase you'd like to target
    #v2f: (n,n) array, the volt to flux matrix
    #v_offsets: (n,) array, collection of voltage offsets
    #SRS_ind: int, labels the qubit being driven
    #Finds the phase for a qubit when one SRS is at a voltage point.
    diags = np.diagonal(v2f)
    return (1/diags[SRS_ind-1])*phase + v_offsets[SRS_ind-1]

# For qubit frequency calibration script
def calibration_slope(file_path, qubit_num):
    full_data = userfuncs.LoadFull(file_path+'.pkl')
    freq_list = []
    flux_list = []
    flux_num = len(full_data[0]['full_fluxes'])
    for i in range(flux_num):
        result = fit_model(full_data[0]['specdata']['xaxis'], full_data[0]['specdata']['mags'][i], 'lorenz')
        freq_list.append(result['center'])
        flux_list.append(full_data[0]['full_fluxes'][i][qubit_num-1])
    popt, pcov = curve_fit(linear_func, freq_list, flux_list, method='lm')
    slope = popt[0]
    return freq_list, flux_list, slope

def linear_func(x, a, b):
    return a*x + b

def calibration_flux(file_path, ref_freq_list, slope, qubit_num):
    full_data = userfuncs.LoadFull(file_path+'.pkl')
    del_f = []
    del_phi = []
    flux_list = full_data[0]['fluxes']
    flux_pts = len(flux_list)
    for i in range(flux_pts):
        result = fit_model(full_data[0]['specdata']['xaxis'], full_data[0]['specdata']['mags'][i], 'lorenz')
        del_f.append(ref_freq_list[i] - result['center'])
        del_phi.append(slope*(del_f[i]))
    avg_del_phi = sum(del_phi)/flux_pts
    avg_del_f = sum(del_f)/flux_pts
    new_flux_start = full_data[0]['full_fluxes'][0][qubit_num-1]+avg_del_phi
    new_flux_stop = full_data[0]['full_fluxes'][flux_pts-1][qubit_num-1]+avg_del_phi
    print(qubit_num-1)
    cali_flux_start = full_data[0]['full_fluxes'][0]
    cali_flux_start[qubit_num-1] = new_flux_start
    cali_flux_stop = full_data[0]['full_fluxes'][flux_pts-1]
    cali_flux_stop[qubit_num-1] = new_flux_stop
    return cali_flux_start, cali_flux_stop, avg_del_f