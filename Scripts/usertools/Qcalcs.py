# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 16:04:59 2020

@author: Kollarlab
"""
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize
import scipy.signal
import scipy.constants

def fit_func(omega,omega_0,delta_omega,Q_i,Q_c,dB_offset):
    Q_0 = 1/(1/Q_i + 1/Q_c)
    num = Q_0/Q_c - 2*1j*Q_0*delta_omega/omega_0
    den = 1 + 2*1j*Q_0*(omega-omega_0)/omega_0
    return 20*np.log10(np.abs(1 - num/den)) + dB_offset

def fit_Q(freqs, mag, center):
    xdata = 2*np.pi*freqs
    ydata  = mag
    
    #f0_guess = start+span/2
    f0_guess = center
    
    popt, pcov = scipy.optimize.curve_fit(fit_func, xdata, ydata,p0=(2*np.pi*f0_guess,\
                                                                     2*np.pi*1e3,1e4,5e3,np.average(ydata)))
    fit_omega_0 = popt[0]
    print("Center Freq (GHz)",fit_omega_0/(2*np.pi*1e9))
    fit_delta_omega = popt[1]
    print("Assymetric (kHz)",fit_delta_omega/(2*np.pi*1e3))
    fit_Q_i = popt[2]
    print("Internal Quality Factor",fit_Q_i)
    fit_Q_c = popt[3]
    print("Coupling Quality Factor",fit_Q_c)
    # Q_c = pi/(2 q_i**2) where q_i = omega C Z0
    # assume Z0 = 50 Ohm
    coupling_C = np.sqrt(np.pi/(2*fit_Q_c))*(1/(fit_omega_0*50))
    print("Coupling capacitor (fF)",coupling_C*1e15)
    Q_loaded = (1/(1/fit_Q_i + 1/fit_Q_c))
    print("Loaded Quality Factor",Q_loaded)
    
    plt.figure(dpi=120)
    plt.plot((xdata-fit_omega_0)/(2*np.pi*1e6),ydata,color='b')
    plt.plot((xdata-fit_omega_0)/(2*np.pi*1e6),fit_func(xdata,*popt),color='r')
    plt.xlabel("Frequency (MHz)")
    plt.ylabel(r"$S_{21}$ (dB)")
    plt.title('F = {}GHz'.format(center/1e9))
    plt.grid(True)
    
    plt.show()
    return fit_omega_0, fit_Q_i, Q_loaded