# -*- coding: utf-8 -*-
"""
Created on Tue Dec  8 15:27:44 2020

@author: Kollarlab
"""
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def expff(x, tau, amp, offset):
    vals = amp*np.exp(- (x)/tau) + offset
    return vals

def exp_cosff(x, tau, amp, offset, freq, phi):
    pi = np.pi
    vals = amp*np.cos(2*pi*freq*x+phi)*np.exp(- (x)/tau) + offset
    return vals

def fit_T1(taus, amps, fit_guess):
    ts = np.linspace(min(taus), max(taus), 10*len(taus))
    
    fit_out, pcov = curve_fit(expff, taus, amps, p0 = fit_guess)
    
    fit_curve = expff(ts, fit_out[0], fit_out[1], fit_out[2])
    plt.plot(taus, amps, 'x')
    plt.plot(ts, fit_curve)
    
    plt.title('T1 fit, T1: {}us'.format(round(fit_out[0], 2)))
    plt.xlabel('Tau (s)')
    plt.ylabel('Amp')
    
    return fit_out[0]

def fit_T2(taus, amps, fit_guess):
    ts = np.linspace(min(taus), max(taus), 10*len(taus))
    
    fit_out, pcov = curve_fit(exp_cosff, taus, amps, p0 = fit_guess)
    
    fit_curve = exp_cosff(ts, fit_out[0], fit_out[1], fit_out[2], fit_out[3], fit_out[4])
    plt.plot(taus*1e6, amps, 'x')
    plt.plot(ts*1e6, fit_curve)
    
    plt.title('T2 fit, T2: {}us, detuning: {} kHz'.format(round(fit_out[0]*1e6, 2), round(fit_out[3]/1e3, 3)))
    plt.xlabel('Tau (us)')
    plt.ylabel('Amp')
    
    return fit_out[0], fit_out[3]