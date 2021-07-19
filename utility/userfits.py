# -*- coding: utf-8 -*-
"""
Created on Tue Dec  8 15:27:44 2020

@author: Kollarlab
"""
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from lmfit.models import LorentzianModel, ConstantModel

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
    plt.figure(62)
    plt.clf()
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
    plt.figure(62)
    plt.clf()
    plt.plot(taus*1e6, amps, 'x')
    plt.plot(ts*1e6, fit_curve)
    
    plt.title('T2 fit, T2: {}us, detuning: {} kHz'.format(round(fit_out[0]*1e6, 2), round(fit_out[3]/1e3, 3)))
    plt.xlabel('Tau (us)')
    plt.ylabel('Amp')

    return fit_out[0], fit_out[3]

def convert_to_linear(yvals):
    linear = 10**(yvals/10)  
    return linear

def fit_lorentzian(freqs, mags, plot=False, fig_num=None):
    line = convert_to_linear(mags)
    data_test = -(line-max(line))
    
    peak = LorentzianModel()
    const = ConstantModel()
    mod = peak + const
    
    center = 0
    sigma = 0
    amp   = 0
    try:
        pars = peak.guess(data_test, x=freqs)
        pars += const.make_params()
        
        out = mod.fit(data_test, pars, x=freqs)
        fit_params = out.best_values
        center = fit_params['center']
        sigma  = fit_params['sigma']
        amp    = fit_params['amplitude']/(np.pi*sigma)
        model  = mod.eval(params=out.params, x=freqs)
        amp_dB = 10*np.log10(min(model)/max(model))
    except:
        print('Fit did not converge, returning zeros')
        return center, sigma
    
    if plot:
        if fig_num is None:
            plt.figure()
        else:
            plt.figure(fig_num)
        plt.clf()
        plt.plot(freqs, data_test)
        plt.plot(freqs, model)
        plt.show()
        
    return center, 2*sigma, amp_dB
