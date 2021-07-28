# -*- coding: utf-8 -*-
"""
Created on Tue Dec  8 15:27:44 2020

@author: Kollarlab
"""
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from lmfit.models import LorentzianModel, ConstantModel

def T1_model(x, tau, amp, offset):
    return amp*np.exp(- (x)/tau) + offset

def T2_model(x, tau, amp, offset, freq, phi):
    pi = np.pi
    return amp*np.cos(2*pi*freq*x+phi)*np.exp(- (x)/tau) + offset

def fit_T1(taus, amps, fit_guess):
    ts = np.linspace(min(taus), max(taus), 10*len(taus))
    fit_curve = np.zeros(len(ts))

    tau = 0
    amp = 0
    offset = 0

    try:
        fit_out, pcov = curve_fit(T1_model, taus, amps, p0 = fit_guess)
        tau, amp, offset, *_ = fit_out
        fit_curve = T1_model(ts, tau, amp, offset)
    except:
        print('Fit did not converge, returning zeros')
        return tau, amp, offset, ts, fit_curve

    return tau, amp, offset, ts, fit_curve

def fit_T2(taus, amps, fit_guess):
    ts = np.linspace(min(taus), max(taus), 10*len(taus))
    fit_curve = np.zeros(len(ts))

    tau = 0
    amp = 0
    offset = 0
    freq = 0
    phi = 0

    try:
        fit_out, pcov = curve_fit(T2_model, taus, amps, p0=fit_guess)
        tau, amp, offset, freq, phi, *_ = fit_out
        fit_curve = T2_model(ts, tau, amp, offset, freq, phi)   
    except:
        print('Fit did not converge, returning zeros')
        return tau, amp, offset, freq, phi, ts, fit_curve 
    
    return tau, amp, offset, freq, phi, ts, fit_curve 

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
