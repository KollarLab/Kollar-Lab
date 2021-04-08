# -*- coding: utf-8 -*-
"""
Created on Fri Apr  2 10:40:36 2021

@author: Kollarlab
"""
import userfuncs
import os
import numpy as np
from lmfit.models import LorentzianModel, ConstantModel
import matplotlib.pyplot as plt

def convert_to_linear(yvals):
    linear = np.zeros(yvals.shape)
    for index, val in enumerate(yvals):
        linear[index] = 10**(val/20)
    
    return linear

saveDir = r'Z:\Data\Fluxonium_Raman\CRF01_A3\Spec\20210406'
filename = 'fluxon_0mV_line_drift_20210406_195005.pkl'
dat = userfuncs.LoadFull(os.path.join(saveDir, filename))

data = dat[0]['full_data']
xaxis = data['xaxis']
yvals = data['mags']

linear = convert_to_linear(yvals)

centers = np.zeros(len(linear))
sigmas = np.zeros(len(linear))

for index, line in enumerate(linear):
    data_test = -(line-max(line))
#    plt.figure()
#    plt.plot(xaxis, data_test)
    
    peak = LorentzianModel()
    const = ConstantModel()
    mod = peak + const
    
    pars = peak.guess(data_test, x=xaxis)
    pars += const.make_params()
    
    out = mod.fit(data_test, pars, x=xaxis)
    fit_params = out.best_values
    centers[index] = fit_params['center']
    sigmas[index]  = fit_params['sigma']
    
#    plt.plot(xaxis, mod.eval(params=out.params, x=xaxis))

powers = dat[0]['powers'] - dat[1]['Qbit_Attenuation']
plt.figure()
plt.subplot(1,2,1)
plt.plot(powers, centers/1e9)
plt.title('Center freq')
plt.xlabel('Power (dBm)')
plt.ylabel('Freq (GHz)')

plt.subplot(1,2,2)
plt.plot(powers, sigmas/1e6)
plt.title('sigma (fit)')
plt.xlabel('Power (dBm)')
plt.ylabel('Freq (MHz)')

plt.suptitle('line drift')