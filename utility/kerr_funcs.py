# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 15:48:20 2024

@author: kollarlab
"""

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import userfuncs
from numpy.polynomial import Polynomial
from utility.userfits_v2 import fit_model, lorenztian_model
from utility.plotting_tools import general_colormap_subplot
import os


################################################################################
#Generic stuff

def lin_fun(x,m,b):
    '''
    lin_fun _summary_

    :param x: _description_
    :type x: _type_
    :param m: _description_
    :type m: _type_
    :param b: _description_
    :type b: _type_
    :return: _description_
    :rtype: _type_
    '''    
    #x, m, b: float
    #You get it
    return m*x + b

def tan_model(y,sigma,offset):
    x = sigma*np.tan((y-offset)*(-1*np.pi/180))
    return(x)

def arctan_model(f_ax,f_c,sigma, offset):
    y = -(180/np.pi) * np.arctan((f_ax-f_c)/sigma) + offset
    return y

def semiclassic(x,x0,E,kappa,nu):
    y = -(E**2/((x-x0)**2 + kappa**2/4))*(1-(2*nu*(x-x0)*E**2)/(((x-x0)**2 + kappa**2/4)**2) - (nu**2*E**4)/(((x-x0)**2 + kappa**2/4)**3))
    return y


##################################################################################
#Fittings to Kerr data
def stark_fit_mod(data_file,save_fig = True,show_fig = False,plot=True):
    '''
    Parameters
    ----------
    data_file : str
        pkl file with stark scan data set
    save_fig : bool
        Whether you want the fitted data set to be saved at the end

    Returns
    -------
    Slope + intercept of fitted AC Stark data set

    '''
    dataset = userfuncs.LoadFull(data_file)

    spec_dat =  dataset[0]['specdata']
    spec_freqs = spec_dat['xaxis']

    stark_pows = dataset[0]['powers']
    lin_pows = 10**((stark_pows-dataset[1]['exp_globals']['CAV_Attenuation'])/10)
    lin_mags = 10**(spec_dat['mags']/20) #Martin says fit function works better if voltages are used instead of pows

    centers = np.zeros(len(stark_pows))
    full_results = []
    

    for pind in range(len(stark_pows)):
        results = fit_model(1e9*spec_freqs,lin_mags[pind],'gauss',plot=show_fig)
        centers[pind] = results['center']/1e9
        full_results.append(results)

    m_guess = (centers[-1] - centers[0])/(lin_pows[-1]-lin_pows[0])
    fit_out1, pcov = curve_fit(lin_fun,lin_pows,centers,p0 = [m_guess,centers[0]])

    if plot:
        plt.figure(5)
        plt.clf()
        plt.plot(lin_pows,centers,'x',label='Data')
        plt.plot(lin_pows,lin_fun(lin_pows,fit_out1[0],fit_out1[1]))
        plt.xlabel('VNA Power (mW)')
        plt.ylabel('Qubit Frequency (GHz)')
        plt.title('Stark Shift Fitting, m=' + str(np.round(fit_out1[0],3))+ ' GHz/mW: ' + data_file[-4])
    
    if save_fig:
        plt.savefig(os.path.join(data_file[:-4] + '_Fitting.png'),dpi=150)
        
    #qfreq = centers[0]
    
    return fit_out1, pcov, centers, full_results, lin_pows


def arctan_kerr(data_trans, data_transoe, data_spec, RF_atten, index=False, plot=True): #This is the final version for analysis
    '''
    kerr_fit _summary_

    :param data_trans: transmission data set
    :type data_trans: dictionary
    :param data_transoe: transmission data set with offset embed correction
    :type data_transoe: dictionary
    :param data_spec: dictionary
    :type data_spec: kerr spec data
    :param RF_atten: attenuation on the drive line
    :type RF_atten: int
    :param index: if extracting from data sets with multiple flux points, index picks out which flux point we are working with, defaults to False
    :type index: bool
    :return: dictionary with key values (various fit parameters)
    :rtype: dict
    '''    
    if index:
        data_trans_mags = data_trans['full_data']['mags'][index]
        data_trans_phases = data_trans['full_data']['phases'][index]
        data_transoe_phases = data_transoe['full_data']['phases'][index]
    else:
        data_trans_mags = data_trans['full_data']['mags'][0]
        data_trans_phases = data_trans['full_data']['phases'][0]
        data_transoe_phases = data_transoe['full_data']['phases'][0]
    
    data_trans_freqs = data_trans['full_data']['xaxis']
    
    zero_point_ind = np.argmax(data_trans_mags)
 
    #if plot==True:
        #path = os.path.join(data_spec['saveDir'],data_spec['filename'])
        #plt.savefig(path+'_Polynomial_Fits.png')

    
    normalize = data_transoe_phases[zero_point_ind] - data_trans_phases[zero_point_ind]
    print(data_transoe_phases[zero_point_ind])
    print(data_trans_phases[zero_point_ind])
    normalized_phases = data_transoe_phases - normalize
    normalized_phases = np.unwrap(normalized_phases,period=360)
    
    if normalized_phases[zero_point_ind]>180:
        normalized_phases = normalized_phases - 360
    elif normalized_phases[zero_point_ind]<-180:
        normalized_phases = normalized_phases + 360
    
    #insert arctan fitting

    
    p0 = [data_trans_freqs[0], 0.0001, normalized_phases[zero_point_ind]]
    
    xaxis_sub = data_trans_freqs[180:300] #This is hard coded and bad
    phase_sub = normalized_phases[180:300]
    
    
    arc_co, blank = curve_fit(arctan_model, xaxis_sub, phase_sub, method='lm', maxfev=10000, p0= p0)
    
    if plot==True:
        plt.figure()
        plt.plot(data_trans_freqs,normalized_phases)
        plt.plot(data_trans_freqs,arctan_model(data_trans_freqs,arc_co[0],arc_co[1],arc_co[2]))
        plt.title('arctan Fit Check')
        
        plt.figure()
        plt.plot(data_trans_freqs,data_trans_phases,label='uncorrected')
        plt.plot(data_trans_freqs,data_transoe_phases,label='corrected')
        plt.legend()
        

    
    # output_dict = {'power_mW': lin_ax, 'phase_shift' : data_spec_phase_min, 
    #                'freq_shift' : monitor_tone,'slope_GHz_mW' : popt[0], 'full_fit':popt,
    #                'sk_freqs' : sk_freqs, 'self_Kerr_GHz_mW' : popt_sk[0], 'full_sk_fit' : popt_sk}
    return arc_co, data_trans_phases[zero_point_ind]

def recursive_mod(deltas,delta0,E,kappa,nu):
    ns = np.zeros(len(deltas))
    full_roots = []
    for d in range(len(deltas)):
        coeffs = [nu**2,2*nu*(deltas[d]-delta0),(deltas[d]-delta0)**2+kappa**2/4,-1*E**2]
        roots = np.roots(coeffs)
        
        real_roots = roots[np.isreal(roots)].real
        
        full_roots.append(real_roots)
        
        ns[d] = real_roots[0]
        # if len(real_roots) > 1:
            # print('More than one solution -- be careful')
            # print(real_roots)
            # ns[d] = real_roots[1]
    ns = -1*ns
    return ns

def g_square_gamma(alpha,qfreq,modefreq,stark_slope):
    y = stark_slope*(qfreq - modefreq)*(alpha - qfreq + modefreq)/(2*alpha)
    return y

def g_square_kerr(gg_gamma,alpha,qfreq,mode1,mode2,mkerr):
    d1 = qfreq - mode1
    d2 = qfreq - mode2
    y = mkerr*((d1**2) * (d2**2) *  (alpha - (d1+d2)))/(2*(alpha*(d1+d2)*gg_gamma))
    return y

def g_square_sk(gg_gamma,alpha,qfreq,mode,skerr):
    d = qfreq - mode
    y = (skerr * d**3 * (alpha - 2*d))/(2*alpha*gg_gamma)
    return y

def g_square_sk_v2(gg_gamma,alpha,qfreq,mode,skerr):
    d = qfreq - mode
    y = (skerr * d**3 * (alpha - 2*d))/(alpha*gg_gamma)
    return y



