# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 15:34:34 2021

@author: Kollarlab
"""

import numpy as np
from scipy.optimize import minimize

def least_squares(energies, anchors):
    return np.square(energies-anchors).sum()
 
def fit_func(EC, EJ, EL, anchors, flux_array):
    (evals_s, _) = flux_sweep(flux_array*2*np.pi, EC, EJ, EL)
    energies = evals_s.T
    dist = least_squares(energies, anchors)
    
flux_vals = np.array([])
g2_vals = np.array([])
e1_vals = np.array([])

anchors = np.array(g2_vals, e1_vals)