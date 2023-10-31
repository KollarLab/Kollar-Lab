# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 18:12:10 2023

@author: Kollarlab
"""

## Quick attempt at unit conversion stuff 
from bidict import bidict
import numpy as np

SI_prefix = bidict({
        'p': -12,
        'n': -9,
        'u': -6,
        'm': -3,
        '' : 0,
        'k': 3,
        'M': 6,
        'G': 9
        })

def set_prefix(val, prefix):
    return val*10**SI_prefix[prefix]

def get_prefix(val, sigfigs=5):
    temp = val
    sign = -1
    exponent = 0
    if val>1:
        temp = 1./val
        sign = 1
        exponent = 0
    while temp<1:
        temp = temp*1e3
        exponent = exponent+3
    if val>1:
        exponent -= 3
    exponent = sign*exponent
    return np.round(val/10**exponent, sigfigs), SI_prefix.inverse[exponent]

def convert_to_prefix(vals, prefix):
    return np.array(vals)/10**SI_prefix[prefix]
            