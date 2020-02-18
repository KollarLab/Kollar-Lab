#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 09:53:45 2018

@author: kollar2
"""

import re
import scipy
import pylab
import numpy
import time

import pickle
import datetime
import os
import sys

import scipy.io as sio
from scipy import signal


from LayoutGenerator5 import PlanarLayout




#####
#settings and data set
#####asy
##########
#load actual data
##########


#ModeFamily = 'HW'
ModeFamily = 'FW'


numSamples = 15


##############################

if ModeFamily == 'FW':
    ####FW C4
    file1Path = r'/volumes/ourphoton/Alicia/Data/Hyperbolic_07_01/pnaxDownConversionSpectra121018/pnaxDownconversionHomodyne_FW_fullRes_C4_20181210_0919.mat'
#    file1Path = r'/volumes/ourphoton/Alicia/Data/Hyperbolic_07_01/pnaxDownConversionSpectra121018/pnaxDownconversionHomodyne_FW_fullRes_C1_20181210_1334.mat'
    Data1Struct = sio.loadmat(file1Path ,struct_as_record=False,squeeze_me=True);
    transmissionVec1_dB = Data1Struct['transmissionSpectrum']
    transmissionFreqs1 = Data1Struct['config'].driveFreqs
    ###temp oevrrride
    transmissionVec1_dB = transmissionVec1_dB-68
    #####
    transmissionVec1 = 10.**(transmissionVec1_dB/10)
    
    
    
#    ####FW  line background
#    file2Path = r'/volumes/ourphoton/Alicia/Data/Hyperbolic_03_02/phaseMeasurements032718/pnaxDownconversionHomodyne_phase_pnaxAcquired_GenDownconversionSpectrum_phase_FW_input_C1_AbeMixer_15K_20180327_0938.mat'
#    Data2Struct = sio.loadmat(file2Path ,struct_as_record=False,squeeze_me=True);
#    transmissionVec2_dB = Data2Struct['transmissionSpectrum']
#    transmissionFreqs2 = Data2Struct['config'].driveFreqs
#    ###temp oevrrride
#    transmissionVec2_dB = transmissionVec2_dB-60
#    #####
#    transmissionVec2 = 10.**(transmissionVec2_dB/10)
#    ######
    
elif ModeFamily == 'HW':
    ####HW B3
    file1Path = r'/volumes/ourphoton/Alicia/Data/Hyperbolic_07_01/directPNAX102618/hyperbolic07_01_B3_20181026_104324.mat'
    Data1Struct = sio.loadmat(file1Path ,struct_as_record=False,squeeze_me=True);
    transmissionVec1_dB = Data1Struct['S21amp']
    transmissionFreqs1 = Data1Struct['S21freqvector']
    ###temp oevrrride
    transmissionVec1_dB = transmissionVec1_dB-40
    #####
    transmissionVec1 = 10.**(transmissionVec1_dB/10)
    ######
    
    
    
#    ####HW B1
#    file2Path = r'/volumes/ourphoton/Alicia/Data/Hyperbolic_07_01/directPNAX102618/hyperbolic07_01_B1_20181026_105403.mat'
#    Data2Struct = sio.loadmat(file2Path ,struct_as_record=False,squeeze_me=True);
#    transmissionVec2_dB = Data2Struct['S21amp']
#    transmissionFreqs2 = Data2Struct['S21freqvector']
#    ###temp oevrrride
#    transmissionVec2_dB = transmissionVec2_dB-40
#    #####
#    transmissionVec2 = 10.**(transmissionVec2_dB/10)
#    ######
else:
    raise ValueError, 'invalid ModeType'

##############################






C1 = 7.28
C2 = 11.73
midFraction = 0.5  #ratio of symmetry and antisymmetry, I'm not sure I've set this quantity right at all
midFraction = C1**2/C2**2
#midFraction = C1/C2
#midFraction = 0.1

###fake symmetric couplers
midFraction = 1



#####
#initial guess at spectrum
#####





if ModeFamily == 'HW':
    
    ##option22, HW parameters, for hyperbolic 0701
    ##under development
#    JJ = -67.0*10**6
#    omega0 = (8.025*10**9 - 2*76.*10**6) - 17*10**6 - 2*JJ
    JJ = -66.0*10**6
    omega0 = (8.025*10**9 - 2*76.*10**6) - 17*10**6 - 2*JJ
    kappa = 0.7*10**6
#    sigma_diag = 6*10**(-4)
#    sigma_J = 0.02
    sigma_diag = 4*10**(-4)
    sigma_J = 0.03
#    sigma_diag = 2*10**(-8)
#    sigma_J = 0.03*10**(-8)
    

    outerAz = 0.*10**6
    innerAz = 0.*10**6
    innerRad = 0.*10**6
    injection_ind = 27
#    extraction_ind = 67
    extraction_ind = 5
#    extraction_ind = 45
    
    detuneSites = False
    detunedSites = [5,27,45,67]
    #siteDetuning = C1/C2*JJ
    siteDetuning = (C2-C1)*JJ/C2
#    siteDetuning = -200.*10**6
    
    defect = False
#    defect = True
    defect_loc = 91
    #defect_loc = TwoAzStart + 66
    defect_mag = -1500*10**6*3
    
    
    
    
    
#        ##option22, HW parameters, for hyperbolic 0701
#    ##under development
#    JJ = -64.0*10**6
#    omega0 = (8.025*10**9 - 2*76.*10**6) - 7*10**6 - 2*JJ
#    kappa = 0.2*10**6
##    sigma_diag = 6*10**(-4)
##    sigma_J = 0.02
##    sigma_diag = 3*10**(-4)
##    sigma_J = 0.02
#    sigma_diag = 2*10**(-8)
#    sigma_J = 0.03*10**(-8)
#    
##    outerAz = -34.*10**6
##    innerAz = 0.*10**6
##    innerRad = -20.*10**6
#    outerAz = 0.*10**6
#    innerAz = -10.*10**6
#    innerRad = 0.*10**6
#    injection_ind = 27
##    extraction_ind = 67
#    extraction_ind = 5
##    extraction_ind = 45
#    
#    detuneSites = False
#    detunedSites = [5,27,45,67]
#    #siteDetuning = C1/C2*JJ
#    siteDetuning = (C2-C1)*JJ/C2
##    siteDetuning = -200.*10**6
#    
#    defect = False
##    defect = True
#    defect_loc = 91
#    #defect_loc = TwoAzStart + 66
#    defect_mag = -1500*10**6*3
    

    FullDetuningControl = False




if ModeFamily == 'FW':
    
#    ##option20, variant on FW data  for hyperbolic 05_03
#    ##under development
#    JJ = -71.8*10**6*2
#    omega0 = (15.737*10**9- 2.5*10**6 - 2*JJ)
#    kappa = 0.7*10**6
##    kappa =1.4*10**6
##    sigma_diag = 3*10**(-4)
##    sigma_J = 0.02
#    sigma_diag = 2*10**(-8)
#    sigma_J = 0.03*10**(-8)
#    
#    outerAz = 0*10**6
#    innerAz = 0.*10**6
#    innerRad = 0*10**6
#    injection_ind = 27
#    #extraction_ind = 67
#    #extraction_ind = 5
#    extraction_ind = 45
#    
#    detuneSites = False
#    detunedSites = [5,27,45,67]
#    #siteDetuning = C1/C2*JJ
#    siteDetuning = (C2-C1)*JJ/C2
#    #siteDetuning = siteDetuning
#    
#    defect = False
##    defect = True
#    #defect_loc = 102
#    defect_loc = 91
#    #defect_loc = TwoAzStart + 66
#    defect_mag = -1500*10**6*3
    
    
    
#    ##option25, from HW parameters, for hyperbolic 0701
#    ##under development
#    JJ = -66.0*10**6*2
#    freqOffset = 0.*10**6;
#    freqOffset = -35.*10**6;
#    omega0 = ((8.025*10**9 - 2*76.*10**6) - 17*10**6 )*2- 2*JJ + freqOffset
#    kappa = 0.7*10**6*2
#    sigma_diag = 6*10**(-4)*2
#    sigma_J = 0.04*2
##    sigma_diag = 2*10**(-8)
##    sigma_J = 0.03*10**(-8)
#    
#
#    outerAz = 0.*10**6*2
#    innerAz = 10.*10**6*2
#    innerRad = 0.*10**6*2
#    injection_ind = 27
##    extraction_ind = 67
#    extraction_ind = 5
##    extraction_ind = 45
#    
#    detuneSites = False
#    detunedSites = [5,27,45,67]
#    #siteDetuning = C1/C2*JJ
#    siteDetuning = (C2-C1)*JJ/C2
##    siteDetuning = -200.*10**6
#    
#    defect = False
##    defect = True
#    defect_loc = 91
#    #defect_loc = TwoAzStart + 66
#    defect_mag = -1500*10**6*3
#    
    
    
    
#    ##option26, from FW parameters, for hyperbolic 0701
#    ##under development
##    JJ = -66.50*10**6*2
#    JJ = -67.1*10**6*2
#    freqOffset = 0.*10**6;
#    freqOffset = -14.*10**6;
#    omega0 = ((8.025*10**9 - 2*76.*10**6) - 17*10**6 )*2- 2*JJ + freqOffset
#    kappa = 0.7*10**6*2
#    sigma_diag = 4*10**(-4)*2
#    sigma_J = 0.03*2
##    sigma_diag = 2*10**(-8)
##    sigma_J = 0.03*10**(-8)
#    
#    FullDetuningControl = False
#    outerAz = 0.*10**6*2
#    innerAz = -5.*10**6*2
##    innerAz = 0.*10**6*2
#    innerRad = 0.*10**6*2
#    injection_ind = 27
##    extraction_ind = 67
#    extraction_ind = 5
##    extraction_ind = 45
#    
#    detuneSites = False
#    detunedSites = [5,27,45,67]
#    #siteDetuning = C1/C2*JJ
#    siteDetuning = (C2-C1)*JJ/C2
##    siteDetuning = -200.*10**6
#    
#    defect = False
##    defect = True
#    defect_loc = 91
#    #defect_loc = TwoAzStart + 66
#    defect_mag = -1500*10**6*3
    
    
    
    
#    ##option27, from FW parameters, for hyperbolic 0701, larger t
#    ##under development
##    JJ = -66.50*10**6*2
#    JJ = -69*10**6*2
#    freqOffset = 0.*10**6;
#    freqOffset = -16.*10**6;
#    omega0 = ((8.025*10**9 - 2*76.*10**6) - 17*10**6 )*2- 2*JJ + freqOffset
#    kappa = 0.7*10**6*2
#    sigma_diag = 6*10**(-4)*2
#    sigma_J = 0.04*2
##    sigma_diag = 2*10**(-8)
##    sigma_J = 0.03*10**(-8)
#    
#    FullDetuningControl = False
#    outerAz = 0.*10**6*2
#    innerAz = 5.*10**6*2
##    innerAz = 0.*10**6*2
#    innerRad = 0.*10**6*2
#    injection_ind = 27
##    extraction_ind = 67
#    extraction_ind = 5
##    extraction_ind = 45
#    
#    detuneSites = False
#    detunedSites = [5,27,45,67]
#    #siteDetuning = C1/C2*JJ
#    siteDetuning = (C2-C1)*JJ/C2
##    siteDetuning = -200.*10**6
#    
#    defect = False
##    defect = True
#    defect_loc = 91
#    #defect_loc = TwoAzStart + 66
#    defect_mag = -1500*10**6*3
    


    
    
#    ##option28, from FW parameters, for hyperbolic 0701, guessing
#    ##under development
##    JJ = -66.50*10**6*2
#    JJ = -68.4*10**6*2
#    freqOffset = 0.*10**6;
#    freqOffset = -20.*10**6;
#    omega0 = ((8.025*10**9 - 2*76.*10**6) - 17*10**6 )*2- 2*JJ + freqOffset
#    kappa = 0.7*10**6*2
#    sigma_diag = 3*10**(-4)*2
#    sigma_J = 0.02*2
##    sigma_diag = 2*10**(-4)*2
##    sigma_J = 0.01*2
##    sigma_diag = 2*10**(-8)
##    sigma_J = 0.03*10**(-8)
#    
#    FullDetuningControl = True
#    ZeroAzDet = 0.*10**6*2
#    OneRadDet = 0.*10**6*2
#    OneAzDet = -10.*10**6*2
#    TwoRadDet = 0.*10**6*2
#    TwoAzDet = 7.*10**6*2
##    outerAz = 0.*10**6*2
##    innerAz = -10.*10**6*2
###    innerAz = 0.*10**6*2
##    innerRad = 0.*10**6*2
#    injection_ind = 27
##    extraction_ind = 67
##    extraction_ind = 5
#    extraction_ind = 45
#    
#    detuneSites = False
#    detunedSites = [5,27,45,67]
#    #siteDetuning = C1/C2*JJ
#    siteDetuning = (C2-C1)*JJ/C2
##    siteDetuning = -200.*10**6
#    
#    defect = False
##    defect = True
#    defect_loc = 91
#    #defect_loc = TwoAzStart + 66
#    defect_mag = -1500*10**6*3
    
    
#    ##option29, from FW parameters, for hyperbolic 0701, experimentation with new BG techniques
#    ##under development
#    JJ = -68.4*10**6*2
#    freqOffset = -20.*10**6;
#    omega0 = ((8.025*10**9 - 2*76.*10**6) - 17*10**6 )*2- 2*JJ + freqOffset
#    kappa = 0.7*10**6*2
#    sigma_diag = 3*10**(-4)*2
#    sigma_J = 0.02*2
##    sigma_diag = 2*10**(-4)*2
##    sigma_J = 0.01*2
##    sigma_diag = 2*10**(-8)
##    sigma_J = 0.03*10**(-8)
#    
#    FullDetuningControl = True
##    ZeroAzDet = -5.*10**6*2
##    OneRadDet = -5.*10**6*2
##    OneAzDet = -7.*10**6*2
##    TwoRadDet = -5.*10**6*2
##    TwoAzDet = 13.*10**6*2
#    ZeroAzDet = -5.*10**6*2
#    OneRadDet = -5.*10**6*2
#    OneAzDet = -7.*10**6*2
#    TwoRadDet = -5.*10**6*2
#    TwoAzDet = 10.*10**6*2
#    injection_ind = 27
##    extraction_ind = 67
##    extraction_ind = 5
#    extraction_ind = 45
#    
#    detuneSites = False
#    detunedSites = [5,27,45,67]
#    #siteDetuning = C1/C2*JJ
#    siteDetuning = (C2-C1)*JJ/C2
##    siteDetuning = -200.*10**6
#    
#    defect = False
##    defect = True
#    defect_loc = 91
#    #defect_loc = TwoAzStart + 66
#    defect_mag = -1500*10**6*3
    
    
    
    
#    ##option30, from FW parameters, for hyperbolic 0701, experimentation with new BG techniques, small t
#    ##under development
#    JJ = -66.50*10**6*2
#    freqOffset = 0.*10**6;
#    omega0 = ((8.025*10**9 - 2*76.*10**6) - 17*10**6 )*2- 2*JJ + freqOffset
#    kappa = 0.7*10**6*2
#    sigma_diag = 3*10**(-4)*2
#    sigma_J = 0.02*2
##    sigma_diag = 2*10**(-4)*2
##    sigma_J = 0.01*2
##    sigma_diag = 2*10**(-8)
##    sigma_J = 0.03*10**(-8)
#    
#    FullDetuningControl = True
##    ZeroAzDet = -5.*10**6*2
##    OneRadDet = -5.*10**6*2
##    OneAzDet = -7.*10**6*2
##    TwoRadDet = -5.*10**6*2
##    TwoAzDet = 13.*10**6*2
#    ZeroAzDet = -5.*10**6*2
#    OneRadDet = -5.*10**6*2
#    OneAzDet = -7.*10**6*2
#    TwoRadDet = -5.*10**6*2
#    TwoAzDet = 10.*10**6*2
#    injection_ind = 27
##    extraction_ind = 67
##    extraction_ind = 5
#    extraction_ind = 45
#    
#    detuneSites = False
#    detunedSites = [5,27,45,67]
#    #siteDetuning = C1/C2*JJ
#    siteDetuning = (C2-C1)*JJ/C2
##    siteDetuning = -200.*10**6
#    
#    defect = False
##    defect = True
#    defect_loc = 91
#    #defect_loc = TwoAzStart + 66
#    defect_mag = -1500*10**6*3
#    

    
    
#    ##option31, from FW parameters, for hyperbolic 0701, experimentation with new BG techniques, wild card
#    ##under development
#    JJ = -67.*10**6*2
#    freqOffset = -5.*10**6;
#    omega0 = ((8.025*10**9 - 2*76.*10**6) - 17*10**6 )*2- 2*JJ + freqOffset
#    kappa = 0.7*10**6*2
#    sigma_diag = 3*10**(-4)*2
#    sigma_J = 0.02*2
##    sigma_diag = 2*10**(-4)*2
##    sigma_J = 0.01*2
##    sigma_diag = 2*10**(-8)
##    sigma_J = 0.03*10**(-8)
#    
#    FullDetuningControl = True
##    ZeroAzDet = -5.*10**6*2
##    OneRadDet = -5.*10**6*2
##    OneAzDet = -7.*10**6*2
##    TwoRadDet = -5.*10**6*2
##    TwoAzDet = 13.*10**6*2
#    ZeroAzDet = -5.*10**6*2
#    OneRadDet = -5.*10**6*2
#    OneAzDet = -7.*10**6*2
#    TwoRadDet = -5.*10**6*2
#    TwoAzDet = 7.*10**6*2
#    injection_ind = 27
##    extraction_ind = 67
##    extraction_ind = 5
#    extraction_ind = 45
#    
#    detuneSites = False
#    detunedSites = [5,27,45,67]
##    detunedSites = [119, 127, 135, 140, 143]
#    #siteDetuning = C1/C2*JJ
##    siteDetuning = (C2-C1)*JJ/C2
#    siteDetuning = -10.*10**6
#    
#    defect = False
##    defect = True
#    defect_loc = 127
#    #defect_loc = TwoAzStart + 66
#    defect_mag = +10.*10**6
    
    
    
#    ##option32, from FW parameters, for hyperbolic 0701, experimentation with new BG techniques, wildcard 2
#    ##under development
##    JJ = -68.2*10**6*2
##    freqOffset = -23.*10**6;
#    JJ = -68.9*10**6*2
#    freqOffset = -30.*10**6;
#    omega0 = ((8.025*10**9 - 2*76.*10**6) - 17*10**6 )*2- 2*JJ + freqOffset
#    kappa = 0.7*10**6*2
#    sigma_diag = 3*10**(-4)*2
#    sigma_J = 0.02*2
##    sigma_diag = 2*10**(-4)*2
##    sigma_J = 0.01*2
##    sigma_diag = 2*10**(-8)
##    sigma_J = 0.03*10**(-8)
#    
#    FullDetuningControl = True
#    ZeroAzDet = 5.*10**6*2
#    OneRadDet = 0.*10**6*2
#    OneAzDet = -7.*10**6*2
#    TwoRadDet = 5.*10**6*2
#    TwoAzDet = 13.*10**6*2
#    injection_ind = 27
##    extraction_ind = 67
##    extraction_ind = 5
#    extraction_ind = 45
#    
#    detuneSites = False
#    detunedSites = [5,27,45,67]
#    #siteDetuning = C1/C2*JJ
#    siteDetuning = (C2-C1)*JJ/C2
##    siteDetuning = -200.*10**6
#    
#    defect = False
##    defect = True
#    defect_loc = 91
#    #defect_loc = TwoAzStart + 66
#    defect_mag = -1500*10**6*3
#    
    
#    ##option33, from FW parameters, for hyperbolic 0701, experimentation with new BG techniques, wildcard 3
#    ##under development
##    JJ = -68.2*10**6*2
##    freqOffset = -23.*10**6;
#    JJ = -68.1*10**6*2
#    freqOffset = -26.*10**6;
#    omega0 = ((8.025*10**9 - 2*76.*10**6) - 17*10**6 )*2- 2*JJ + freqOffset
#    kappa = 0.7*10**6*2
#    sigma_diag = 3*10**(-4)*2
#    sigma_J = 0.02*2
##    sigma_diag = 2*10**(-4)*2
##    sigma_J = 0.01*2
##    sigma_diag = 2*10**(-8)
##    sigma_J = 0.03*10**(-8)
#    
#    FullDetuningControl = True
##    ZeroAzDet = 5.*10**6*2
##    OneRadDet = 0.*10**6*2
##    OneAzDet = 5.*10**6*2
##    TwoRadDet = 5.*10**6*2
##    TwoAzDet = 15.*10**6*2
#    #
##    ZeroAzDet = 5.*10**6*2
##    OneRadDet = 10.*10**6*2
##    OneAzDet = 3.*10**6*2
##    TwoRadDet = 5.*10**6*2
##    TwoAzDet = 15.*10**6*2
#    # 
##    ZeroAzDet = 5.*10**6*2
##    OneRadDet = 10.*10**6*2
##    OneAzDet = 5.*10**6*2
##    TwoRadDet = 5.*10**6*2
##    TwoAzDet = 15.*10**6*2
#    #
##    ZeroAzDet = 5.*10**6*2
##    OneRadDet = 0.*10**6*2
##    OneAzDet = -7.*10**6*2
##    TwoRadDet = 5.*10**6*2
##    TwoAzDet = 15.*10**6*2
#    #
#    ZeroAzDet = 5.*10**6*2
#    OneRadDet = -3.*10**6*2
#    OneAzDet = -7.*10**6*2
#    TwoRadDet = 5.*10**6*2
#    TwoAzDet = 15.*10**6*2
#    
#    injection_ind = 27
##    extraction_ind = 67
##    extraction_ind = 5
#    extraction_ind = 45
#    
#    detuneSites = False
#    detunedSites = [5,27,45,67]
#    #siteDetuning = C1/C2*JJ
#    siteDetuning = (C2-C1)*JJ/C2
##    siteDetuning = -200.*10**6
#    
#    defect = False
##    defect = True
#    defect_loc = 91
#    #defect_loc = TwoAzStart + 66
#    defect_mag = -1500*10**6*3
    

#    ##option33, from FW parameters, for hyperbolic 0701, experimentation with new BG techniques, wildcard 4
#    ##under development
##    JJ = -68.2*10**6*2
##    freqOffset = -23.*10**6;
#    JJ = -67.2*10**6*2
#    freqOffset = -15.*10**6;
##    JJ = -69.0*10**6*2
##    freqOffset = -35.*10**6;
#    omega0 = ((8.025*10**9 - 2*76.*10**6) - 17*10**6 )*2- 2*JJ + freqOffset
#    kappa = 0.7*10**6*2
#    sigma_diag = 3*10**(-4)*2
#    sigma_J = 0.02*2
##    sigma_diag = 2*10**(-4)*2
##    sigma_J = 0.01*2
##    sigma_diag = 2*10**(-8)
##    sigma_J = 0.03*10**(-8)
#    
#    FullDetuningControl = True
#    ZeroAzDet = 5.*10**6*2
#    OneRadDet = 0.*10**6*2
#    OneAzDet = 20.*10**6*2
#    TwoRadDet = 5.*10**6*2
#    TwoAzDet = 10.*10**6*2
#    
#    injection_ind = 27
##    extraction_ind = 67
##    extraction_ind = 5
#    extraction_ind = 45
#    
#    detuneSites = False
#    detunedSites = [5,27,45,67]
#    #siteDetuning = C1/C2*JJ
#    siteDetuning = (C2-C1)*JJ/C2
##    siteDetuning = -200.*10**6
#    
#    defect = False
##    defect = True
#    defect_loc = 91
#    #defect_loc = TwoAzStart + 66
#    defect_mag = -1500*10**6*3
    
    
    
    ##option35, from FW parameters, for hyperbolic 0701, experimentation with new BG techniques, wildcard 5
    ##under development
#    JJ = -68.2*10**6*2
#    freqOffset = -23.*10**6;
    JJ = -68.1*10**6*2
    freqOffset = -26.*10**6;
    omega0 = ((8.025*10**9 - 2*76.*10**6) - 17*10**6 )*2- 2*JJ + freqOffset
    kappa = 0.7*10**6*2
    sigma_diag = 3*10**(-4)*2
    sigma_J = 0.02*2
#    sigma_diag = 2*10**(-4)*2
#    sigma_J = 0.01*2
#    sigma_diag = 2*10**(-8)
#    sigma_J = 0.03*10**(-8)
    
    FullDetuningControl = True
    ZeroAzDet = 5.*10**6*2
    OneRadDet = -5.*10**6*2
    OneAzDet = -7.*10**6*2
    TwoRadDet = 5.*10**6*2 #8 seems very similar
    TwoAzDet = 15.*10**6*2
    # base
#    ZeroAzDet = 5.*10**6*2
#    OneRadDet = -3.*10**6*2
#    OneAzDet = -7.*10**6*2
#    TwoRadDet = 5.*10**6*2
#    TwoAzDet = 15.*10**6*2
    
    injection_ind = 27
#    extraction_ind = 67
#    extraction_ind = 5
    extraction_ind = 45
    
    detuneSites = False #true is also vaguely reasonable.
#    detunedSites = [5,27,45,67]
    detunedSites = [5,6,27,28,45,46,67,68]
#    detunedSites = [5,4,27,26,45,44,67,66]
    #siteDetuning = C1/C2*JJ
#    siteDetuning = (C2-C1)*JJ/C2
    siteDetuning = +10.*10**6
    
    defect = False
#    defect = True
    defect_loc = 91
    #defect_loc = TwoAzStart + 66
    defect_mag = -1500*10**6*3




##horsing around with the square lattice
#injection_ind = 0
#extraction_ind = 10
#injection_ind = 0
#extraction_ind = 8







#######
#load in a tiling, or make one
#######

t0 = time.time()

test = PlanarLayout(gon = 7, vertex = 3, side =1, radius_method = 'lin', modeType = ModeFamily)
#test = PlanarLayout(gon = 4, vertex = 4, side =1, radius_method = 'lin', modeType = ModeFamily)
test.populate(maxItter = 2)

t_elapsed = time.time()- t0

print 'generate time = ' + str(t_elapsed)



#####
#find all the classes of points
#####
ZeroAzStart = test.get_SDindex(0, 0, az = True)

OneAzStart = test.get_SDindex(0, 1, az = True)
OneRadStart = test.get_SDindex(0, 1, az = False)

TwoAzStart = test.get_SDindex(0, 2, az = True)
TwoRadStart = test.get_SDindex(0, 2, az = False)

boundaries = numpy.asarray([ZeroAzStart, OneAzStart, OneRadStart, TwoAzStart, TwoRadStart, test.SDx.size])



#####
#helper functions
#####
def calculate_in_out_states(couplingMode = 'mid', injection_ind = 9, extraction_ind = 33):
    if couplingMode == 'mid':
        injection_state = test.build_local_state_az(injection_ind)
    if couplingMode == 'end':
        numPoints = test.SDpoints[test.itter].shape[0]
        ind1 = numpy.mod(injection_ind, numPoints)
        ind2 = numpy.mod(injection_ind + 1, numPoints)
        injection_state = test.build_local_state_az(ind1) + test.build_local_state_az(ind2)
        injection_state = injection_state/numpy.linalg.norm(injection_state)
    
    if couplingMode == 'mid':
        extraction_state = test.build_local_state_az(extraction_ind)
    if couplingMode == 'end':
        numPoints = test.SDpoints[test.itter].shape[0]
        ind1 = extraction_ind
        ind2 = numpy.mod(extraction_ind + 1, numPoints)
        extraction_state = test.build_local_state_az(ind1) + test.build_local_state_az(ind2)
        extraction_state = extraction_state/numpy.linalg.norm(extraction_state)
        
    return [injection_state, extraction_state]
        
def simulate_overlaps(injection_state, extraction_state, eigMat ):
    '''gets powers, or at least absolute values '''
    injection_amps = numpy.dot(numpy.transpose(numpy.conj(eigMat)), injection_state)
    extraction_amps = numpy.dot(numpy.transpose(numpy.conj(eigMat)), extraction_state)
    injection_abs = numpy.abs(injection_amps)
    extraction_abs = numpy.abs(extraction_amps)

    trans_abs = injection_abs*extraction_abs
    return [trans_abs, injection_abs, extraction_abs]

def simulate_overlap_amps(injection_state, extraction_state, eigMat):

    #####
    injection_amps = numpy.dot(numpy.transpose(numpy.conj(eigMat)), injection_state)
    extraction_amps = numpy.dot(numpy.transpose(numpy.conj(eigMat)), extraction_state)

    trans_amps = injection_amps*numpy.conj(extraction_amps)


    return [trans_amps, injection_amps, extraction_amps]

def Lorr(xs, x0, Gamma, norm = 'L2'):
    '''Power Lorentzian'''
    numerator = Gamma/2./numpy.pi
    denominator = (xs-x0)**2 + (Gamma/2.)**2
    lorentz = numerator/denominator
    if norm == 'supnorm':
#        return lorentz/numpy.max(lorentz)
        return lorentz*numpy.pi*Gamma/2.
    else: #L2
        return lorentz
    
def Lorr_amp(xs, x0, Gamma, norm = 'L2'):
    '''amplitude Lorentzian '''
    numerator = numpy.sqrt(Gamma/2./numpy.pi)
    denominator = (xs-x0) + (Gamma/2.)*1j
    lorentz = numerator/denominator
    if norm == 'supnorm':
#        return lorentz/numpy.max(numpy.abs(lorentz))
        return lorentz* numpy.sqrt(numpy.pi*Gamma/2.)
    else: #L2
        return lorentz
    

def simulate_transmission(xs, sites, eigs, eigMat, kappa, couplingMode = 'mid', norm = 'L2'):
    '''gets powers, or at least abosulte values '''
    injection_ind = sites[0]
    extraction_ind = sites[1]
    [injection_state, extraction_state] = calculate_in_out_states(couplingMode, injection_ind, extraction_ind )
    
    [trans_abs, injection_abs, extraction_abs] = simulate_overlaps(injection_state, extraction_state, eigMat)
    
    simulated_transmission = numpy.zeros(len(xs))
    for ind in range(0,len(eigs)):
        tempLorr = Lorr(xs, eigs[ind], kappa, norm = norm)* trans_abs[ind]
        simulated_transmission = simulated_transmission + tempLorr
    
    return simulated_transmission

def simulate_transmission_amps(xs, sites, eigs, eigMat, kappa, couplingMode = 'mid', norm = 'L2'):
    injection_ind = sites[0]
    extraction_ind = sites[1]
    [injection_state, extraction_state] = calculate_in_out_states(couplingMode, injection_ind, extraction_ind )
    
    [trans_amps, injection_amps, extraction_amps] = simulate_overlap_amps(injection_state, extraction_state, eigMat)
    
#    trans_amps = numpy.abs(trans_amps)
    
    simulated_transmission_amp = numpy.zeros(len(xs))
    for ind in range(0,len(eigs)):
        tempLorr = Lorr_amp(xs, eigs[ind], kappa, norm = norm)* trans_amps[ind]
        simulated_transmission_amp = simulated_transmission_amp + tempLorr

    return simulated_transmission_amp

def simulate_device_transmission_amps(xs, couplings, sites, eigs, eigMat, kappa, couplingMode = 'mid', norm = 'L2'):
    injection_ind = sites[0]
    extraction_ind = sites[1]
    #there are four ways in and out for a given coupler
    [injection_state1, extraction_state1] = calculate_in_out_states(couplingMode = couplingMode, injection_ind = injection_ind, extraction_ind  = extraction_ind)
    [injection_state2, extraction_state2] = calculate_in_out_states(couplingMode = couplingMode, injection_ind = injection_ind+1, extraction_ind  = extraction_ind)
    [injection_state3, extraction_state3] = calculate_in_out_states(couplingMode = couplingMode, injection_ind = injection_ind, extraction_ind  = extraction_ind+1)
    [injection_state4, extraction_state4] = calculate_in_out_states(couplingMode = couplingMode, injection_ind = injection_ind+1, extraction_ind  = extraction_ind+1)
    
    [trans_coefs1, injection_amps1, extraction_amps1] = simulate_overlap_amps(injection_state1, extraction_state1, eigMat )
    [trans_coefs2, injection_amps2, extraction_amps2] = simulate_overlap_amps(injection_state2, extraction_state2, eigMat )
    [trans_coefs3, injection_amps3, extraction_amps3] = simulate_overlap_amps(injection_state3, extraction_state3, eigMat )
    [trans_coefs4, injection_amps4, extraction_amps4] = simulate_overlap_amps(injection_state4, extraction_state4, eigMat )
    
    simTransAmp1 = simulate_transmission_amps(xs, [injection_ind, extraction_ind], eigs, eigMat, kappa, couplingMode = couplingMode, norm = norm)
    simTransAmp2 = simulate_transmission_amps(xs, [injection_ind+1, extraction_ind], eigs, eigMat, kappa, couplingMode = couplingMode, norm = norm)
    simTransAmp3 = simulate_transmission_amps(xs, [injection_ind, extraction_ind+1], eigs, eigMat, kappa, couplingMode = couplingMode, norm = norm)
    simTransAmp4 = simulate_transmission_amps(xs, [injection_ind+1, extraction_ind+1], eigs, eigMat, kappa, couplingMode = couplingMode, norm = norm)

    transAmpOut = couplings[0]*simTransAmp1 + couplings[1]*simTransAmp2 + couplings[2]*simTransAmp3 + couplings[3]*simTransAmp4
    transCoefsOut = couplings[0]*trans_coefs1 + couplings[1]*trans_coefs2 + couplings[2]*trans_coefs3 + couplings[3]*trans_coefs4
    
#    if testMode:
#        return [trans_coefs1, trans_coefs2, trans_coefs3, trans_coefs4]

    return transAmpOut, transCoefsOut







def simulate_transmission_amps_withKappas(xs, sites, eigs, eigMat, kappa_vec, couplingMode = 'mid', norm = 'L2'):
    injection_ind = sites[0]
    extraction_ind = sites[1]
    [injection_state, extraction_state] = calculate_in_out_states(couplingMode, injection_ind, extraction_ind )
    
    [trans_amps, injection_amps, extraction_amps] = simulate_overlap_amps(injection_state, extraction_state, eigMat)
    
#    trans_amps = numpy.abs(trans_amps)
    
    simulated_transmission_amp = numpy.zeros(len(xs))
    for ind in range(0,len(eigs)):
        tempLorr = Lorr_amp(xs, eigs[ind], kappa_vec[ind], norm = norm)* trans_amps[ind]
        simulated_transmission_amp = simulated_transmission_amp + tempLorr

    return simulated_transmission_amp



def simulate_device_transmission_amps_variableKappa(xs, couplings, sites, eigs, eigMat, J, couplingMode = 'mid', norm = 'L2', returnKappas = False, kappa = 0.7*10**6):
    injection_ind = sites[0]
    extraction_ind = sites[1]
    #there are four ways in and out for a given coupler
    [injection_state1, extraction_state1] = calculate_in_out_states(couplingMode = couplingMode, injection_ind = injection_ind, extraction_ind  = extraction_ind)
    [injection_state2, extraction_state2] = calculate_in_out_states(couplingMode = couplingMode, injection_ind = injection_ind+1, extraction_ind  = extraction_ind)
    [injection_state3, extraction_state3] = calculate_in_out_states(couplingMode = couplingMode, injection_ind = injection_ind, extraction_ind  = extraction_ind+1)
    [injection_state4, extraction_state4] = calculate_in_out_states(couplingMode = couplingMode, injection_ind = injection_ind+1, extraction_ind  = extraction_ind+1)
    
    [trans_coefs1, injection_amps1, extraction_amps1] = simulate_overlap_amps(injection_state1, extraction_state1, eigMat )
    [trans_coefs2, injection_amps2, extraction_amps2] = simulate_overlap_amps(injection_state2, extraction_state2, eigMat )
    [trans_coefs3, injection_amps3, extraction_amps3] = simulate_overlap_amps(injection_state3, extraction_state3, eigMat )
    [trans_coefs4, injection_amps4, extraction_amps4] = simulate_overlap_amps(injection_state4, extraction_state4, eigMat )
    
    totalInCoupling = injection_amps1 + injection_amps2 + injection_amps3+ injection_amps4
    totalOutCoupling = extraction_amps1 + extraction_amps2 + extraction_amps3+ extraction_amps4
    kappaIn = numpy.abs(J) * totalInCoupling*numpy.conj(totalInCoupling)
    kappaOut = numpy.abs(J) * totalOutCoupling*numpy.conj(totalOutCoupling)
#    kappa_vec = 0.7*10**6 + (kappaIn + kappaOut)*2
#    kappa_vec = 0.7*10**6 + (kappaIn + kappaOut)*2*0.01 #gonna incorrectly assume that the other two ports have the same coefficients
    kappa_vec = kappa*1.0 + (kappaIn + kappaOut)*8*0.01 #was in use up to 5 pm 12/18/18
#    kappa_vec = kappa*1.0 + (kappaIn + kappaOut)*4*0.01

#    kappaIn = numpy.abs(J) * numpy.abs(totalInCoupling)
#    kappaOut = numpy.abs(J) * numpy.abs(totalOutCoupling )  
#    kappa_vec = 0.7*10**6 + (kappaIn + kappaOut)*2*0.003
#    kappa_vec = 1.4*10**6 * numpy.ones(len(eigs))# placeholder
    
    simTransAmp1 = simulate_transmission_amps_withKappas(xs, [injection_ind, extraction_ind], eigs, eigMat, kappa_vec, couplingMode = couplingMode, norm = norm)
    simTransAmp2 = simulate_transmission_amps_withKappas(xs, [injection_ind+1, extraction_ind], eigs, eigMat, kappa_vec, couplingMode = couplingMode, norm = norm)
    simTransAmp3 = simulate_transmission_amps_withKappas(xs, [injection_ind, extraction_ind+1], eigs, eigMat, kappa_vec, couplingMode = couplingMode, norm = norm)
    simTransAmp4 = simulate_transmission_amps_withKappas(xs, [injection_ind+1, extraction_ind+1], eigs, eigMat, kappa_vec, couplingMode = couplingMode, norm = norm)

    transAmpOut = couplings[0]*simTransAmp1 + couplings[1]*simTransAmp2 + couplings[2]*simTransAmp3 + couplings[3]*simTransAmp4
    transCoefsOut = couplings[0]*trans_coefs1 + couplings[1]*trans_coefs2 + couplings[2]*trans_coefs3 + couplings[3]*trans_coefs4
    
#    if testMode:
#        return [trans_coefs1, trans_coefs2, trans_coefs3, trans_coefs4]
    
#    kappa_vec = totalInCoupling
    if returnKappas:
        return transAmpOut, transCoefsOut, kappa_vec
    else:
        return transAmpOut, transCoefsOut



####
#set up frequency sweep
####
freq_range = 4
freq_res = 0.05

freqs = scipy.arange(-freq_range, freq_range, freq_res) + freq_res/2.
freq_bins = scipy.arange(-freq_range, freq_range+freq_res, freq_res)

[fullDOS, bins_out] = numpy.histogram(test.Es, freq_bins)


#####
#et up eigenenergies
#####
if FullDetuningControl:
    systematicOffsets = numpy.asarray([ZeroAzDet, OneAzDet,OneRadDet, TwoAzDet, TwoRadDet])
else:
    systematicOffsets = numpy.asarray([innerAz, innerAz,innerRad, outerAz, innerRad])
#systematicOffsets = numpy.asarray([5, 5,10, 0, 0])*10**6



realisticFreqs = scipy.arange(-6*JJ, 4*JJ, JJ/1000) + omega0


for dind in range(0,numSamples):
    
    
    matSize = test.H.shape[0]
    Hmat_all = omega0*numpy.random.normal(0,sigma_diag,(matSize,matSize))*numpy.identity(matSize)\
                 - JJ*numpy.random.normal(0,sigma_J,(matSize,matSize))*test.H \
                 - JJ*test.H\
                 + omega0*numpy.identity(matSize)
    
    #systematically ofset different resonator shapes
    for sind in range(0, len(boundaries)-1):
        startInd = boundaries[sind]
        stopInd = boundaries[sind+1]
        for site in range(startInd, stopInd):
            Hmat_all[site, site] = Hmat_all[site, site] + systematicOffsets[sind]
    
    
    
    #selectively detune some sites
    if detuneSites:
        print 'detuning sites'
        detunedIndices = numpy.zeros(len(detunedSites))
        for ind in range(0,len(detunedSites)):
            detunedIndices[ind] = int(test.get_SDindex(detunedSites[ind], 2, az = True))
        for ind in range(0, len(detunedSites)):
            site = int(detunedIndices[ind])
            Hmat_all[site,site] = Hmat_all[site,site] + siteDetuning ####!!!!!!mess up one lattice site
            
    
    if defect:
        Hmat_all[defect_loc, defect_loc] = Hmat_all[defect_loc, defect_loc] + defect_mag
    
    
    #fix the Hailtonian and make it actually Hermitian. Fixes errors from the random noise
    for rowdx in range(0, Hmat_all.shape[0]):
        for coldx in range(rowdx+1, Hmat_all.shape[0]):
            #print str(rowdx) + ' , ' +str(coldx)
            Hmat_all[coldx, rowdx] = Hmat_all[rowdx, coldx]
      
     
    
    [Es_all, Psis_all] = scipy.linalg.eigh(Hmat_all)
    Eorder_all = numpy.argsort(Es_all)
    
    Es_perfect = omega0 - test.Es*JJ
    
    #couplin ratios for different setups
    #proper couplings at the input and output ports of the lattice  is weak +, strong-, weak + , strong -
    if ModeFamily == 'FW':
    #    symCouplings = numpy.sqrt(2)*numpy.asarray([1,1,1,1])
    #    asymCouplings = numpy.asarray([1,0,0,0])
    #    genCouplings = numpy.asarray([1.,numpy.sqrt(midFraction),1*numpy.sqrt(midFraction),midFraction])/numpy.sqrt(1 + midFraction)
        symCouplings = numpy.sqrt(2)*numpy.asarray([1,1,1,1])
        asymCouplings = numpy.asarray([0,0,0,1])    
        genCouplings = numpy.asarray([midFraction,numpy.sqrt(midFraction),numpy.sqrt(midFraction),1])/numpy.sqrt(1 + midFraction)
    elif ModeFamily == 'HW':
    #    symCouplings = numpy.sqrt(2)*numpy.asarray([1,-1,-1,1])
    #    asymCouplings = numpy.asarray([1,0,0,0])
    #    genCouplings = numpy.asarray([1.,-numpy.sqrt(midFraction),-1*numpy.sqrt(midFraction),midFraction])/numpy.sqrt(1 + midFraction)
        symCouplings = numpy.sqrt(2)*numpy.asarray([1,-1,-1,1])
        asymCouplings = numpy.asarray([0,0,0,1])    
        genCouplings = numpy.asarray([midFraction,-numpy.sqrt(midFraction),-numpy.sqrt(midFraction),1])/numpy.sqrt(1 + midFraction)
    else:
        raise ValueError, 'needs to be FW or HW ModeFamily'
        
    
    
    
    
    
    
    
    
    
    
    #####
    #general transmission figure.
    #####
    
    [injection_state1, extraction_state1] = calculate_in_out_states(couplingMode = 'mid', injection_ind = injection_ind, extraction_ind  = extraction_ind)
    [injection_state2, extraction_state2] = calculate_in_out_states(couplingMode = 'mid', injection_ind = injection_ind+1, extraction_ind  = extraction_ind+1)
    injection_state = (injection_state1+injection_state2)/numpy.sqrt(2)
    extraction_state = (extraction_state1+extraction_state2)/numpy.sqrt(2)
    
    
    #disorder free transmission
    #[trans_coefs, injection_amp, extraction_amp] = simulate_overlap_amps(injection_state, extraction_state, test.Psis )
    
    #symmetricTrans_amp, symmetricTrans_coefs = simulate_device_transmission_amps(realisticFreqs, symCouplings, [injection_ind, extraction_ind], Es_perfect, test.Psis, kappa)
    #asymmetricTrans_amp, asymmetricTrans_coefs = simulate_device_transmission_amps(realisticFreqs, asymCouplings, [injection_ind, extraction_ind], Es_perfect, test.Psis, kappa)
    #generalTrans_amp, generalTrans_coefs = simulate_device_transmission_amps(realisticFreqs, genCouplings, [injection_ind, extraction_ind], Es_perfect, test.Psis, kappa)
    symmetricTrans_amp, symmetricTrans_coefs = simulate_device_transmission_amps_variableKappa(realisticFreqs, symCouplings, [injection_ind, extraction_ind], Es_perfect, test.Psis, JJ)
    asymmetricTrans_amp, asymmetricTrans_coefs = simulate_device_transmission_amps_variableKappa(realisticFreqs, asymCouplings, [injection_ind, extraction_ind], Es_perfect, test.Psis, JJ)
    generalTrans_amp, generalTrans_coefs = simulate_device_transmission_amps_variableKappa(realisticFreqs, genCouplings, [injection_ind, extraction_ind], Es_perfect, test.Psis, JJ)
    
    
    
    #convert to transmmitted Power
    symmetricTrans = symmetricTrans_amp* numpy.conj(symmetricTrans_amp)
    asymmetricTrans = asymmetricTrans_amp* numpy.conj(asymmetricTrans_amp)
    generalTrans = generalTrans_amp* numpy.conj(generalTrans_amp)
    
    
    
    
    
    
    
    
    ######
    #compare the different transmissions
    ######
    DOStrans_amp = numpy.zeros(len(realisticFreqs))
    for ind in range(0,len(Es_perfect)):
        tempLorr = Lorr_amp(realisticFreqs, Es_perfect[ind], 10**6, norm = 'L2')* 1
        DOStrans_amp = DOStrans_amp + tempLorr
    DOStrans = DOStrans_amp*numpy.conj(DOStrans_amp)
    
    temp1 = simulate_transmission_amps(realisticFreqs, [injection_ind, extraction_ind], Es_perfect, test.Psis, kappa, couplingMode = 'mid', norm = 'L2')
    temp2 = temp1*numpy.conj(temp1)
    
    temp3 = simulate_transmission(realisticFreqs, [injection_ind, extraction_ind], Es_perfect, test.Psis, kappa, couplingMode = 'mid', norm = 'L2')
    temp3 = temp3/100
    
    
    
    
    
    
    
    
    
    
    ####look at the effect of smaple disorder
    #symmetricTrans_amp_disorder, symmetricTrans_coefs_disorder = simulate_device_transmission_amps(realisticFreqs, symCouplings, [injection_ind, extraction_ind], Es_all, Psis_all, kappa)
    #asymmetricTrans_amp_disorder, asymmetricTrans_coefs_disorder = simulate_device_transmission_amps(realisticFreqs, asymCouplings, [injection_ind, extraction_ind], Es_all, Psis_all, kappa)
    #generalTrans_amp_disorder, generalTrans_coefs_disorder = simulate_device_transmission_amps(realisticFreqs, genCouplings, [injection_ind, extraction_ind], Es_all, Psis_all, kappa)
    symmetricTrans_amp_disorder, symmetricTrans_coefs_disorder = simulate_device_transmission_amps_variableKappa(realisticFreqs, symCouplings, [injection_ind, extraction_ind], Es_all, Psis_all, JJ, kappa = kappa)
    asymmetricTrans_amp_disorder, asymmetricTrans_coefs_disorder = simulate_device_transmission_amps_variableKappa(realisticFreqs, asymCouplings, [injection_ind, extraction_ind], Es_all, Psis_all, JJ, kappa = kappa)
    generalTrans_amp_disorder, generalTrans_coefs_disorder = simulate_device_transmission_amps_variableKappa(realisticFreqs, genCouplings, [injection_ind, extraction_ind], Es_all, Psis_all, JJ, kappa = kappa)
    
    #make a stransmission simulation that can use empirical background
    generalTrans_amp_disorder_forFiltering, generalTrans_coefs_disorder_forFiltering = simulate_device_transmission_amps_variableKappa(transmissionFreqs1, genCouplings, [injection_ind, extraction_ind], Es_all, Psis_all, JJ, kappa = kappa)
    
    #generalTrans_amp_disorder, generalTrans_coefs_disorder = simulate_device_transmission_amps_variableKappa(realisticFreqs, genCouplings, [injection_ind, extraction_ind], Es_all, Psis_all, JJ, norm = 'supnorm')
    tempTrans_amp_disorder, tempTrans_coefs_disorderm, tempKappas = simulate_device_transmission_amps_variableKappa(realisticFreqs, genCouplings, [injection_ind, extraction_ind], Es_all, Psis_all, JJ, returnKappas = True, kappa = kappa)
    
    
    #convert to transmitted power
    symmetricTrans_disorder = symmetricTrans_amp_disorder* numpy.conj(symmetricTrans_amp_disorder)
    asymmetricTrans_disorder = asymmetricTrans_amp_disorder* numpy.conj(asymmetricTrans_amp_disorder)
    generalTrans_disorder = generalTrans_amp_disorder* numpy.conj(generalTrans_amp_disorder)
    
    
    
    
    
    ######
    #compare the different transmissions, with disorder
    ######
    DOStrans_amp_disorder = numpy.zeros(len(realisticFreqs))
    for ind in range(0,len(Es_perfect)):
        tempLorr = Lorr_amp(realisticFreqs, Es_all[ind], kappa, norm = 'L2')* 1
        DOStrans_amp_disorder = DOStrans_amp_disorder + tempLorr
    DOStrans_disorder = DOStrans_amp_disorder*numpy.conj(DOStrans_amp_disorder)
    
    temp4 = simulate_transmission_amps(realisticFreqs, [injection_ind, extraction_ind], Es_all, Psis_all, 10**6, couplingMode = 'mid', norm = 'L2')
    temp5 = temp4*numpy.conj(temp4)
    
    temp6 = simulate_transmission(realisticFreqs, [injection_ind, extraction_ind], Es_all, Psis_all, 10**6, couplingMode = 'mid', norm = 'L2')
    temp6 = temp6/100
    
    
    
    
    
    
    
    
    ######
    #plot data and simulation together
    ######
    amp1 = 2*10**4;
    amp2 = 1*10**4;
    
    #amp1 = 2*10**7;
    #amp2 = 1*10**7;
    
    
    
    #######
    #plot in log scale
    ######
    if ModeFamily == 'HW':
        tempamp1 = 1.*10**-0
        #tempamp2 = 1.*10**3
        tempamp2 = 3.*10**4
    else:
        tempamp1 = 1.*10**-0
        #tempamp2 = 1.*10**3
        tempamp2 = 0.5*10**3
    
#    BGpower = 10.**(-10)
    BGpower = 2*10.**(-10)
    
#    #general simulation with incoherent background.
    generalTrans_disorder_dB_incoherentBG = numpy.log10(generalTrans_disorder*tempamp2 + BGpower)*10 
#    generalTrans_disorder_dB = generalTrans_disorder_dB_incoherentBG
    
    #general simulation adding coherent background to get fano shapes, hoepfully
    
    BGamp = numpy.sqrt(BGpower)
    ampCoefficient = numpy.sqrt(tempamp2)
    simulatedField = generalTrans_amp_disorder*ampCoefficient + BGamp
    simulatedPower = simulatedField*numpy.conj(simulatedField)
    generalTrans_disorder_dB_coherentBG = numpy.log10(simulatedPower)*10
    
    
    filterWidth = 30
#    filterWidth = 60
    lowPassedData1 = scipy.ndimage.filters.gaussian_filter(transmissionVec1_dB, filterWidth)
    filteredData1 = transmissionVec1_dB - lowPassedData1 + numpy.mean(lowPassedData1)
    
    lowPassedData1_lin = 10**(lowPassedData1/10)
    lowPassedData1_lin_amp = numpy.sqrt(lowPassedData1_lin)
    
    
    simulatedField = generalTrans_amp_disorder_forFiltering*ampCoefficient + lowPassedData1_lin_amp
    simulatedPower = simulatedField*numpy.conj(simulatedField)
    generalTrans_disorder_dB_filterBG = numpy.log10(simulatedPower)*10
    
    
    #choose your option
#    generalTrans_disorder_dB = generalTrans_disorder_dB_incoherentBG
#    plotFreqs = realisticFreqs
#    generalTrans_disorder_dB = generalTrans_disorder_dB_coherentBG
#    plotFreqs = realisticFreqs
    generalTrans_disorder_dB = generalTrans_disorder_dB_filterBG
    plotFreqs = transmissionFreqs1
    
    
    
    
    
#    tempamp3 = 1.*10**-1
#    DOStrans_disorder_dB2 = numpy.log10(DOStrans_disorder*tempamp3 + 10**-10)*10
    #DOStrans_disorder_dB2 = numpy.log10(DOStrans_disorder*tempamp3 + 10**-8.5)*10
    
    
    if dind == 0:
        #plot single disorder realization
        pylab.figure(9)
        pylab.clf()
        #ax= pylab.subplot(1,2,1)
        ax= pylab.subplot(1,1,1)
        ###old figure
        pylab.plot(transmissionFreqs1/10**9, transmissionVec1_dB, 'mediumblue', linewidth=1, label = 'expteriment')
        pylab.plot(realisticFreqs/10**9,generalTrans_disorder_dB_incoherentBG, 'firebrick', linewidth=0.5, label = 'old theory')
        pylab.plot(transmissionFreqs1/10**9,lowPassedData1, 'deepskyblue', linewidth=1, label = 'low passed data')
        pylab.plot(transmissionFreqs1/10**9,generalTrans_disorder_dB_filterBG, 'green', linewidth=1, label = 'theory with empirical background')
        
        
        
        pylab.xlabel('Frequency (GHz)')
        pylab.ylabel('Transmission (dB)')
        pylab.title('Experiment v Theory')
        ax.legend()
        
        titleStr = 'Device Transmission      input and output sites: ' + str(injection_ind) + ' , ' + str(extraction_ind) + '     midFraction = ' + str(midFraction)
        #titleStr = 'Device Transmission compared to theory : '
        pylab.suptitle(titleStr)
        
        pylab.show()
        
        
        
        pylab.figure(10)
        pylab.clf()
        #ax= pylab.subplot(1,2,1)
        ax= pylab.subplot(1,1,1)
        ###old figure
#        pylab.plot(transmissionFreqs1/10**9, transmissionVec1_dB, 'mediumblue', linewidth=1, label = 'expteriment')
        pylab.plot(realisticFreqs/10**9,generalTrans_disorder_dB_incoherentBG, 'mediumblue', linewidth=0.5, label = 'theory with incoherent coherent BG')
        pylab.plot(realisticFreqs/10**9,generalTrans_disorder_dB_coherentBG, 'firebrick', linewidth=0.5, label = 'theory with flat coherent BG')
#        pylab.plot(transmissionFreqs1/10**9,lowPassedData1, 'deepskyblue', linewidth=1, label = 'low passed data')
        pylab.plot(transmissionFreqs1/10**9,generalTrans_disorder_dB_filterBG, 'deepskyblue', linewidth=1, label = 'theory with empirical coherent BG')
        
        pylab.xlabel('Frequency (GHz)')
        pylab.ylabel('Transmission (dB)')
        pylab.title('Different theories')
        ax.legend()
        
        titleStr = 'Device Transmission      input and output sites: ' + str(injection_ind) + ' , ' + str(extraction_ind) + '     midFraction = ' + str(midFraction)
        #titleStr = 'Device Transmission compared to theory : '
        pylab.suptitle(titleStr)
        
        pylab.show()
        
        
        
        pylab.figure(11)
        pylab.clf()
        #ax= pylab.subplot(1,2,1)
        ax= pylab.subplot(1,1,1)
        ###old figure
        pylab.plot(transmissionFreqs1/10**9,numpy.real(generalTrans_amp_disorder_forFiltering), 'mediumblue', linewidth=0.5, label = 'real part')
        pylab.plot(transmissionFreqs1/10**9,numpy.imag(generalTrans_amp_disorder_forFiltering), 'firebrick', linewidth=0.5, label = 'imaginary part')
        
        pylab.xlabel('Frequency (GHz)')
        pylab.ylabel('Transmission (linear ?)')
        pylab.title('Simulated transmission')
        ax.legend()
        
        titleStr = 'Device Transmission      input and output sites: ' + str(injection_ind) + ' , ' + str(extraction_ind) + '     midFraction = ' + str(midFraction)
        #titleStr = 'Device Transmission compared to theory : '
        pylab.suptitle(titleStr)
        
        pylab.show()
        
        
        
        #plot single disorder realization
        pylab.figure(12)
        pylab.clf()
        #ax= pylab.subplot(1,2,1)
        ax= pylab.subplot(1,1,1)
        ###old figure
        pylab.plot(transmissionFreqs1/10**9, transmissionVec1_dB, 'mediumblue', linewidth=1, label = 'expteriment')
        pylab.plot(transmissionFreqs1/10**9,generalTrans_disorder_dB_filterBG, 'firebrick', linewidth=0.75, label = 'theory with empirical background')
        
        
        
        pylab.xlabel('Frequency (GHz)')
        pylab.ylabel('Transmission (dB)')
        pylab.title('Experiment v Theory')
        ax.legend()
        
        titleStr = 'Device Transmission      input and output sites: ' + str(injection_ind) + ' , ' + str(extraction_ind) + '     midFraction = ' + str(midFraction)
        #titleStr = 'Device Transmission compared to theory : '
        pylab.suptitle(titleStr)
        
        pylab.show()
        
        
        
        
        
        #plot together, plot data and a few things the first time
        ##############
        fig = pylab.figure(1)
        pylab.clf()
        ax= pylab.subplot(1,1,1)
        ###old figure
        pylab.plot(transmissionFreqs1/10**9, transmissionVec1_dB, 'b', linewidth=1, label = 'experiment', alpha = 1)
        pylab.xlabel('Frequency (GHz)')
        pylab.ylabel('Transmission (dB)')
        pylab.title('Experiment v Theory')
#        ax.legend()
        
        titleStr = 'Device Transmission      input and output sites: ' + str(injection_ind) + ' , ' + str(extraction_ind) + '     midFraction = ' + str(midFraction)
        #titleStr = 'Device Transmission compared to theory : '
        pylab.suptitle(titleStr)
     
    pylab.figure(1)
    ax= pylab.subplot(1,1,1)
    #add the subsequent disorders  
    if numSamples == 1:
        pylab.plot(plotFreqs/10**9,generalTrans_disorder_dB, 'firebrick', linewidth=0.5, alpha = 1, label = 'theory')
    else:
        if dind == (numSamples-1):
            pylab.plot(plotFreqs/10**9,generalTrans_disorder_dB, 'firebrick', linewidth=0.5, alpha = 5./numSamples, label = 'theory')
        else:
            pylab.plot(plotFreqs/10**9,generalTrans_disorder_dB, 'firebrick', linewidth=0.5, alpha = 5./numSamples)
    pylab.show()  
    
    print dind
 
pylab.figure(1)
ax= pylab.subplot(1,1,1)   
ax.legend()
if ModeFamily == 'FW':   
    ax.set_xlim([15.5,16.5]) 
if ModeFamily == 'HW':   
    ax.set_xlim([7.8,8.3]) 
        
fig.set_size_inches([10.9, 5.4])  
#fig.savefig('disorderSimTrans.png', transparent= False, dpi = 200)
        

#alternate, wider version 
ax.set_xlim([15.68,16.5])    
ax.set_ylim([-120, -61]) 
fig.set_size_inches([8.76, 3.03]) 
#fig.set_size_inches([13.85, 3.98]) 
#fig.set_size_inches([10.35, 3.34])

pylab.tight_layout()
#fig.set_size_inches([13.85, 3.98]) 
#fig.set_size_inches([10.35, 3.34])
#fig.savefig('disorderSimTrans5.png', transparent= False, dpi = 200)     















arr

















