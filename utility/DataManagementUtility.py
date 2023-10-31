#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 16 15:39:15 2023

@author: kollar2
"""

import re
import random
# from scipy import *
import scipy
# import pylab
import matplotlib.pyplot as plt
import numpy as np
import time



import pickle
import datetime
import os
import sys

import copy



class DataManagementUtility(object):
    def __init__(self, master_data_folder):
        '''
        object to hold all the annoying path twidling 
        
        pass it the folder where wehave set up our data management library
        
        '''
        
        #Data storage structure
        self.communal_plot_repo_path =  master_data_folder
        
        #place for compiled plots
        self.plot_subfolder = os.path.join(self.communal_plot_repo_path,r'Compiled Plots')
        self.check_and_create(self.plot_subfolder)
        
        #raw data location
        self.raw_data_subfolder = os.path.join(self.communal_plot_repo_path,r'Raw Data Files') #intended file location, doesn't currently match Alicia's computer
        self.check_and_create(self.raw_data_subfolder)
        
        #location to save the curated data
        self.replotter_subfolder = os.path.join(self.communal_plot_repo_path,r'Curated Data')
        self.check_and_create(self.replotter_subfolder)
        
        return
    
    def check_and_create(self, path):
        '''See if a directory path exists. If not, create it. '''
        if os.path.exists(path):
            pass
        else:
            os.mkdir(path)
        return
    
    def autosave_directory(self, script_name, base_path = ''):
        '''create a directory automatically named after a script
        
        defaults to putting it in the plot subfolder
        '''
        
        if base_path == '':
            base_path = self.plot_subfolder 
    
        if script_name[-3:] == '.py':
            #I'm a python file
            script_name = script_name[:-3]
        elif script_name[-6:] == '.ipynb':
            #I'm an ipython notebook
            script_name = script_name[:-6]
        
        
        save_path = os.path.join(base_path, script_name)
        
        if os.path.exists(save_path):
            pass
        else:
            os.mkdir(save_path)
        return save_path

    def make_data_collection(self, collection_name):
        subgroupfolder = os.path.join(self.replotter_subfolder, collection_name) #folder for a group of related things
 
        self.check_and_create(subgroupfolder)
        
        return(subgroupfolder)

    def path_for_data_collection(self, collection_name):
        subgroupfolder = os.path.join(self.replotter_subfolder, collection_name) #folder for a group of related things
        
        if os.path.exists(subgroupfolder):
            return subgroupfolder
        else:
            raise ValueError('Data Collection with this name does not seem to exist.')
        
        