# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 13:46:45 2020

@author: Kollarlab
"""
import matplotlib.pyplot as plt
import numpy as np

import multiprocessing
multiprocessing.freeze_support() # <- may be required on windows

def plot(datax, datay, name):
    x = datax
    y = datay
    plt.scatter(x, y, label=name)
    plt.legend()
    plt.show()

def run_plot(name):
    
def multiP():
    for i in range(2):
        print('fds')
        p = multiprocessing.Process(target=plot, args=(i, i, i))
        p.start()

if __name__ == "__main__": 
    input('Value: ') 
    multiP()
    input('Value: ')