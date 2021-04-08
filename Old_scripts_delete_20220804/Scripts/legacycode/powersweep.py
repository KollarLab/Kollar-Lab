# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 16:11:08 2020

@author: Kollarlab
"""

import numpy
import matplotlib.pylab as plt

powers = numpy.linspace(-25,-10,2)

settings = vna.trans_default_settings()
settings['averages'] = 5
settings['frequency'] = 7.173e9
settings['span'] = 25e6 
settings['sweep_points'] = 901
settings['ifBW'] = 100
settings['channel'] = 1

mags, phases, axes = vna.power_sweep(settings, powers)
    

fig = plt.figure(dpi=150)
XX,YY = np.meshgrid(f1[0],powers)
im = plt.pcolormesh(XX/1e9,YY,m1,cmap="viridis")
plt.xlabel("Frequency (GHz)")
plt.ylabel(r"Power (dBm)")
ax = plt.gca()
ax.xaxis.set_major_locator(plt.MaxNLocator(5))
ax.yaxis.set_major_locator(plt.MaxNLocator(5))
axins = inset_axes(ax,
                    width="20%",  # width = 50% of parent_bbox width
                    height="2%",  # height : 5%
                    loc='upper right',
                    bbox_to_anchor=(0.0, 0., 1.0, 1.05),
                    bbox_transform=ax.transAxes,
                    borderpad=0,)
cbar = fig.colorbar(im,cax=axins,orientation="horizontal",ticks=plt.MaxNLocator(2))
cbar.set_label(r"$S_{21}$ (dB)",labelpad = -10,y=1,x=-0.35)
axins.xaxis.set_ticks_position("top")
axins.tick_params(labelsize=9)