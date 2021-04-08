# -*- coding: utf-8 -*-
"""
Created on Fri Dec 18 13:51:28 2020

@author: Kollarlab
"""

dat = vna.get_channel_data(1)
plt.figure(14)
plt.clf()
plt.plot(dat['xaxis']/1e9, dat['mag'])
plt.xlabel('Frequency (GHz)')
plt.show()

tips = tip(multiple = True)


#plt.figure(14)
#plt.clf()
#plt.plot(mags[0])
#plt.plot(trans_mags[0])
#plt.show()