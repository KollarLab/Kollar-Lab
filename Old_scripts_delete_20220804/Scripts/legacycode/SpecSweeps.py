# -*- coding: utf-8 -*-
"""
Created on Sun Nov  8 18:31:45 2020

@author: Kollarlab
"""
import userfuncs

saveDir = r'Z:\Data\HouckTaTransmon\20201109'

stamp = userfuncs.timestamp()

name = 'RFscanBroad'

scanname = name + '_' + stamp

settings = {}
settings['channel'] = 1
settings['averages'] = 80
settings['measurement'] = 'S21'
settings['meas_form'] = 'MLOG'
settings['start'] = 4.9e9
settings['stop'] = 5.2e9
settings['sweep_points'] = 1001
settings['RFpower'] = -5
settings['RFport'] = 3
settings['Mport'] = 2
settings['CAVport'] = 1
settings['CAVpower'] = -15
settings['CAVfreq'] = 7.173378e9
settings['ifBW'] = 500

powers = numpy.linspace(-15,-5,25)

HWattenuation = -30
initsettings = settings

mags = []
phases = []
for power in powers:
    print('Power: {}, final power: {}'.format(power, powers[-1]))
    initsettings['RFpower'] = power
    vna.spec_meas(settings)
    m = vna.get_trace(1,'spec')
    p = vna.get_trace(1,'phase')
    mags.append(m)
    phases.append(p)
    
freqs = vna.get_channel_axis(1)

fig2 = plt.figure(8)
plt.clf()

mags2 = numpy.asarray(mags)
phases2 = numpy.asarray(phases)
freqs2 = np.zeros(len(freqs) + 1)
fdiff = freqs[1]-freqs[0]
powers2 = np.zeros(len(powers)+1)
pdiff = powers[1] - powers[0]

freqs2[0:-1] = freqs-fdiff/2
freqs2[-1] = freqs[-1] + fdiff/2

powers2[0:-1] = powers-pdiff/2
powers2[-1] = powers[-1] + pdiff/2

ax = plt.subplot(1,2,1)
XX2,YY2 = np.meshgrid(freqs2,powers2+HWattenuation)
plt.pcolormesh(XX2/1e9,YY2,mags2, shading = 'nearest')



#ax.set_yaxis_ticks(powers, powers)

plt.xlabel("Frequency (GHz)")
plt.ylabel(r"Power (dBm)")
plt.title('S21 mag, {}'.format(scanname))
##cbar = fig.colorbar(im,cax=axins,orientation="horizontal",ticks=plt.MaxNLocator(2))
#cbar = plt.colorbar()
##cbar.set_label(r"$S_{21}$ (dB)",labelpad = -10,y=1,x=-0.35)


ax = plt.gca()
ax.xaxis.set_major_locator(plt.MaxNLocator(5))
ax.yaxis.set_major_locator(plt.MaxNLocator(5))
axins = inset_axes(ax,
                    width="20%",  # width = 50% of parent_bbox width
                    height="2%",  # height : 5%
                    loc='upper right',
                    bbox_to_anchor=(-0.02, -0.12, 1.0, 1.05),
                    bbox_transform=ax.transAxes,
                    borderpad=0,)
cbar = fig.colorbar(im,cax=axins,orientation="horizontal",ticks=plt.MaxNLocator(2))
cbar.set_label(r"$S_{21}$ (dB)",labelpad = -10,y=1,x=-0.35)
axins.xaxis.set_ticks_position("top")
axins.tick_params(labelsize=9)

ax = plt.subplot(1,2,2)
XX2,YY2 = np.meshgrid(freqs2,powers2+HWattenuation)
plt.pcolormesh(XX2/1e9,YY2,phases2, shading = 'nearest')

plt.xlabel("Frequency (GHz)")
plt.ylabel(r"Power (dBm)")
plt.title('S21 phase, {}'.format(scanname))

ax = plt.gca()
ax.xaxis.set_major_locator(plt.MaxNLocator(5))
ax.yaxis.set_major_locator(plt.MaxNLocator(5))
axins = inset_axes(ax,
                    width="20%",  # width = 50% of parent_bbox width
                    height="2%",  # height : 5%
                    loc='upper right',
                    bbox_to_anchor=(-0.02, -0.12, 1.0, 1.05),
                    bbox_transform=ax.transAxes,
                    borderpad=0,)
cbar = fig.colorbar(im,cax=axins,orientation="horizontal",ticks=plt.MaxNLocator(2))
cbar.set_label(r"$S_{21}$ (dB)",labelpad = -10,y=1,x=-0.35)
axins.xaxis.set_ticks_position("top")
axins.tick_params(labelsize=9)

plt.show()

userfuncs.SaveFull(saveDir, scanname, ['mags', 'phases', 'freqs', 'powers'], locals(), expsettings=initsettings)
plt.savefig(os.path.join(saveDir, scanname+'.png'), dpi = 150)


saveDir = r'Z:\Data\HouckTaTransmon\20201109'

stamp = userfuncs.timestamp()

name = 'RFscanFine'

scanname = name + '_' + stamp

settings = {}
settings['channel'] = 1
settings['averages'] = 80
settings['measurement'] = 'S21'
settings['meas_form'] = 'MLOG'
settings['start'] = 5.05e9
settings['stop'] = 5.1e9
settings['sweep_points'] = 1001
settings['RFpower'] = -5
settings['RFport'] = 3
settings['Mport'] = 2
settings['CAVport'] = 1
settings['CAVpower'] = -15
settings['CAVfreq'] = 7.173378e9
settings['ifBW'] = 500

powers = numpy.linspace(-15,-5,25)

HWattenuation = -30
initsettings = settings

mags = []
phases = []
for power in powers:
    print('Power: {}, final power: {}'.format(power, powers[-1]))
    initsettings['RFpower'] = power
    vna.spec_meas(settings)
    m = vna.get_trace(1,'spec')
    p = vna.get_trace(1,'phase')
    mags.append(m)
    phases.append(p)
    
freqs = vna.get_channel_axis(1)

fig2 = plt.figure(9)
plt.clf()

mags2 = numpy.asarray(mags)
phases2 = numpy.asarray(phases)
freqs2 = np.zeros(len(freqs) + 1)
fdiff = freqs[1]-freqs[0]
powers2 = np.zeros(len(powers)+1)
pdiff = powers[1] - powers[0]

freqs2[0:-1] = freqs-fdiff/2
freqs2[-1] = freqs[-1] + fdiff/2

powers2[0:-1] = powers-pdiff/2
powers2[-1] = powers[-1] + pdiff/2

ax = plt.subplot(1,2,1)
XX2,YY2 = np.meshgrid(freqs2,powers2+HWattenuation)
plt.pcolormesh(XX2/1e9,YY2,mags2, shading = 'nearest')



#ax.set_yaxis_ticks(powers, powers)

plt.xlabel("Frequency (GHz)")
plt.ylabel(r"Power (dBm)")
plt.title('S21 mag, {}'.format(scanname))
##cbar = fig.colorbar(im,cax=axins,orientation="horizontal",ticks=plt.MaxNLocator(2))
#cbar = plt.colorbar()
##cbar.set_label(r"$S_{21}$ (dB)",labelpad = -10,y=1,x=-0.35)


ax = plt.gca()
ax.xaxis.set_major_locator(plt.MaxNLocator(5))
ax.yaxis.set_major_locator(plt.MaxNLocator(5))
axins = inset_axes(ax,
                    width="20%",  # width = 50% of parent_bbox width
                    height="2%",  # height : 5%
                    loc='upper right',
                    bbox_to_anchor=(-0.02, -0.12, 1.0, 1.05),
                    bbox_transform=ax.transAxes,
                    borderpad=0,)
cbar = fig.colorbar(im,cax=axins,orientation="horizontal",ticks=plt.MaxNLocator(2))
cbar.set_label(r"$S_{21}$ (dB)",labelpad = -10,y=1,x=-0.35)
axins.xaxis.set_ticks_position("top")
axins.tick_params(labelsize=9)

ax = plt.subplot(1,2,2)
XX2,YY2 = np.meshgrid(freqs2,powers2+HWattenuation)
plt.pcolormesh(XX2/1e9,YY2,phases2, shading = 'nearest')

plt.xlabel("Frequency (GHz)")
plt.ylabel(r"Power (dBm)")
plt.title('S21 phase, {}'.format(scanname))

ax = plt.gca()
ax.xaxis.set_major_locator(plt.MaxNLocator(5))
ax.yaxis.set_major_locator(plt.MaxNLocator(5))
axins = inset_axes(ax,
                    width="20%",  # width = 50% of parent_bbox width
                    height="2%",  # height : 5%
                    loc='upper right',
                    bbox_to_anchor=(-0.02, -0.12, 1.0, 1.05),
                    bbox_transform=ax.transAxes,
                    borderpad=0,)
cbar = fig.colorbar(im,cax=axins,orientation="horizontal",ticks=plt.MaxNLocator(2))
cbar.set_label(r"$S_{21}$ (dB)",labelpad = -10,y=1,x=-0.35)
axins.xaxis.set_ticks_position("top")
axins.tick_params(labelsize=9)

plt.show()

userfuncs.SaveFull(saveDir, scanname, ['mags', 'phases', 'freqs', 'powers'], locals(), expsettings=initsettings)
plt.savefig(os.path.join(saveDir, scanname+'.png'), dpi = 150)
########################################################################################################
#stamp = userfuncs.timestamp()
#
#name = 'CAVpower_scan'
#
#CAVscanname = name + '_' + stamp
#
#powers = numpy.linspace(-20,-10,50)
#
#initsettings = settings
#
#mags = []
#
#for power in powers:
#    initsettings['CAVpower'] = power
#    vna.spec_meas(settings)
#    m = vna.get_trace(1,'spec')
#    mags.append(m)
#    
#freqs = vna.get_channel_axis(1)
#
#fig = plt.figure(dpi=150)
#XX,YY = np.meshgrid(freqs,powers)
#im = plt.pcolormesh(XX/1e9,YY,mags,cmap="viridis")
#plt.xlabel("Frequency (GHz)")
#plt.ylabel(r"Power (dBm)")
#plt.title('Cavity power scan, {}'.format(CAVscanname))
#ax = plt.gca()
#ax.xaxis.set_major_locator(plt.MaxNLocator(5))
#ax.yaxis.set_major_locator(plt.MaxNLocator(5))
#axins = inset_axes(ax,
#                    width="20%",  # width = 50% of parent_bbox width
#                    height="2%",  # height : 5%
#                    loc='lower right',
#                    bbox_to_anchor=(-0.02, -0.12, 1.0, 1.05),
#                    bbox_transform=ax.transAxes,
#                    borderpad=0,)
#cbar = fig.colorbar(im,cax=axins,orientation="horizontal",ticks=plt.MaxNLocator(2))
#cbar.set_label(r"$S_{21}$ (dB)",labelpad = -10,y=1,x=-0.35)
#axins.xaxis.set_ticks_position("top")
#axins.tick_params(labelsize=9)
#
#userfuncs.SaveFull(saveDir, CAVscanname, ['mags', 'freqs', 'powers'], locals(), expsettings=initsettings)
#
#plt.savefig(os.path.join(saveDir, CAVscanname+'.png'), dpi = 150)