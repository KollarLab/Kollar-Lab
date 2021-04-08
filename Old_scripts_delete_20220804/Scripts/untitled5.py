# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 16:53:28 2021

@author: Kollarlab
"""

filedir = r'Z:\Data\HouckDualHangerFluxonium\drainage_T1\20210224'
file1 = 'qubit_off_20210224_171654.pkl'
file2 = 'qubit_on_5us_overlap_20210224_171448.pkl'

raw_data1 = userfuncs.LoadFull(os.path.join(filedir, file1))
raw_data2 = userfuncs.LoadFull(os.path.join(filedir, file2))

amp1 = raw_data1[0]['amp']
amp2 = raw_data2[0]['amp']
time = raw_data1[0]['xaxis_us']

plt.plot(time, amp2-amp1)
plt.xlabel('Time (us)')
plt.ylabel('Diff (V)')
plt.title('Difference between qubit pulse and no qubit pulse')