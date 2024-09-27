# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 10:16:40 2023

@author: kyebchoo
"""

import numpy as np
import scipy as sci
import scipy.signal as sig
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os
import math
import pandas as pd
import pyvisa
import time
from matplotlib import cbook
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

#%%

# rm = pyvisa.ResourceManager()
# print(rm.list_resources())
# my_instrument = rm.open_resource('USB0::0x0AAD::0x01D6::110244::INSTR')
# my_instrument.query('*IDN?')

# #%%

# # my_instrument.write("WGENerator:FREQuency 20")
# # my_instrument.write("MEASurement1:RESult:AVG?")
# # print(my_instrument.read())

# my_instrument.write("WGENerator:VOLTage 4")
# my_instrument.write("WGENerator:VOLTage:OFFSet 0")

# my_instrument.write("CHANnel1:OFFSet 0")
# my_instrument.write("CHANnel2:OFFSet 0")

# my_instrument.write("CHANnel1:SCALe 750e-3")
# my_instrument.write("CHANnel2:SCALe 750e-3")

# my_instrument.write("TIMebase:SCALe " + str(1/20))

# #%%
# start_freq = 10000 #Hz
# end_freq   = 150000 #Hz
# pause_time = 5 # do not change
# # resolution = 100
# wavecount = 50
# intervals = 1000

# amplitude1 = []
# stddev1 = []

# amplitude2 = []
# stddev2 = []

# amplituderatio = []
# stddevratio = []

# phase = []
# stddevphase = []

# frequency = []
# stddevfrequency = []

# print("Estimated runtime:" +str((end_freq  - start_freq)*pause_time/intervals/60/60) + "hours")

# #%%
# for i in range(int(start_freq/intervals), int(end_freq/intervals) + 1):
#     my_instrument.write("WGENerator:FREQuency " + str(i*intervals))
#     my_instrument.write("TIMebase:SCALe " + str(1/(i*intervals)))
#     my_instrument.write("MEASurement1:STATistics:RESet")
#     my_instrument.write("MEASurement2:STATistics:RESet")
#     my_instrument.write("MEASurement3:STATistics:RESet")
#     my_instrument.write("MEASurement4:STATistics:RESet")
#     pause = True
#     while pause == True:
#         time.sleep(pause_time)
#         my_instrument.write("MEASurement3:RESult:WFMCount?")
#         wcount = int(my_instrument.read())
#         if wcount >= wavecount:
#             pause = False
#             pass
#         elif wcount < wavecount:
#             pass
#         pass
    
#     my_instrument.write("MEASurement1:RESult:AVG?")
#     amp1 = float(my_instrument.read())
#     amplitude1.append(amp1)
#     my_instrument.write("MEASurement1:RESult:STDDev?")
#     readstddev1 = float(my_instrument.read())
#     stddev1.append(readstddev1)
    
#     my_instrument.write("MEASurement2:RESult:AVG?")
#     amp2 = float(my_instrument.read())
#     amplitude2.append(amp2)
#     my_instrument.write("MEASurement2:RESult:STDDev?")
#     readstddev2 = float(my_instrument.read())
#     stddev2.append(readstddev2)
    
#     amplituderatio.append(amp2/amp1)
#     stddevratio = (amp2/amp1)*(np.sqrt(readstddev1**2 + readstddev2**2))
    
#     my_instrument.write("MEASurement3:RESult:AVG?")
#     frequency.append(float(my_instrument.read()))
#     my_instrument.write("MEASurement3:RESult:STDDev?")
#     stddevfrequency.append(float(my_instrument.read()))
    
#     my_instrument.write("MEASurement4:RESult:AVG?")
#     phase.append(float(my_instrument.read()))
#     my_instrument.write("MEASurement4:RESult:STDDev?")
#     stddevphase.append(float(my_instrument.read()))
    
#     percent = str(round((i - (start_freq/intervals))*100/((end_freq/intervals) - (start_freq/intervals)), 3))
#     print("%s percent completed." % (percent))
    
#     pass

plt.figure()
plt.plot(frequency, amplitude1)
plt.ylim(1, 4)
plt.show()

plt.figure()
plt.plot(frequency, amplitude2)
plt.ylim(0, 2)
plt.show()

plt.figure()
plt.plot(frequency, amplituderatio)
plt.ylim(0, 1)
plt.show()

plt.figure()
plt.plot(frequency, phase)
plt.ylim(-300, 300)
plt.show()

#%%

d = {'Amplitude 1': amplitude1, 'Amplitude 2': amplitude2, 'Amplitude Ratio': amplituderatio, 'Phase': phase, 'Frequency': frequency}
df = pd.DataFrame(data = d)
data_array = np.array(df)
amp1, amp2, amprat, phs, freq = data_array[:, 0], data_array[:, 1], data_array[:, 2], data_array[:, 3], data_array[:, 4]

count = 0
for i in range(0, len(freq)):
    if amp2[i - count] > 2:
        amp1 = np.delete(amp1, i - count)
        amp2 = np.delete(amp2, i - count)
        amprat = np.delete(amprat, i - count)
        phs = np.delete(phs, i - count)
        freq = np.delete(freq, i - count)
        count += 1
        pass
    
d2 = {'Amplitude 1': amp1, 'Amplitude 2': amp2, 'Amplitude Ratio': amprat, 'Phase': phs, 'Frequency': freq*2*np.pi}
df2 = pd.DataFrame(data = d2)

plt.figure()
plt.scatter(freq, phs)
plt.show()

phsnew = []
last = 'negative'
current = 'negative'
cycle = 0
for i in range(0, len(freq)):
    if phs[i] < 0:
        current = 'negative'
        pass
    if phs[i] >= 0:
        current = 'positive'
        pass
    if current == 'negative' and past == 'positive':
        cycle += 1
        pass
    phsnew.append()
    last = current
    pass

        
    






























