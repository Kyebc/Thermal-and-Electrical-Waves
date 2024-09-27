# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 22:27:03 2023

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
from matplotlib import cbook
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
#%%
'''
importing data from excel
'''

df = pd.read_excel("Combined2.xlsx", engine = "openpyxl")
df = df.sort_values(by = ['frequency'])
data_array = np.array(df)
datind, amp1, amp2, freq, phs, datset = data_array[:, 0], data_array[:, 1], data_array[:, 2], data_array[:, 3], data_array[:, 4], data_array[:, 5]
# amprat = amp1/amp2
amprat = amp2/amp1

'''
data contained in datind, amp1, amp2, freq, phs, datset, ratamp respectively
'''
#%%

# plt.figure()
# plt.plot(np.arange(0, len(freq)), freq)
# plt.show()

hold = 0
counter = 0


in_freq = []
in_ratio = []
in_phase = []
in_amp1 = []
in_amp2 = []

out_freq = []
out_ratio = []
out_phase = []
out_amp1 = []
out_amp2 = []

'''
'''
in_flag = False
out_flag = False
error_tolerance = 5.0 # in degrees (angles)

in_stat_amp = []
in_stat_freq = []
in_stat_phs = []
out_stat_amp = []
out_stat_freq = []
out_stat_phs = []


in_stat_amp_val = []
in_stat_amp_err = []
in_stat_freq_val = []
in_stat_freq_err = []
in_stat_phs_val = []
in_stat_phs_err = []
out_stat_amp_val = []
out_stat_amp_err = []
out_stat_freq_val = []
out_stat_freq_err = []
out_stat_phs_val = []
out_stat_phs_err = []
in_count = 1
out_count = 0

for i in range(0, len(datind)):
    if phs[i] <= (-180 + error_tolerance) or phs[i] >= (180 - error_tolerance):
        out_flag = True
        out_freq.append(freq[i]*2*np.pi)
        out_ratio.append(amprat[i])
        out_phase.append(phs[i])
        out_amp1.append(amp1[i])
        out_amp2.append(amp2[i])
        out_stat_amp.append(amprat[i])
        out_stat_freq.append(freq[i])
        # out_stat_phs.append(phs[i])
        if phs[i] > 0:
            out_stat_phs.append(phs[i])
            pass
        elif phs[i] <= 0:
            out_stat_phs.append(360 + phs[i])
        if in_flag == True:
            in_stat_amp_val.append(np.mean(in_stat_amp))
            in_stat_amp_err.append(np.std(in_stat_amp))
            in_stat_freq_val.append(np.mean(in_stat_freq))
            in_stat_freq_err.append(np.std(in_stat_freq))
            in_stat_phs_val.append((((np.mean(np.array(in_stat_phs) + 360)) - 360) + in_count*360)/360)
            in_stat_phs_err.append(np.std(in_stat_phs)/360)
            in_count += 1
            in_stat_amp = []
            in_stat_freq = []
            in_stat_phs = []
            in_flag = False
            pass
        pass
    else:
        # if out_flag == True:
        #     print("out")
        #     out_stat_amp_val.append(np.mean(out_stat_amp))
        #     out_stat_amp_err.append(np.std(out_stat_amp))
        #     out_stat_freq_val.append(np.mean(out_stat_freq))
        #     out_stat_freq_err.append(np.std(out_stat_freq))
        #     out_stat_amp = []
        #     out_stat_freq = []
        #     out_flag = False
        pass

    if phs[i] >= (-1*error_tolerance) and phs[i] <= (error_tolerance): 
        in_flag = True
        in_freq.append(freq[i]*2*np.pi)
        in_ratio.append(amprat[i])
        in_phase.append(phs[i])
        in_amp1.append(amp1[i])
        in_amp2.append(amp2[i])
        in_stat_amp.append(amprat[i])
        in_stat_freq.append(freq[i])
        in_stat_phs.append(phs[i])
        if out_flag == True:
            out_stat_amp_val.append(np.mean(out_stat_amp))
            out_stat_amp_err.append(np.std(out_stat_amp))
            out_stat_freq_val.append(np.mean(out_stat_freq))
            out_stat_freq_err.append(np.std(out_stat_freq))
            out_stat_phs_val.append(((np.mean(np.array(out_stat_phs))) + (out_count*360))/360)
            out_stat_phs_err.append(np.std(out_stat_phs)/360)
            out_count += 1
            out_stat_amp = []
            out_stat_freq = []
            out_stat_phs = []
            out_flag = False
            pass
        pass
    else:
        # if in_flag == True:
        #     print("in")
        #     in_stat_amp_val.append(np.mean(in_stat_amp))
        #     in_stat_amp_err.append(np.std(in_stat_amp))
        #     in_stat_freq_val.append(np.mean(in_stat_freq))
        #     in_stat_freq_err.append(np.std(in_stat_freq))
        #     in_stat_amp = []
        #     in_stat_freq = []
        #     in_flag = False
        pass
    pass


# plt.figure()
# plt.scatter(in_stat_phs_val, np.array(in_stat_freq_val)*2*np.pi)
# plt.scatter(out_stat_phs_val, np.array(out_stat_freq_val)*2*np.pi)
# plt.title("")
# plt.xlabel("")
# plt.ylabel("")
# plt.show()

# plt.figure()
# plt.scatter(np.array(in_stat_freq_val)*2*np.pi, in_stat_phs_val)
# plt.scatter(np.array(out_stat_freq_val)*2*np.pi, out_stat_phs_val)
# plt.title("")
# plt.xlabel("")
# plt.ylabel("")
# plt.show()

#%%

'''Task 3.7'''

# plt.figure()
# plt.scatter(np.array(in_freq)/1000, in_ratio, color = "black")
# plt.xlabel("Angular Frequency (\u03C9)/ $\mathregular{krad^{-1}}$")
# plt.ylabel("Amplitude Ratio (In)")
# plt.title("Amplitude (In Phases) Against Frequency")
# plt.grid()
# plt.show()

# #%%

# plt.figure()
# plt.scatter(np.array(out_freq)/1000, out_ratio)
# plt.xlabel("Angular Frequency (\u03C9)/ $\mathregular{krad^{-1}}$")
# plt.ylabel("Amplitude Ratio (Out)")
# plt.title("Amplitude (Out Phases) Against Frequency")
# plt.grid()
# plt.show()

# #%%

# plt.figure()
# plt.scatter(np.array(freq)*2*np.pi/1000, amprat, marker = ".", linewidths = 0.01)

# plt.xlabel("Angular Frequency (\u03C9)/ $\mathregular{krad^{-1}}$")
# plt.ylabel("Amplitude Ratio (Out/In)")
# plt.title("Ratio of Output Over Input (Amplitude) Against Frequency")
# plt.grid()
# plt.show()

#%%

# output amplitude against frequency

# last_four_freq = []
# last_four_amp = []
# last_four_freq.append(in_stat_freq_val[-1]*2*np.pi/1000)
# last_four_amp.append(in_stat_amp_val[-1])
# last_four_freq.append(in_stat_freq_val[-2]*2*np.pi/1000)
# last_four_amp.append(in_stat_amp_val[-2])
# last_four_freq.append(out_stat_freq_val[-1]*2*np.pi/1000)
# last_four_amp.append(out_stat_amp_val[-1])
# last_four_freq.append(out_stat_freq_val[-2]*2*np.pi/1000)
# last_four_amp.append(out_stat_amp_val[-2])

# def line(x, m, c):
#     return m*x + c

# popt, pcov = curve_fit(line, last_four_freq, last_four_amp)

''''''

# freqcombined = in_stat_freq_val + out_stat_freq_val
# ampcombined = in_stat_amp_val +out_stat_amp_val

# def exp(x, alpha, beta, gamma):
#     return alpha - beta*np.exp(gamma*x)

# popt, pcov = curve_fit(exp, freqcombined, ampcombined, p0 = [1.932, ])

# fitfreq = []
# fitmag = []
# fitmagmin = []
# fitmagmax = []
# a = popt[0]
# b = popt[1]
# c = popt[2]
# amin = popt[0] - np.sqrt(pcov[0, 0])
# amax = popt[0] + np.sqrt(pcov[0, 0])
# bmax = popt[1] - np.sqrt(pcov[1, 1])
# bmin = popt[1] + np.sqrt(pcov[1, 1])
# cmax = popt[2] - np.sqrt(pcov[2, 2])
# cmin = popt[2] + np.sqrt(pcov[2, 2])
# for i in range(-100, 1100):
#     fitfreq.append(i)
#     fitmag.append(exp(i, a, b, c))
#     fitmagmin.append(exp(i, amin, bmax, cmax))
#     fitmagmax.append(exp(i, amax, bmin, bmin))
    

''''''

freqcombined = in_stat_freq_val[-3: -1] + out_stat_freq_val[-3: -1]
ampcombined = in_stat_amp_val[-3: -1] + out_stat_amp_val[-3: -1]

# freqcombined = freq[-3999: -1]
# ampcombined = amprat[-3999: -1]

offset = min(freqcombined)
freqcombined = freqcombined - offset

def line(x, m, c):
    return m*x + c

popt, pcov = curve_fit(line, freqcombined, ampcombined)

m = popt[0]
c = popt[1]
mmin = popt[0] - np.sqrt(pcov[0, 0])
mmax = popt[0] + np.sqrt(pcov[0, 0])
cmin = popt[1] - np.sqrt(pcov[1, 1])
cmax = popt[1] + np.sqrt(pcov[1, 1])

fitfreq = []
fitmag = []
fitmagmin = []
fitmagmax = []

for i in range(int(min(freqcombined)), 200000):
    fitfreq.append(i + offset)
    fitmag.append(line(i, m, c))
    fitmagmin.append(line(i, mmin, cmin))
    fitmagmax.append(line(i, mmax, cmax))


cutoff = -1*(c/m) + offset
cutoffmin = -1*(cmin/mmin) + offset
cutoffmax = -1*(cmax/mmax) + offset

diff = (cutoffmax - cutoffmin)/2

cutoff = cutoff*2*np.pi/1000
diff = diff*2*np.pi/1000
''''''

fig, ax = plt.subplots()
plt.style.use('default')

ax.scatter(np.array(freq)*2*np.pi/1000, amprat, marker = ".", s = 1, label = "All Measurements")
ax.errorbar(np.array(in_stat_freq_val)*2*np.pi/1000, in_stat_amp_val, yerr = np.array(in_stat_amp_err), xerr = np.array(in_stat_freq_err)/1000, fmt = "none", capsize = 2, color = 'red')
ax.errorbar(np.array(out_stat_freq_val)*2*np.pi/1000, out_stat_amp_val, yerr = out_stat_amp_err, xerr = np.array(out_stat_freq_err)/1000, fmt = "none", capsize = 2, color = 'red')
ax.scatter(np.array(in_stat_freq_val)*2*np.pi/1000, in_stat_amp_val, marker = "o", linewidths = 0.05, label = "In Phase")
ax.scatter(np.array(out_stat_freq_val)*2*np.pi/1000, out_stat_amp_val, marker = "o", linewidths = 0.05, label = "Out Phase", color = "purple")


plt.xlim(-20, 950)
plt.ylim(-0.1, 1.1)
plt.xlabel("Angular Frequency (\u03C9)/ $\mathregular{k(rad)s^{-1}}$")
plt.ylabel("Amplitude Ratio (Out/In)")
plt.title("Ratio of Output and Input Amplitude Against Frequency")
plt.grid()

ax.plot(np.array(fitfreq)*2*np.pi/1000, fitmag, linestyle = "dashed", color = 'black', label = "Linear Projection")
# ax.plot(np.array(fitfreq)*2*np.pi/1000, fitmagmin)
# ax.plot(np.array(fitfreq)*2*np.pi/1000, fitmagmax)
ax.fill_between(np.array(fitfreq)*2*np.pi/1000, fitmagmin, fitmagmax, color = "grey", alpha = 0.25, label = "One Standard Deviation")
ax.text(650, 0.02, "$\omega_c$= %s $\pm$ %s" % (int(cutoff), int(diff)))
ax.legend(loc = 1)


plt.style.use('default')
axins = zoomed_inset_axes(ax, 50, loc = 3) # zoom = 6
axins.grid()
mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")

axins.scatter(np.array(freq)*2*np.pi/1000, amprat, marker = ".", s = 1, label = "All Measurements")
axins.errorbar(np.array(in_stat_freq_val)*2*np.pi/1000, in_stat_amp_val, yerr = np.array(in_stat_amp_err), xerr = np.array(in_stat_freq_err)/1000, fmt = "none", capsize = 2, color = 'red')
axins.errorbar(np.array(out_stat_freq_val)*2*np.pi/1000, out_stat_amp_val, yerr = out_stat_amp_err, xerr = np.array(out_stat_freq_err)/1000, fmt = "none", capsize = 2, color = 'red')
axins.scatter(np.array(in_stat_freq_val)*2*np.pi/1000, in_stat_amp_val, marker = "o", linewidths = 0.05, label = "In Phase")
axins.scatter(np.array(out_stat_freq_val)*2*np.pi/1000, out_stat_amp_val, marker = "o", linewidths = 0.05, label = "Out Phase", color = "purple")

axins.set_xlim(518.2, 522.5) # Limit the region for zoom
axins.set_ylim(0.490, 0.499)
axins.yaxis.tick_right()
axins.xaxis.tick_top()
axins.text(521, 0.498, "50x")
axins.set_facecolor('white')


plt.style.use('default')
axins2 = zoomed_inset_axes(ax, 3, loc = 7) # zoom = 6
axins2.grid()
mark_inset(ax, axins2, loc1=3, loc2=4, fc="none", ec="0.5")

axins2.scatter(np.array(freq)*2*np.pi/1000, amprat, marker = ".", s = 1, label = "All Measurements")
axins2.errorbar(np.array(in_stat_freq_val)*2*np.pi/1000, in_stat_amp_val, yerr = np.array(in_stat_amp_err), xerr = np.array(in_stat_freq_err)/1000, fmt = "none", capsize = 2, color = 'red')
axins2.errorbar(np.array(out_stat_freq_val)*2*np.pi/1000, out_stat_amp_val, yerr = out_stat_amp_err, xerr = np.array(out_stat_freq_err)/1000, fmt = "none", capsize = 2, color = 'red')
axins2.scatter(np.array(in_stat_freq_val)*2*np.pi/1000, in_stat_amp_val, marker = "o", linewidths = 0.05, label = "In Phase")
axins2.scatter(np.array(out_stat_freq_val)*2*np.pi/1000, out_stat_amp_val, marker = "o", linewidths = 0.05, label = "Out Phase", color = "purple")
axins2.plot(np.array(fitfreq)*2*np.pi/1000, fitmag, linestyle = "dashed", color = 'black', label = "Linear Projection")
axins2.fill_between(np.array(fitfreq)*2*np.pi/1000, fitmagmin, fitmagmax, color = "grey", alpha = 0.25, label = "One Standard Deviation")

axins2.set_xlim(840, 890) # Limit the region for zoom
axins2.set_ylim(-0.02, 0.09)
# axins2.yaxis.tick_right()
# axins2.xaxis.tick_top()
axins2.text(521, 0.498, "50x")
axins2.set_facecolor('white')


plt.show()

#%%

#Input amplitude against frequency

fig, ax = plt.subplots()
plt.style.use('default')

ax.scatter(np.array(freq)*2*np.pi/1000, amp1, marker = ".", s = 1, label = "All Datapoints")
ax.scatter(np.array(in_freq)/1000, in_amp1, marker = "o", linewidths = 0.05, label = "In Phase Datapoints")
ax.scatter(np.array(out_freq)/1000, out_amp1, marker = "o", linewidths = 0.05, label = "Out Phase Datapoints", color = 'purple')
ax.legend(loc = 2)

plt.xlim(-20, 920)
plt.ylim(1.750, 3.500)
plt.xlabel("Angular Frequency (\u03C9)/ $\mathregular{k(rad)s^{-1}}$")
plt.ylabel("Amplitude/ V")
plt.title("Input Amplitude Against Frequency")
ax.grid()


plt.style.use('default')
axins = zoomed_inset_axes(ax, 9, loc = 6) # zoom = 6
axins.grid()
mark_inset(ax, axins, loc1=1, loc2=4, fc="none", ec="0.5")

axins.scatter(np.array(freq)*2*np.pi/1000, amp1, marker = ".", s = 1, label = "All Datapoints")
axins.scatter(np.array(in_freq)/1000, in_amp1, marker = "o", linewidths = 0.05, label = "In Phase Datapoints")
axins.scatter(np.array(out_freq)/1000, out_amp1, marker = "o", linewidths = 0.05, label = "Out Phase Datapoints", color = 'purple')

axins.set_xlim(625, 645) # Limit the region for zoom
axins.set_ylim(2.475, 2.545)
axins.yaxis.tick_right()
# axins.xaxis.tick_top()
axins.text(628, 2.53, "9x")
axins.set_facecolor('white')

plt.show()

#%%

# output amplitude against frequency

# last_four_freq = []
# last_four_amp = []
# last_four_freq.append(in_freq[-1]/1000)
# last_four_amp.append(in_amp2[-1])
# last_four_freq.append(in_freq[-2]/1000)
# last_four_amp.append(in_amp2[-2])
# last_four_freq.append(out_freq[-1]/1000)
# last_four_amp.append(out_amp2[-1])
# last_four_freq.append(out_freq[-2]/1000)
# last_four_amp.append(out_amp2[-2])

# def line(x, m, c):
#     return m*x + c

# fitfreq = []
# fitmag = []
# fitmagmin = []
# fitmagmax = []
# popt, pcov = curve_fit(line, last_four_freq, last_four_amp)
# m = popt[0]
# c = popt[1]
# mmin = popt[0] - np.sqrt(pcov[0, 0])
# mmax = popt[0] + np.sqrt(pcov[0, 0])
# cmax = popt[1] - np.sqrt(pcov[1, 1])
# cmin = popt[1] + np.sqrt(pcov[1, 1])
# for i in range((int(min(last_four_freq))), 1100):
#     fitfreq.append((min(last_four_freq)))
#     fitmag.append(line(i, m, c))
#     fitmagmin.append(line(i, mmin, cmin))
#     fitmagmax.append(line(i, mmax, cmax))
    



fig, ax = plt.subplots()
plt.style.use('default')
ax.scatter(np.array(freq)*2*np.pi/1000, amp2, marker = ".", s = 1, label = "All Datapoints")
ax.scatter(np.array(in_freq)/1000, in_amp2, marker = "o", linewidths = 0.05, label = "In Phase Datapoints")
ax.scatter(np.array(out_freq)/1000, out_amp2, marker = "o", linewidths = 0.05, label = "Out Phase Datapoints", color = "purple")
ax.legend(loc = 1)

# ax.plot(fitfreq, fitmag)
# ax.plot(fitfreq, fitmagmin)
# ax.plot(fitfreq, fitmagmax)

plt.xlim(-20, 920)
plt.ylim(-0.10, 2.0)
ax.grid()
plt.xlabel("Angular Frequency (\u03C9)/ $\mathregular{k(rad)s^{-1}}$")
plt.ylabel("Amplitude/ V")
plt.title("Output Amplitude Against Frequency")

plt.style.use('default')
axins = zoomed_inset_axes(ax, 9, loc = 3) # zoom = 6
axins.grid()
mark_inset(ax, axins, loc1=1, loc2=2, fc="none", ec="0.5")
axins.scatter(np.array(freq)*2*np.pi/1000, amp2, marker = ".", s = 1, label = "All Datapoints")
axins.scatter(np.array(in_freq)/1000, in_amp2, marker = "o", linewidths = 0.05, label = "In Phase Datapoints")
axins.scatter(np.array(out_freq)/1000, out_amp2, marker = "o", linewidths = 0.05, label = "Out Phase Datapoints", color = "purple")
axins.set_xlim(445, 505) # Limit the region for zoom
axins.set_ylim(1.170, 1.260)
axins.yaxis.tick_right()
axins.xaxis.tick_top()
axins.text(453, 1.185, "x9")
axins.set_facecolor('white')

plt.show()

#%%

'''Task 3.8'''

#3.8a

positive = False
phs2 = []
freq2 = []
loop = 0
countdown = 10
count = countdown
count2 = -countdown
for i in range(0, len(freq)):
    if phs[i] > 0:
        if count2 >= 0:
            phs2.append(phs[i]/360 + loop - 1)
            pass
        else:
            phs2.append(phs[i]/360 + loop)
            pass
        freq2.append(freq[i]*2*np.pi)
        positive = True
        count -= 1
    if phs[i] <= 0:
        if count <= 0:
            loop += 1
            count2 = countdown
            pass
        phs2.append(phs[i]/360 + loop)
        freq2.append(freq[i]*2*np.pi)
        positive = False
        count = countdown
        count2 -= 1
        pass

def func(x, m, c):
    return m*x + c

fit_range = 0.50
popt, pcov = curve_fit(func, np.array(freq2[0:int(fit_range*len(freq2))]), np.array(phs2[0:int(fit_range*len(freq2))]))

yfit = []
ylow = []
yhigh = []
for i in range(0, len(freq2)):
    yfit.append(func(freq2[i], popt[0], popt[1]))
    ylow.append(-999)
    yhigh.append(999)
    pass


# plt.figure()
# plt.grid()
# plt.scatter(np.array(freq2)/1000, np.array(phs2)*2*np.pi/40, label = "All Datapoints", s = 1)
# plt.scatter(np.array(in_stat_freq_val)*2*np.pi/1000, np.array(in_stat_phs_val)*2*np.pi/40, label = "In Phase")
# plt.errorbar(np.array(in_stat_freq_val)*2*np.pi/1000, np.array(in_stat_phs_val)*2*np.pi/40, np.array(in_stat_phs_err)*2*np.pi/1000, np.array(in_stat_freq_err)*2*np.pi/40, capsize = 2, fmt = 'none', color = 'red')
# plt.scatter(np.array(out_stat_freq_val)*2*np.pi/1000, np.array(out_stat_phs_val)*2*np.pi/40, label = "Out Phase", color = "purple")
# plt.errorbar(np.array(out_stat_freq_val)*2*np.pi/1000, np.array(out_stat_phs_val)*2*np.pi/40, np.array(out_stat_phs_err)*2*np.pi/1000, np.array(out_stat_freq_err)*2*np.pi/40, capsize = 2, fmt = 'none', color = 'red')
# plt.plot(np.array(freq2)/1000, np.array(yfit)*2*np.pi/40, linestyle = "dashed", color = "black", alpha = 0.50, label = "Linear Relationship")
# plt.fill_between([-999, int(freq2[int(fit_range*len(freq2))])/1000], [-999, -999], [999, 999], color = "grey", alpha = 0.25, label = "Linear Region")
# plt.title("Angular Frequency Against Phase")
# plt.ylabel("Phase (k)/ $\mathregular{section^{-1}}$")
# plt.xlabel("Angular Frequency (\u03C9)/ $\mathregular{krad^{-1}}$")
# plt.xlim(-40, 950)
# plt.ylim(-0.1, 2.75)
# plt.legend()
# plt.show()

#%%
plt.style.use('default')
fig, ax = plt.subplots()
ax.grid()
ax.scatter(np.array(freq2)/1000, np.array(phs2)*2*np.pi/40, label = "All Datapoints", s = 1)
ax.scatter(np.array(in_stat_freq_val)*2*np.pi/1000, np.array(in_stat_phs_val)*2*np.pi/40, label = "In Phase")
# ax.errorbar(np.array(in_stat_freq_val)*2*np.pi/1000, np.array(in_stat_phs_val)*2*np.pi/40, np.array(in_stat_phs_err)*2*np.pi/40, np.array(in_stat_freq_err)*2*np.pi/1000, capsize = 2, fmt = 'none', color = 'red')
ax.scatter(np.array(out_stat_freq_val)*2*np.pi/1000, np.array(out_stat_phs_val)*2*np.pi/40, label = "Out Phase", color = "purple")
# ax.errorbar(np.array(out_stat_freq_val)*2*np.pi/1000, np.array(out_stat_phs_val)*2*np.pi/40, np.array(out_stat_phs_err)*2*np.pi/40, np.array(out_stat_freq_err)*2*np.pi/1000, capsize = 2, fmt = 'none', color = 'red')
ax.plot(np.array(freq2)/1000, np.array(yfit)*2*np.pi/40, linestyle = "dashed", color = "black", alpha = 0.50, label = "Linear Relationship")
ax.fill_between([-999, int(freq2[int(fit_range*len(freq2))])/1000], [-999, -999], [999, 999], color = "grey", alpha = 0.25, label = "Linear Region")
plt.xlim(-20, 900)
plt.ylim(-0.10, 3.70)
plt.title("Phase Against Frequency")
plt.ylabel("Phase (k)/ $\mathregular{k section^{-1}}$")
plt.xlabel("Angular Frequency (\u03C9)/ $\mathregular{k(rad)s^{-1}}$")

axins = zoomed_inset_axes(ax, 200, loc = 2) # zoom = 6
axins.grid()
mark_inset(ax, axins, loc1=1, loc2=3, fc="none", ec="0.5")
axins.scatter(np.array(freq2)/1000, np.array(phs2)*2*np.pi/40, label = "All Datapoints", s = 1)
axins.scatter(np.array(in_stat_freq_val)*2*np.pi/1000, np.array(in_stat_phs_val)*2*np.pi/40, label = "In Phase")
axins.errorbar(np.array(in_stat_freq_val)*2*np.pi/1000, np.array(in_stat_phs_val)*2*np.pi/40, np.array(in_stat_phs_err)*2*np.pi/40, np.array(in_stat_freq_err)*2*np.pi/1000, capsize = 2, fmt = 'none', color = 'red')
axins.scatter(np.array(out_stat_freq_val)*2*np.pi/1000, np.array(out_stat_phs_val)*2*np.pi/40, label = "Out Phase", color = "purple")
axins.errorbar(np.array(out_stat_freq_val)*2*np.pi/1000, np.array(out_stat_phs_val)*2*np.pi/40, np.array(out_stat_phs_err)*2*np.pi/40, np.array(out_stat_freq_err)*2*np.pi/1000, capsize = 2, fmt = 'none', color = 'red')
axins.plot(np.array(freq2)/1000, np.array(yfit)*2*np.pi/40, linestyle = "dashed", color = "black", alpha = 0.50, label = "Linear Relationship")
# axins.fill_between([-999, int(freq2[int(fit_range*len(freq2))])/1000], [-999, -999], [999, 999], color = "grey", alpha = 0.25)
axins.set_xlim(392.1, 393.7) # Limit the region for zoom
axins.set_ylim(0.860, 0.868)
axins.yaxis.tick_right()
axins.text(392.4, 0.867, "Zoomed in 200x")
axins.set_facecolor('white')

plt.xticks(visible = True)  # Not present ticks
plt.yticks(visible = True)

ax.legend(loc = 1)
plt.draw()
plt.show()

#%%

# #3.8b

# gradient = []
# freq_g = []
# gradienterr = []
# freq_gerr = []

# xpool = []
# ypool = []
# newx = []
# newxerr = []
# newy = []
# newyerr = []
# for i in range(0, len(freq2)):
#     if i%500 != 0:
#         xpool.append(phs2[i])
#         ypool.append(freq2[i])
#         pass
#     else:
#         newx.append(np.mean(xpool))
#         newy.append(np.mean(ypool))
#         newxerr.append(np.std(xpool))
#         newyerr.append(np.std(ypool))
#         xpool = []
#         ypool = []
        
# for i in range(0, len(newy) - 1):
#     dellx = newx[i + 1] - newx[i]
#     delly = newy[i + 1] - newy[i]
#     if dellx == 0:
#         pass
#     else:
#         gradient.append(delly/dellx)
#         freq_g.append((newy[i] + newy[i + 1])/2)
#         pass
#     pass

# phsv = (np.array(freq2)/(2*np.pi))/np.array(phs2)
# # phsv = (np.array(freq2))/np.array(phs2)

# plt.figure()
# plt.grid()
# plt.scatter(np.array(freq_g)/1000, np.array(gradient)/(2*np.pi), label = "Group Velocity")
# plt.scatter(np.array(freq2)/1000, phsv, label = "Phase Velocity")
# plt.title("Group Velocity Against Frequency")
# plt.xlabel("Angular Frequency (\u03C9)/ $\mathregular{krad^{-1}}$")
# plt.ylabel("Velocity")
# plt.ylim(bottom = 0)
# plt.legend()
# plt.show()

#%%

'''3.8b'''

new_phs = np.concatenate((np.array(in_stat_phs_val)*2*np.pi/40, np.array(out_stat_phs_val)*2*np.pi/40), axis=None)
new_phs_err = np.concatenate((np.array(in_stat_phs_err)*2*np.pi/40, np.array(out_stat_phs_err)*2*np.pi/40), axis=None)
new_freq = np.concatenate((np.array(in_stat_freq_val)*2*np.pi/1000, np.array(out_stat_freq_val)*2*np.pi/1000), axis=None)
new_freq_err = np.concatenate((np.array(in_stat_freq_err)*2*np.pi/1000, np.array(out_stat_freq_err)*2*np.pi/1000), axis=None)

d = {'Phase': new_phs, 'Frequency': new_freq, 'Phase Error': new_phs_err, 'Frequency Error': new_freq_err}
df2 = pd.DataFrame(data = d)
df2 = df2.sort_values(by = ['Frequency'])
data_array = np.array(df2)
phs3, freq3, phs3err, freq3err = data_array[:, 0], data_array[:, 1], data_array[:, 2], data_array[:, 3]

# plt.figure()
# plt.scatter(phs3, freq3)
# plt.show()

gradient = []
freq_g = []
gradienterr = []
freq_gerr = []

for i in range(0, len(freq3) - 1):
    dellx = phs3[i + 1] - phs3[i]
    dellxerr = np.sqrt(phs3err[i + 1]**2 + phs3err[i]**2)
    delly = freq3[i + 1] - freq3[i]
    dellyerr = np.sqrt(freq3err[i + 1]**2 + freq3err[i]**2)
    gradient.append(delly/dellx)
    gradienterr.append((delly/dellx)*np.sqrt((dellxerr/dellx)**2 + (dellyerr/delly)**2))
    # gradienterr.append((delly/dellx)*np.sqrt((dellxerr)**2 + (dellyerr)**2))
    freq_g.append((freq3[i + 1] + freq3[i])/2)
    freq_gerr.append(np.sqrt(freq3err[i + 1]**2 + freq3err[i]**2))
    pass

phsv = np.array(freq3)/np.array(phs3)
phsverr = np.array(phsv)*np.sqrt((np.array(freq3err)/np.array(freq3))**2 + (np.array(phs3err)/np.array(phs3))**2)



# def exp(x, a, b, c):
#     return a + b*np.exp(c*x)

def exp1(x, b, c): #group vel
    return 485 + b*np.exp(c*x)

def exp2(x, b, c):
    return 472 + b*np.exp(c*x)

popt1, pcov1 = curve_fit(exp1, freq_g, gradient, p0 = [-15, 0.0018]) #group vel
popt2, pcov2 = curve_fit(exp2, freq3, phsv, p0 = [-12, 0.0025])





x1 = []
y1 = []
x1err = []
y1min = []
y1max = []
x2 = []
y2 = []
x2err = []
y2min = []
y2max = []
for i in range(0, 1000):
    x1.append(i)
    x2.append(i)
    y1.append(exp1(i, popt1[0], popt1[1]))
    y2.append(exp2(i, popt2[0], popt2[1]))
    if i >= int(freq2[int(fit_range*len(freq2))])/1000:
        x1err.append(i)
        y1min.append(exp1(i, popt1[0] - np.sqrt(pcov1[0, 0]), popt1[1] + np.sqrt(pcov1[1, 1])))
        y1max.append(exp1(i, popt1[0] + np.sqrt(pcov1[0, 0]), popt1[1] - np.sqrt(pcov1[1, 1])))
        x2err.append(i)
        y2min.append(exp2(i, popt2[0] - np.sqrt(pcov2[0, 0]), popt2[1] + np.sqrt(pcov2[1, 1])))
        y2max.append(exp2(i, popt2[0] + np.sqrt(pcov2[0, 0]), popt2[1] - np.sqrt(pcov2[1, 1])))
    pass


plt.figure()
plt.grid()

plt.errorbar(np.array(freq_g), np.array(gradient), gradienterr, freq_gerr, capsize = 2.5, fmt = 'none', color = 'black')
plt.scatter(np.array(freq_g), np.array(gradient), label = "Group Velocity", color = 'black')
plt.errorbar(np.array(freq3), phsv, np.array(phsverr), np.array(freq3err), capsize = 2.5, fmt = 'none', color = 'darkblue')
plt.scatter(np.array(freq3), phsv, label = "Phase Velocity", color = 'darkblue')

plt.plot(x1, y1, linestyle = 'dashed', color = 'red', label = "Group Velocity Fit")
plt.fill_between(x1err, y1min, y1max, color = "red", alpha = 0.10, label = "1$\sigma$")
# plt.plot(x1err, y1min)
# plt.plot(x1err, y1max)

plt.plot(x2, y2, linestyle = 'dashed', color = 'magenta', label = "Phase Velocity Fit")
plt.fill_between(x2err, y2min, y2max, color = "magenta", alpha = 0.10, label = "1$\sigma$")
# plt.plot(x2err, y2min)
# plt.plot(x2err, y2max)

plt.fill_between([-999, int(freq2[int(fit_range*len(freq2))])/1000], [-999, -999], [999, 999], color = "grey", alpha = 0.25, label = "Linear Region")
plt.title("Group and Phase Velocity Against Frequency")
plt.xlabel("Angular Frequency (\u03C9)/ $\mathregular{k(rad)s^{-1}}$")
plt.ylabel("Velocity/ m sections $s^{-1}$")
plt.xlim(-10, 950)
plt.ylim(200, 500)
plt.legend()
plt.show()





