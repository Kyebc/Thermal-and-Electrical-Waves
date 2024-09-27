# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 10:45:26 2023

@author: kyebchoo
"""

import numpy as np
import scipy as sci
import scipy.signal as sig
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os
#%%

signal_amplitude = 50

def square_signal(reference_data_x_array, period):
    xsignal = reference_data_x_array
    ysignal = signal_amplitude*(1 + sig.square(xsignal*2*np.pi/period))
    return ysignal

xarray = []
yarray = []
period = 0.1
def plot_graph(file_name, period_sig, file_name2 = "False"):
    x1, y1 = np.loadtxt("%s" % (file_name), unpack = True, skiprows = 3)
    global xarray, yarray, period
    xarray = x1
    yarray = y1
    period = period_sig
    plt.figure()
    ysquare1 = square_signal(x1, period_sig)
    if file_name2 != "False":
        plt.title("Signal Plot")
        label1 = "%s" % (file_name)
        label2 = "%s" % (file_name2)
    else:
        plt.title("File name: \"%s\"" % (file_name))
        label1 = "Measurement"
    plt.plot(x1/10, ysquare1, color = "blue", label = "Signal")
    plt.plot(x1/10, y1, label = "%s" % (label1), color = "red")
    if file_name2 != "False":
        x2, y2 = np.loadtxt("%s" % (file_name2), unpack = True, skiprows = 3)
        plt.plot(x2/10, y2, label = "%s" % (label2), color = "darkgreen")
    plt.fill_between(x1, ysquare1, -99, color = 'blue', alpha = .5)
    plt.ylim(-10, 110)
    plt.xlim(0, x1[-1])
    plt.ylabel("Amplitude [a.u.]")
    plt.xlabel("Time/seconds")
    plt.legend(loc = "upper right")
    plt.grid()
    plt.show()

def get_mode(t, n = 0):
    return (200*(n%2)/np.pi/n)*np.sin(2*np.pi*t*n/period)

def plot_fundamental(n = 0):
    x = []
    y = []
    for i in range(0, len(xarray)):
        x.append(xarray[i])
        y.append(50 + get_mode(xarray[i], n))
    plt.figure()
    plt.plot(xarray, yarray, color = "red", label = "Inner measurement")
    plt.plot(x, y, color = "blue", label = "Signal harmonic: %s" % (n))
    plt.grid()
    plt.legend(loc = "upper right")
    plt.xlabel("Time/seconds")
    plt.ylabel("Amplitude [a.u.]")
    plt.title("Harmonic %s with measurement" % (n))
    plt.show


thermal_diffusivity = 0
def get_thermal_diffusivity(material): # get data from the internet
    data_entry = [["PTFE", 0.124],     # J. Blumm; A. Lindemann; M. Meyer; C. Strasser (2011). "Characterization of PTFE Using Advanced Thermal Analysis Technique". International Journal of Thermophysics.
                  ["PVC", 0.008],      # Jim Wilson (August 2007). "Materials Data"
                  ["Water", 0.143],    # J. Blumm; A. Lindemann (2003–2007). "Characterization of the thermophysical properties of molten polymers and liquids using the flash technique"
                  ["Glass", 0.34],
                  ["Steel", 18.8],     # Lienhard, John H. Lienhard, John H. (2019). A Heat Transfer Textbook (5th ed.)
                  ["Air", 19]          # Jim Wilson (August 2007). "Materials Data"
                  ]
    for i in range (0, len(data_entry)):
        if material == data_entry[i][0]:
            global thermal_diffusivity
            thermal_diffusivity = data_entry[i][1]
            return data_entry[i][1]
    else:
        raise Exception("Could not extract data from list")
        os.exit()

def sine(x, omega, A, b, c):
    out = A*np.sin(omega*x - b) + c
    return out

pcov = []
popt = []
def cal_amplitude(x = [], y = []):
    if x == []:
        x = xarray
    if y == []:
        y = yarray
    global popt, pcov
    popt, pcov = curve_fit(sine, x, y, [2*np.pi/2400, 10, 2*np.pi/1200, 50])
    print(popt)
    print(pcov)
    return popt[1], np.sqrt(pcov[1, 1])

transmission_factor = []
def calc_transmission_factor(omega = 0, delta_r = 0, amplitude_inner = 0, amplitude_outer = 0, err_inner = 0, err_outer = 0):
    global transmission_factor
    if thermal_diffusivity != 0:
        if omega != 0 and delta_r != 0:
            D = thermal_diffusivity
            gamma = np.exp(-np.sqrt(omega/(2*D))*delta_r)
            egamma = 0 # fill in later
            transmission_factor = [gamma, egamma]
            return gamma, egamma
        else:
            raise Exception("Omega and delta r not defined.")
            os.exit()
    elif amplitude_inner != 0 and amplitude_outer != 0:
        gamma = amplitude_inner/amplitude_outer
        egamma = np.sqrt(err_inner**2 + err_outer**2)
        transmission_factor = [gamma, egamma]
        return gamma, egamma
    else:
        raise Exception("Could not calculate transmission factor.")
        os.exit()

phase_lag = []
def calc_phase_lag(omega = 0, delta_r = 0, phase_inner = 0, phase_outer = 0, err_inner = 0, err_outer = 0):
    global phase_lag
    if thermal_diffusivity != 0:
        if omega != 0 and delta_r != 0:
            D = thermal_diffusivity
            phi = (np.sqrt(omega/(2*D)))*delta_r
            ephi = 0
            phase_lag = [phi, ephi]
            return phi, ephi
        else:
            raise Exception("Omega and delta r not defined.")
            os.exit()
    elif omega == 0 and delta_r == 0:
        phi = phase_inner - phase_outer
        ephi = np.sqrt(err_inner**2 + err_outer**2)
        phase_lag = [phi, ephi]
        return phi, ephi
    else:
        raise Exception("Could not calculate phase lag.")
        os.exit()
        
def calc_thermal_diffusivity(choice = "", omega = 2*np.pi/period, delta_r = 0, unc_delta_r = 0):
    if transmission_factor != 0 and choice == "TF":
        if omega != 0 and delta_r != 0:
            Dtf = (omega*(delta_r**2)/(2*((np.log(transmission_factor[0]))**2)))
            eDtf = Dtf*np.sqrt(2*((unc_delta_r/delta_r)**2) + 2*((transmission_factor[1]/transmission_factor[0])**2)) #assuming that gamma and omega constant
            return Dtf, eDtf
        else:
            raise Exception("Could not calculate thermal diffusivity from transmission factor.")
            os.exit()
    elif phase_lag != 0 and choice == "PL":
        if omega != 0 and delta_r != 0:
            Dpl = (omega*(delta_r**2)/(2*((phase_lag[0])**2)))
            eDpl = Dpl*np.sqrt(2*(unc_delta_r/delta_r) + 2*((phase_lag[1]/phase_lag[0])**2))
            return Dpl, eDpl
        else:
            raise Exception("Could not calculate thermal diffusivity from phase lag.")
            os.exit()
    else:
        raise Exception("Transmission factor not defined.")
        os.exit()


    

    
#%%
'''Task 2.3 a)'''
plot_graph("thermal_4min_a.txt", 240)

#%%
'''Task 2.3 b)'''
plot_fundamental(1)
Aopt, Asig = cal_amplitude()
calc_transmission_factor(omega = 0, delta_r = 0, amplitude_inner = Aopt, amplitude_outer = 200/np.pi, err_inner = Asig, err_outer = 0)
calc_phase_lag(omega = 0, delta_r = 0, phase_inner = popt[2], phase_outer = 0, err_inner = np.sqrt(pcov[2,2]), err_outer = 0) # if prefer phase lag = pi, change popt[2]
Dtf = calc_thermal_diffusivity(choice = "TF", omega = 2*np.pi/2400, delta_r = 7.62, unc_delta_r = 0.05)
Dpl = calc_thermal_diffusivity(choice = "PL", omega = 2*np.pi/2400, delta_r = 7.62, unc_delta_r = 0.05)
print("Transmission factor (val, err) :", transmission_factor)
print("Phase lag (val, err)           :", phase_lag)
print("Dtf (val, err)                 :", Dtf)
print("Dpl (val, err)                 :", Dpl)


#%%

# x = []
# y = []
# for i in range (0, len(xarray)):
#     y.append(sine(xarray[i], popt[0], popt[1], popt[2], popt[3]))
#     x.append(xarray[i])

# plt.figure()
# plt.plot(x, y)
# plt.show()












