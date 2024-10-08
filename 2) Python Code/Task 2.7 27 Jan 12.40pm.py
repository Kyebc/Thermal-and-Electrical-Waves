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
import math
#%%
def sine(x, A, omega, b, c):
    out = A*np.sin(omega*x - b) + c
    return out

#%%
class Thermal():
    
    def __init__(self,
                 signal_amplitude = 0,
                 period = 0,
                 thermal_diffusivity = 0
                 ):
        self.__sigamp = signal_amplitude
        self.__xarray = [[], []]
        self.__yarray = [[], []]
        self.__period = period
        self.__diffty = 0
        self.__trsnmf = list()
        self.__phslag = list()
        self.__dffytf = list()
        self.__dffypl = list()
        self.__sinopt = list()
        self.__sincov = list()
        self.__sigarr = list()
        self.__fnames = ["file 1 not loaded", "file 2 not loaded"]
        self.__sigmod = 0
        self.__deltar = (7.6/1000, 0.1/1000)
    
    def __repr__(self):
        return ("\n%s \nSignal Amplitude     = %g\nSignal Period        = %g\nThermal Diffusivity  = %g\nTransmission Factor  = %s\nPhase Lag            = %s\nDelta r              = %s\n\nCalculated Thermal Diffusivity: \n-transmission factor = %s \n-phase lag           = %s \n\nFile names           : %s, %s\nX Array Values: \n%s\n%s \nY Array Values: \n%s\n%s \nEND OF DOCUMENTATION\n "\
                % ("Details: ", self.__sigamp, self.__period, self.__diffty, self.__trsnmf, self.__phslag, self.__deltar, self.__dffytf, self.__dffypl, self.__fnames[0], self.__fnames[1], self.__xarray[0], self.__xarray[1], self.__yarray[0], self.__yarray[1]))
    
    def __str__(self):      
        return ("\n%s\n%g\n%g\n%g\n%s\n%s\n%s\n%s\n%s\n%s,%s\n%s\n%s\n%s\n%s\n"\
                % ("DETAILS ", self.__sigamp, self.__period, self.__diffty, self.__trsnmf, self.__phslag, self.__deltar, self.__dffytf, self.__dffypl, self.__fnames[0], self.__fnames[1], self.__xarray[0], self.__xarray[1], self.__yarray[0], self.__yarray[1]))
    
    def return_full_documentation(self):
        print(self.__repr__)
        print(self.get_fit_optimal())
        print(self.get_fit_uncertainty())
        pass
    
    def get_signal_amplitude(self):
        return self.__sigamp
    
    def set_signal_amplitude(self, val):
        self.__sigamp = val
        pass
    
    def get_xarray(self, index = 0):
        return self.__xarray[index]
    
    def set_xarray(self, array, index = 0):
        self.__xarray[index] = array
        pass
    
    def get_yarray(self, index = 0):
        return self.__yarray[index]
    
    def set_yarray(self, array, index = 0):
        self.__yarray[index] = array
        pass
    
    def get_period(self):
        return self.__period
    
    def set_period(self, val):
        self.__period = val
        pass
    
    def get_thermal_diffusivity(self):
        return self.__diffty
    
    def set_thermal_diffusivity(self, val):
        self.__diffty = val
        pass
    
    def get_transmission_factor(self):
        return self.__trsnmf
    
    def set_transmission_factor(self, val, err):
        self.__trsnmf = val, err
        pass
    
    def get_phase_lag(self):
        return self.__phslag
    
    def set_phase_lag(self, val, err):
        self.__phslag = val, err
        pass
    
    def get_diffusivity_tf(self):
        return self.__dffytf
    
    def set_diffusivity_tf(self, val, err):
        self.__dffytf = val, err
        pass
    
    def get_diffusivity_pl(self):
        return self.__dffypl
    
    def set_diffusivity_pl(self, val, err):
        self.__dffypl = val, err
        pass
    
    def get_fit_optimal(self, parameter = "None"):
        if parameter == "None":
            print("Optimal fit for \"A*sin(wx + b) + c\": \nA = %g \nw = %g \nb = %g \nc = %g\n " % (self.__sinopt[0], self.__sinopt[1], self.__sinopt[2], self.__sinopt[3]))
        elif parameter == "A":
            return self.__sinopt[0]
        elif parameter == "w":
            return self.__sinopt[1]
        elif parameter == "b":
            return self.__sinopt[2]
        elif parameter == "c":
            return self.__sinopt[3]
        else:
            raise Exception("Could not retrieve optimal fit parameters.")
            os.exit()
        pass
    
    def set_fit_optimal(self, val, parameter = "None"):
        if parameter == "None":
            self.__sinopt    = val
        elif parameter == "A":
            self.__sinopt[0] = val
        elif parameter == "w":
            self.__sinopt[1] = val
        elif parameter == "b":
            self.__sinopt[2] = val
        elif parameter == "c":
            self.__sinopt[3] = val
        else:
            raise Exception("Error setting optimal fit parameters.")
        pass
    
    def get_fit_uncertainty(self, parameter = "None"):
        A = np.sqrt(self.__sincov[0][0])
        w = np.sqrt(self.__sincov[1][1])
        b = np.sqrt(self.__sincov[2][2])
        c = np.sqrt(self.__sincov[3][3])
        if parameter == "None":
            print("The uncertainties are: \nA = %s\nw = %s\nb = %s\nc = %s\n " %(A, w, b, c))
        elif parameter == "A":
            return A
        elif parameter == "w":
            return w
        elif parameter == "b":
            return b
        elif parameter == "c":
            return c
        else:
            raise Exception("Error setting optimal fit parameters.")
        pass
        
    def set_fit_uncertainty(self, val):
        self.__sincov = val
        pass
    
    def get_file_name(self, index = 0):
        return self.__fnames[index]
    
    def set_file_name(self, name, index = 0):
        self.__fnames[index] = name
        pass
    
    def get_mode_number(self):
        return self.__sigmod
    
    def set_mode_number(self, val):
        self.__sigmod = val
        pass
    
    def get_delta_r(self):
        return self.__deltar
        
    def set_delta_r(self, val, err):
        self.__deltar = val, err
        pass
    
    def load_data(self, file_name, period = 0, dataset = 0):
        x, y = np.loadtxt("%s" % (file_name), unpack = True, skiprows = 3)
        if dataset > 1 or period == 0:
            raise Exception("Could not load data. Make sure period is define and dataset is less than '1'.")
            os.exit()
        if period != self.get_period() and self.get_period() != 0:
            raise Exception("Could not load data of two different periods.")
            os.exit()
        self.set_xarray(x/10, dataset)
        self.set_yarray(y, dataset)
        self.set_period(period)
        self.set_file_name(file_name, dataset)
        pass
    
    def load_thermal_diffusivity(self, material):
        data_entry = [["PTFE",  0.124*(10**-6)], # J. Blumm; A. Lindemann; M. Meyer; C. Strasser (2011). "Characterization of PTFE Using Advanced Thermal Analysis Technique". International Journal of Thermophysics.
                      ["PVC",   0.008*(10**-6)], # Jim Wilson (August 2007). "Materials Data"
                      ["Water", 0.143*(10**-6)], # J. Blumm; A. Lindemann (2003–2007). "Characterization of the thermophysical properties of molten polymers and liquids using the flash technique"
                      ["Glass",  0.34*(10**-6)],
                      ["Steel",  18.8*(10**-6)], # Lienhard, John H. Lienhard, John H. (2019). A Heat Transfer Textbook (5th ed.)
                      ["Air",      19*(10**-6)]  # Jim Wilson (August 2007). "Materials Data"
                      ]
        for i in range (0, len(data_entry)):
            if material == data_entry[i][0]:
                self.set_thermal_diffusivity(data_entry[i][1])
        return self.get_thermal_diffusivity()
    
    def get_measurement_amplitude(self):
        val = (max(self.get_yarray()) - min(self.get_yarray()))/2
        return val
    
    def find_average(self):
        val = sum(self.get_yarray())/len(self.get_yarray())
        return val
    
    def fit_sine(self, est = []):
        if est == []:
            inicond = [self.get_measurement_amplitude(), 2*np.pi/self.get_period(), np.pi, self.find_average()]
        else:
            inicond = est
        popt, pcov = curve_fit(sine, self.get_xarray(), self.get_yarray(), p0 = inicond, bounds = (0, [self.get_signal_amplitude(), 1, 2*np.pi, self.get_signal_amplitude()*2]))
        self.set_fit_optimal(popt)
        self.set_fit_uncertainty(pcov)
        pass
    
    def generate_square_signal(self):
        xarray = (self.get_xarray()*2*np.pi)/(self.get_period())
        ysignal = self.get_signal_amplitude()*(1 + sig.square(xarray))
        self.__sigarr = ysignal
        return ysignal
    
    def plot_graph(self, data2 = False, ylower = - 10, yupper = 110):
        plt.figure()
        if data2 == True:
            plt.title("Signal Plots")
            label1 = "%s" % (self.get_file_name(0))
            label2 = "%s" % (self.get_file_name(1))
        elif data2 == False:
            plt.title("Plot: \"%s\"" % (self.get_file_name(0)))
            label1 = "Measurement"
        plt.plot(self.get_xarray(), self.generate_square_signal(), color = "blue", label = "Signal")
        plt.plot(self.get_xarray(), self.get_yarray(), label = "%s" % (label1), color = "red")
        if data2 == True:
            plt.plot(self.get_xarray(1), self.get_yarray(1), label = "%s" % (label2), color = "darkgreen")
        plt.fill_between(self.get_xarray(), self.generate_square_signal(), -999, color = 'blue', alpha = .5)
        plt.ylim(ylower, yupper)
        plt.xlim(self.get_xarray()[0], self.get_xarray()[-1])
        plt.ylabel("Amplitude [a.u.]")
        plt.xlabel("Time/seconds")
        plt.legend(loc = "upper right")
        plt.grid()
        plt.show()
        pass
    
    def initialization(self, file_name, period, signal_amplitude, plot = False):
        self.load_data(file_name, period)
        self.set_signal_amplitude(signal_amplitude)
        if plot == True:
            self.plot_graph()
        pass
    
    def get_mode(self, t):
        n = self.get_mode_number()
        return (200*(n%2)/(np.pi*n))*np.sin(2*np.pi*t*n/self.get_period())
    
    def plot_mode(self, n = 0, plotmode = True):
        self.set_mode_number(n)
        if plotmode == True:            
            x = []
            y = []
            y2 = []
            for i in range(0, len(self.get_xarray())):
                x.append(self.get_xarray()[i])
                y.append(self.get_signal_amplitude() + self.get_mode(self.get_xarray()[i]))
                y2.append(sine(self.get_xarray()[i], self.get_fit_optimal(parameter = "A"), self.get_fit_optimal(parameter = "w"), self.get_fit_optimal(parameter = "b"), self.get_fit_optimal(parameter = "c")))
            plt.figure()
            plt.plot(self.get_xarray(), self.get_yarray(), color = "red", label = "Inner measurement")
            plt.plot(x, y, color = "blue", label = "Signal harmonic: %s" % (n))
            # plt.plot(x, y2, color = "green", label = "Sine fit", linestyle = "dashed") # to plot comparison, uncomment to activate
            plt.grid()
            plt.legend(loc = "upper right")
            plt.xlabel("Time/seconds")
            plt.ylabel("Amplitude [a.u.]")
            plt.title("Harmonic (%s) with Measurement: %s" % (n, self.get_file_name()))
            plt.show
        pass    
    
    def calculate_phase_lag(self):
        phi1 = self.get_fit_optimal("b")
        phi2 = 0
        err1 = self.get_fit_uncertainty("b")
        val = phi1 - phi2
        err = err1
        self.set_phase_lag(val, err)
        pass
        
    def calculate_transmission_factor(self):
        n = self.get_mode_number()
        a1 = self.get_fit_optimal("A")
        a2 = (200*(n%2)/(np.pi*n))
        er1 = self.get_fit_uncertainty("A")
        val = a1/a2
        err = val*np.sqrt((er1/a1)**2)
        self.set_transmission_factor(val, err)
        pass
    
    def get_omega(self):
        val = 2*np.pi*self.get_mode_number()/self.get_period()
        return val
        
    def calculate_diffusivity_TF(self):
        val = self.get_omega()*((self.get_delta_r()[0])**2)/(2*((np.log(self.get_transmission_factor()[0]))**2))
        err = val*(np.sqrt(2*(self.get_delta_r()[1]/self.get_delta_r()[0]) + 2*(self.get_transmission_factor()[1]/self.get_transmission_factor()[0])))
        self.set_diffusivity_tf(val, err)
        pass
        
    def calculate_diffusivity_PL(self):
        val = self.get_omega()*((self.get_delta_r()[0])**2)/(2*((self.get_phase_lag()[0])**2))
        err = val*(np.sqrt(2*(self.get_delta_r()[1]/self.get_delta_r()[0]) + 2*(self.get_phase_lag()[1]/self.get_phase_lag()[0])))
        self.set_diffusivity_pl(val, err)
        pass
    
    def calculate_diffusivities(self, mode = 0, out = False, est = [], plotmode = True):
        self.fit_sine(est)
        self.plot_mode(mode, plotmode)
        self.calculate_transmission_factor()
        self.calculate_diffusivity_TF()
        self.calculate_phase_lag()
        self.calculate_diffusivity_PL()
        if out == True:
            print("Calculated values of diffusivities for '%s' are:\nTransmission Factor = %s\nPhase Lag           = %s\n " % (self.get_file_name(), self.get_diffusivity_tf(), self.get_diffusivity_pl()))
        pass
    
    def calculate_a_n(self, measurement_mode = 1):
        integral = 0
        for i in range(0, len(self.get_xarray())):
            integral += (self.get_yarray()[i] * np.cos(2*np.pi*measurement_mode*self.get_xarray()[i]/self.get_period()))
        norm = (2/len(self.get_xarray()))
        val = norm * integral
        return val
    
    def calculate_b_n(self, measurement_mode = 1):
        integral = 0
        for i in range(0, len(self.get_xarray())):
            integral += (self.get_yarray()[i] * np.sin(2*np.pi*measurement_mode*self.get_xarray()[i]/self.get_period())) #correct the term in phi
        norm = (2/len(self.get_xarray()))
        val = norm * integral
        return val
    
    def calculate_beta_n(self, n = 1):
        a = self.calculate_a_n(n)
        b = self.calculate_b_n(n)
        beta = np.sqrt(a**2 + b**2)
        return beta
    
    def calculate_phase_n(self, n = 1):
        a = self.calculate_a_n(n)
        b = self.calculate_b_n(n)
        if b != 0:
            # val = -1*np.arctan(a/b)
            val = -1*math.atan2(a, b)
            if val < 0:
                val = 2*np.pi + val
            elif val > 2*np.pi:
                val = val - 2*np.pi
        elif a == 0 or b == 0:
            val = 0
       
        return val
    
    def get_reduced_amplitude_phase(self, n = 1):
        val = []
        beta = self.calculate_beta_n(n)
        dphi = self.calculate_phase_n(n)
        val.append(beta)
        val.append(dphi)
        return val
    
    def calculate_fourier_diffusivity(self, n = 1, option = "tf"):
        beta, dphi = self.get_reduced_amplitude_phase(n)[0], self.get_reduced_amplitude_phase(n)[1]
        omega = 2*np.pi*n/self.get_period()
        if option == "tf":
            amplitude_inner = beta
            amplitude_outer = (200*(n%2)/(np.pi*n))
            if amplitude_outer != 0:
                transmission_factor = amplitude_inner/amplitude_outer
                val = omega*((self.get_delta_r()[0])**2)/(2*((np.log(transmission_factor))**2))
                err = val*(np.sqrt(2*(self.get_delta_r()[1]/self.get_delta_r()[0])))
                return val, err
        elif option == "pl":
            if n%2 != 0:
                phase_inner = dphi
                phase_outer = 0
                phase_lag = phase_inner - phase_outer
                val = omega*((self.get_delta_r()[0])**2)/(2*(phase_lag**2))
                err = val*(np.sqrt(2*((self.get_delta_r()[1]/self.get_delta_r()[0])**2)))
                return val, err
        else:
            raise Exception("Invalid option.")
            os.exit()
            pass
        pass
    
    # def plot_fourier(self, n = 1):
    #     x = []
    #     y = []
    #     for i in range(0, (self.get_period()*10)+1):
    #         pool = 0
    #         for j in range(0, n):
    #             amplitude, phase = self.get_reduced_amplitude_phase(j)
    #             omega = 2*np.pi*j/self.get_period()
    #             pool += sine(self.get_xarray()[i], amplitude, omega, phase, 0)
    #             pass
    #         x.append(self.get_xarray()[i])
    #         y.append(pool)
    #         print(i, pool) #
    #         pass
    #     plt.figure()
    #     plt.plot(self.get_xarray()[:(self.get_period()*10)], self.get_yarray()[:(self.get_period()*10)], label = "Measurement", color = "blue")
    #     plt.plot(x, y, label = "Fourier fit to n = %s" % (n), color = "orange", linestyle = "dashed")
    #     plt.show()
    #     pass
    
    def plot_fourier(self, n = 1):
        length = len(self.get_xarray())
        x = self.get_xarray()
        y = np.zeros(length)
        for i in range(0, n + 1):
            omega = 2*np.pi*i/self.get_period()
            amplitude, phase = self.get_reduced_amplitude_phase(i)
            for j in range(0, length):
                y[j] += sine(self.get_xarray()[j], amplitude, omega, phase, 0)
                pass
            pass
        y += self.get_reduced_amplitude_phase(0)[0]/2
        plt.figure()
        plt.title("Fourier fit comparison to measured signal")
        plt.xlabel("Time/seconds")
        plt.ylabel("Signal/[a.u.]")
        plt.plot(self.get_xarray(), self.get_yarray(), label = "Measurement", color = "blue")
        plt.plot(x, y, label = "Fourier fit to n = %s" % (n), color = "orange", linestyle = "dashed")
        plt.grid()
        plt.legend(loc = "upper right")
        plt.show()
        pass
    
             
        

#%%

def plot_D(array, option = "tf"):
    x = []
    y = []
    yerr = []
    for i in range(0, len(arraylist)):
        y.append(eval(arraylist[i] + ".get_diffusivity_" + option + "()[0]"))
        yerr.append(eval(arraylist[i] + ".get_diffusivity_" + option + "()[1]"))
        x.append(eval(arraylist[i] + ".get_period()"))
        pass
    plt.figure()
    plt.plot(x, y, color = "blue")
    plt.errorbar(x, y, yerr, ecolor = "orange", capsize = 8)
    plt.title("D_" + option + " against period")
    plt.xlabel("Period/seconds")
    plt.ylabel("Thermal Diffusivity/[a.u.]")
    plt.show()
    pass



#%%
'''Task 2.3 a)'''
fourA = Thermal()
fourA.initialization("thermal_4min_a.txt", 240, 50, plot = False)

'''Optional to plot two data for comparison'''
# fourA.load_data("thermal_4min_b.txt", 240, dataset = 1)
# fourA.plot_graph(data2 = True)

'''Task 2.3 b)'''
# fourA.set_delta_r(7.62/1000, 0.05/1000)
fourA.calculate_diffusivities(1, out = False, plotmode = False)

#%%
'''Task 2.4'''
oneA = Thermal()
oneA.initialization("thermal_1min_a.txt", 60, 50, plot = False)
# strong uniform transient

oneB = Thermal()
oneB.initialization("thermal_1min_b.txt", 60, 50, plot = False)
# smaller waveform

twoA = Thermal()
twoA.initialization("thermal_2min_a.txt", 120, 50, plot = False)
# good data

twoB = Thermal()
twoB.initialization("thermal_2min_b.txt", 120, 50, plot = False)
# strong transient

six = Thermal()
six.initialization("thermal_6min.txt", 360, 50, plot = False)
# triangular waveform

eight = Thermal()
eight.initialization("thermal_8min.txt", 480, 50, plot = False)
# slightly wavelike waveform

sixteen = Thermal()
sixteen.initialization("thermal_16min.txt", 960, 50, plot = False)
# wavelike waveform

#%%
'''Task 2.5 a)'''
oneA.calculate_diffusivities(1, out = False, plotmode = False)
# oneB.calculate_diffusivities(1, out = True, plotmode = True)
eight.calculate_diffusivities(1, out = False, plotmode = False)

'''Task 2.5 b)'''
twoA.calculate_diffusivities(1, out = False, plotmode = False)
# twoB.calculate_diffusivities(1, out = True, plotmode = False)
six.calculate_diffusivities(1, out = False, plotmode = False)
sixteen.calculate_diffusivities(1, out = False, plotmode = False)

'''Task 2.5 c)'''
arraylist = ["oneA", "twoA", "fourA", "six", "eight", "sixteen"]
plot_D(arraylist, "tf")
plot_D(arraylist, "pl")

#%%
'''Task 2.6 a)'''
print(fourA.get_reduced_amplitude_phase(n = 0))
print(fourA.get_reduced_amplitude_phase(n = 1))
print(fourA.get_reduced_amplitude_phase(n = 2))
print(fourA.get_reduced_amplitude_phase(n = 3))

'''Task 2.6 b)'''
for i in range(1, 5):
    print(eval("fourA.calculate_fourier_diffusivity(n = " + str(i) + ", option = \"tf\")"))

for i in range(1, 5):
    print(eval("fourA.calculate_fourier_diffusivity(n = " + str(i) + ", option = \"pl\")"))

for i in range(1, 5):
    print(eval("fourA.get_reduced_amplitude_phase(n = " + str(i) + ")"))

#%%
'''Task 2.7 a)'''
for i in range(1, 5):
    print(eval("oneA.calculate_fourier_diffusivity(n = " + str(i) + ", option = \"tf\")"))

for i in range(1, 5):
    print(eval("oneA.calculate_fourier_diffusivity(n = " + str(i) + ", option = \"pl\")"))

for i in range(1, 5):
    print(eval("twoA.calculate_fourier_diffusivity(n = " + str(i) + ", option = \"tf\")"))

for i in range(1, 5):
    print(eval("twoA.calculate_fourier_diffusivity(n = " + str(i) + ", option = \"pl\")"))

for i in range(1, 5):
    print(eval("six.calculate_fourier_diffusivity(n = " + str(i) + ", option = \"tf\")"))

for i in range(1, 5):
    print(eval("six.calculate_fourier_diffusivity(n = " + str(i) + ", option = \"pl\")"))

for i in range(1, 5):
    print(eval("eight.calculate_fourier_diffusivity(n = " + str(i) + ", option = \"tf\")"))

for i in range(1, 5):
    print(eval("eight.calculate_fourier_diffusivity(n = " + str(i) + ", option = \"pl\")"))

for i in range(1, 5):
    print(eval("sixteen.calculate_fourier_diffusivity(n = " + str(i) + ", option = \"tf\")"))

for i in range(1, 5):
    print(eval("sixteen.calculate_fourier_diffusivity(n = " + str(i) + ", option = \"pl\")"))




