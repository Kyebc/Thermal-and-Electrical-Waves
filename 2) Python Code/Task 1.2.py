# -*- coding: utf-8 -*-
"""
Created on Mon Jan  9 09:59:51 2023

@author: kyebchoo
"""

import numpy as np
import scipy as sci
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#%%
'''To plot the square function'''
period = 240 #tau
amplitude = 100
resolution = 1
cycle = 2
xnum = int(period * resolution)*cycle

tau = []
for i in range (0, xnum):
    tau.append(i/resolution)

taut = []
for i in range (0, xnum):
    taut.append(i/resolution/period)

dataset = []
for i in range (0, xnum):
    j = i%(period * resolution)
    if j < period * resolution/2:
        dataset.append(amplitude)
    elif j >= period * resolution/2:
        dataset.append(0)
    else:
        print("Error")

#%%

n0 = []
def funcn0(x):
    return amplitude/2
for i in range (0, xnum):
    n0.append(funcn0(i/resolution))
# plt.figure()
# plt.plot(tau, n0, color = "magenta", linestyle = "dashed", alpha = 0.5, label = "n = 0")
# plt.show()

n1 = []
def funcn1(x):
    return (200/np.pi)*np.sin(2*np.pi*x/period)
for i in range (0, xnum):
    n1.append(funcn1(i/resolution))
# plt.figure()
# plt.plot(tau, n1, color = "green", linestyle = "dashed", alpha = 0.5, label = "n = 1")
# plt.show()

n2 = []
def funcn2(x):
    return 0
for i in range (0, xnum):
    n2.append(funcn2(i/resolution))
# plt.figure()
# plt.plot(tau, n2, color = "yellow", linestyle = "dashed", alpha = 0.5, label = "n = 2")
# plt.show()

n3 = []
def funcn3(x):
    return (200/np.pi/3)*np.sin(2*np.pi*x*3/period)
for i in range (0, xnum):
    n3.append(funcn3(i/resolution))
# plt.figure()
# plt.plot(tau, n3, color = "orange", linestyle = "dashed", alpha = 0.5, label = "n = 3")
# plt.show()

nsum = []
for i in range (0, xnum):
    nsum.append(funcn0(i/resolution) + funcn1(i/resolution) + funcn2(i/resolution) + funcn3(i/resolution))
# plt.figure()
# plt.plot(tau, nsum, color = "blue", linestyle = "dashed", alpha = 0.5, label = "Fourier sum")
# plt.show()

# ntest = []
def funcnn(x, n):
    return (200*(n%2)/np.pi/n)*np.sin(2*np.pi*x*n/period)
# for i in range (0, xnum):
#     ntest.append(funcnn(i,5))
# plt.figure()
# plt.plot(tau, ntest, color = "orange", linestyle = "dashed", alpha = 0.5, label = "n = 6")
# plt.show()
def sumfuncn(n, xnum):
    hold = []
    for i in range (0, xnum):
        pool = 100/2
        for j in range (1, n + 1):
            pool += funcnn(i/resolution, j)
        hold.append(pool)
    return hold


# plt.figure()
# plt.plot(tau, dataset)
# # plt.scatter(tau, dataset, marker = ".", c = "r")
# plt.plot(tau, n0, color = "magenta", linestyle = "dashed", alpha = 0.5, label = "n = 0")
# plt.plot(tau, n1, color = "green", linestyle = "dashed", alpha = 0.5, label = "n = 1")
# plt.plot(tau, n2, color = "yellow", linestyle = "dashed", alpha = 0.5, label = "n = 2")
# plt.plot(tau, n3, color = "orange", linestyle = "dashed", alpha = 0.5, label = "n = 3")
# plt.plot(tau, nsum, color = "blue", linestyle = "dashed", alpha = 0.5, label = "Fourier sum")
# plt.title("Task 1.2: Fourier Series")
# plt.xlabel("Tau/ [a.u.]")
# plt.ylabel("Amplitude/ [a.u.]")
# plt.legend()
# plt.show()

plt.figure()
plt.plot(taut, dataset)
# plt.scatter(taut, dataset, marker = ".", c = "r")
plt.plot(taut, n0, color = "magenta", linestyle = "dashed", alpha = 0.5, label = "n = 0")
plt.plot(taut, n1, color = "green", linestyle = "dashed", alpha = 0.5, label = "n = 1")
plt.plot(taut, n2, color = "yellow", linestyle = "dashed", alpha = 0.5, label = "n = 2")
plt.plot(taut, n3, color = "orange", linestyle = "dashed", alpha = 0.5, label = "n = 3")
plt.plot(taut, nsum, color = "blue", linestyle = "dashed", alpha = 0.5, label = "Fourier sum")
plt.title("Task 1.2: Fourier Series")
plt.xlabel("t/Tau")
plt.ylabel("Amplitude/ [a.u.]")
plt.legend(loc = 1)
plt.show()

n = 3
n2 = 20
n3 = 500
array = sumfuncn(n, xnum)
array2 = sumfuncn(n2, xnum)
array3 = sumfuncn(n3, xnum)
plt.figure()
plt.plot(taut, dataset)
plt.plot(taut, array, color = "green", linestyle = "dashed", alpha = 0.5, label = ("Fourier sum n = %s" % (n)))
plt.plot(taut, array2, color = "magenta", linestyle = "dashed", alpha = 0.5, label = ("Fourier sum n = %s" % (n2)))
plt.plot(taut, array3, color = "orange", linestyle = "dashed", alpha = 0.5, label = ("Fourier sum n = %s" % (n3)))
plt.title("Task 1.2: Fourier Series")
plt.xlabel("t/Tau")
plt.ylabel("Amplitude/ [a.u.]")
plt.legend(loc = 1)
plt.show()