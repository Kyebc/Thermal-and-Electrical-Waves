# -*- coding: utf-8 -*-
"""
Created on Mon Jan  9 14:36:17 2023

@author: kyebchoo
"""


import numpy as np
import scipy as sci
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#%%
array = np.loadtxt("Sine_low.txt", unpack = True, skiprows = 1)

high = len(array[0,])
pool = 0
for i in range (0, high - 1):
    pool += ((array[0, i + 1] - array[0, i]) * (array[1, i]))
out  = round(pool,5)

plt.figure()
plt.bar(x = array[0,], height = array[1,], color = "none", edgecolor = "blue", width = (array[0,0]-array[0,1]))
plt.plot(array[0,], array[1,], color = "red", linestyle = "dashed")
plt.text(x = 220, y = 0.8, s = "Integrated area: %s" % (out))
plt.title("Low Resolution Sine")
plt.xlabel("x [A.U.]")
plt.ylabel("y [A.U.]")
plt.show()
print(out)

#%%
array = np.loadtxt("Semi_low.txt", unpack = True, skiprows = 1) 

high = len(array[0,])
pool = 0
for i in range (0, high - 1):
    pool += ((array[0, i + 1] - array[0, i]) * (array[1, i]))
out  = round(pool,5)

plt.figure()
plt.bar(x = array[0,], height = array[1,], color = "none", edgecolor = "blue", width = (array[0,0]-array[0,1]))
plt.plot(array[0,], array[1,], color = "red", linestyle = "dashed")
plt.text(x = 1.3, y = 1, s = "Integrated area: %s" % (out))
plt.title("Low Resolution Semicircle")
plt.xlabel("x [A.U.]")
plt.ylabel("y [A.U.]")
plt.show()
print(out)

#%%

array = np.loadtxt("Sine_high.txt", unpack = True, skiprows = 1)

high = len(array[0,])
pool = 0
for i in range (0, high - 1):
    pool += ((array[0, i + 1] - array[0, i]) * (array[1, i]))
out  = round(pool,5)

plt.figure()
plt.bar(x = array[0,], height = array[1,], color = "none", edgecolor = "blue", width = (array[0,0]-array[0,1]))
plt.plot(array[0,], array[1,], color = "red", linestyle = "dashed")
plt.text(x = 220, y = 0.8, s = "Integrated area: %s" % (out))
plt.title("High Resolution Sine")
plt.xlabel("x [A.U.]")
plt.ylabel("y [A.U.]")
plt.show()
print(out)

#%%
array = np.loadtxt("Semi_high.txt", unpack = True, skiprows = 1) 

high = len(array[0,])
pool = 0
for i in range (0, high - 1):
    pool += ((array[0, i + 1] - array[0, i]) * (array[1, i]))
out = round(pool, 5)

plt.figure()
plt.bar(x = array[0,], height = array[1,], color = "none", edgecolor = "blue", width = (array[0,0]-array[0,1]))
plt.plot(array[0,], array[1,], color = "red", linestyle = "dashed")
plt.text(x = 1.3, y = 1, s = "Integrated area: %s" % (out))
plt.title("High Resolution Semicircle")
plt.xlabel("x [A.U.]")
plt.ylabel("y [A.U.]")
plt.show()
print(out)

#%%

array = np.loadtxt("Semi_low.txt", unpack = True, skiprows = 1) 
resolution = 100

high = len(array[0,])
pool = 0
width = (array[0, 1] - array[0, 0])
for i in range (0, high - 1):
    pool += (width * ((array[1, i] + array[1, i + 1])/2))
out  = round(pool,5)

x = []
y = []

for i in range (0, high - 1):
    for j in range (0, resolution):
        x.append(array[0,i] + j*width/resolution)
        if j == 0:
            y.append(0)
        else:
            y.append(array[1, i] + j*((array[1, i + 1] - array[1, i])/resolution))

plt.figure()
plt.plot(x, y)
plt.text(x = 1.3, y = 1, s = "Integrated area: %s" % (out))
plt.title("Low Resolution Semicircle (Trapezium)")
plt.xlabel("x [A.U.]")
plt.ylabel("y [A.U.]")
plt.show()
print(out)
