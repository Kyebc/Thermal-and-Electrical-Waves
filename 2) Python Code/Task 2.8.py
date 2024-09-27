# -*- coding: utf-8 -*-
"""
Created on Sun Feb  5 13:43:27 2023

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
#%%

fourAperiod = [240/1, 240/3, 240/5, 240/7]
fourAtf = [1.0145548805372494e-07, 9.34440478686844e-08, 9.110924575817298e-08, 7.842674461695625e-08]
fourApl = [8.456367349156566e-08, 7.4908130085491e-08, 3.82336698182747e-06, 1.0905471282460851e-06]

oneAperiod = [60/1, 60/3, 60/5, 60/7]
oneAtf = [1.3384209390559388e-07, 1.1275545855025833e-07, 2.4978154438480627e-07, 3.2341207242314596e-07]
oneApl = [1.0851826190662406e-07, 7.06236346085454e-06, 1.641237038196702e-05, 3.37757758782897e-06]

twoAperiod = [120/1, 120/3, 120/5, 120/7]
twoAtf = [9.405769131678028e-08, 9.701212571923659e-08, 1.4600143833724374e-07, 2.057775204191902e-07]
twoApl = [7.68696272335453e-08, 8.079839318609473e-06, 1.9559670598105451e-07, 2.7698863855057637e-07]

sixperiod = [360/1, 360/3, 360/5, 360/7]
sixtf = [1.0564861293192075e-07, 9.769209545227503e-08, 9.604788375147718e-08, 1.0190935996382789e-07]
sixpl = [8.715911486880231e-08, 7.220194316453117e-08, 6.954185782288691e-08, 4.880345766206419e-06]

eightperiod = [480/1, 480/3, 480/5, 480/7]
eighttf = [1.2485238102942335e-07, 1.2525283065566156e-07, 1.0935922771427224e-07, 1.0163185517197822e-07]
eightpl = [1.1130492238690211e-07, 9.153951656981335e-08, 8.383584731671468e-08, 7.419449553923085e-08]

sixteenperiod = [960/1, 960/3, 960/5, 960/7]
sixteentf = [9.187947925671679e-08, 7.14134495682934e-08, 6.096990089333795e-08, 1.1530584883813208e-07]
sixteenpl = [1.1144488043444643e-07, 9.823100976585516e-08, 4.297716876049622e-08, 4.0772352299299524e-08]

#%%

plt.figure()
plt.grid()
plt.scatter(fourAperiod, fourAtf, label = "$D_{TF}$")
plt.scatter(fourAperiod, fourApl, label = "$D_{PL}$")
plt.axhline(1.24e-7, linestyle = "dashed", label = "Actual Diffusivity")
plt.title("Diffusivity for Dataset: Four A")
plt.xlabel("Period/ seconds")
plt.ylabel("Calculated Diffusivity/ $m^{2} s^{-1}$")
plt.legend()
plt.show

plt.figure()
plt.grid()
plt.scatter(oneAperiod, oneAtf, label = "$D_{TF}$")
plt.scatter(oneAperiod, oneApl, label = "$D_{PL}$")
plt.axhline(1.24e-7, linestyle = "dashed", label = "Actual Diffusivity")
plt.title("Diffusivity for Dataset: One A")
plt.xlabel("Period/ seconds")
plt.ylabel("Calculated Diffusivity/ $m^{2} s^{-1}$")
plt.legend()
plt.show

plt.figure()
plt.grid()
plt.scatter(twoAperiod, twoAtf, label = "$D_{TF}$")
plt.scatter(twoAperiod, twoApl, label = "$D_{PL}$")
plt.axhline(1.24e-7, linestyle = "dashed", label = "Actual Diffusivity")
plt.title("Diffusivity for Dataset: Two A")
plt.xlabel("Period/ seconds")
plt.ylabel("Calculated Diffusivity/ $m^{2} s^{-1}$")
plt.legend()
plt.show

plt.figure()
plt.grid()
plt.scatter(sixperiod, sixtf, label = "$D_{TF}$")
plt.scatter(sixperiod, sixpl, label = "$D_{PL}$")
plt.axhline(1.24e-7, linestyle = "dashed", label = "Actual Diffusivity")
plt.title("Diffusivity for Dataset: Six")
plt.xlabel("Period/ seconds")
plt.ylabel("Calculated Diffusivity/ $m^{2} s^{-1}$")
plt.legend()
plt.show

plt.figure()
plt.grid()
plt.scatter(eightperiod, eighttf, label = "$D_{TF}$")
plt.scatter(eightperiod, eightpl, label = "$D_{PL}$")
plt.axhline(1.24e-7, linestyle = "dashed", label = "Actual Diffusivity")
plt.title("Diffusivity for Dataset: Eight")
plt.xlabel("Period/ seconds")
plt.ylabel("Calculated Diffusivity/ $m^{2} s^{-1}$")
plt.legend()
plt.show


plt.figure()
plt.grid()
plt.scatter(sixteenperiod, sixteentf, label = "$D_{TF}$")
plt.scatter(sixteenperiod, sixteenpl, label = "$D_{PL}$")
plt.axhline(1.24e-7, linestyle = "dashed", label = "Actual Diffusivity")
plt.title("Diffusivity for Dataset: Sixteen")
plt.xlabel("Period/ seconds")
plt.ylabel("Calculated Diffusivity/ $m^{2} s^{-1}$")
plt.legend()
plt.show

#%%
`
plt.figure()
plt.grid()
plt.scatter



