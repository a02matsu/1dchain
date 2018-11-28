#!/usr/bin/env python3
# coding: UTF-8

from scipy import stats  
import numpy as np  
import matplotlib.pylab as plt
from scipy import optimize
import sys


NBALL="3"
non_linear="0.000"
Tdiff="1"
Tav="10^3"
#ver="0001"
###################
nbin=25
data=2
###################

#FILE="OUTPUT/SVXX_N" + NBALL.zfill(3) +"NL" + non_linear + "t" + Tdiff.zfill(3) + "D" + Tav.zfill(3) + "_" + ver

args = sys.argv
FILE = args[1]
eigenvals = np.loadtxt(FILE)
dist = np.exp(eigenvals[:,data]) - np.exp(eigenvals[:,data+1] )


#FILE2="OUTPUT/SVXX_N" + NBALL.zfill(3) +"NL" + non_linear + "t" + Tdiff.zfill(3) + "D010" + "_" + ver
#eigenvals2 = np.loadtxt(FILE)
#dist2 = np.exp(eigenvals[:,data]) - np.exp(eigenvals[:,data+1] )
#
#FILE3="OUTPUT/SVXX_N" + NBALL.zfill(3) +"NL" + non_linear + "t" + Tdiff.zfill(3) + "D100" + "_" + ver
#eigenvals3 = np.loadtxt(FILE)
#dist3 = np.exp(eigenvals[:,data]) - np.exp(eigenvals[:,data+1] )
#
def fit_func_GOE(x, b):
    return 2 * b * x * np.exp(-b * x**2)

#def fit_func_GUE(x, c):
#    return 4 * c**(3/2) / np.sqrt(np.pi) * x**2 * np.exp(-c * x**2)

#def fit_func_GSE(x, c):
#    return 8 * c**(5/2) / (3*np.sqrt(np.pi)) * x**4 * np.exp(-c * x**2)

## Histram
#data_label="N=" + NBALL + "," +  \
        #"NL=" + non_linear + "," +  \
        #"Tdiff=" + Tdiff + "," +  \
        #"Tav=" + Tav
data_label=FILE
parameter_initial = np.array([1./(dist[0]**2)]) # initial parameters
bin_heights, bin_borders, _ = plt.hist(dist, bins=nbin, label=data_label, normed=True)
bin_centers = bin_borders[:-1] + np.diff(bin_borders) / 2

## fit
popt, pcov = optimize.curve_fit(fit_func_GOE, bin_centers, bin_heights , p0=parameter_initial)
#popt2, pcov2 = optimize.curve_fit(fit_func_GUE, bin_centers, bin_heights , p0=parameter_initial)
x_interval_for_fit = np.linspace(0.0, bin_borders[-1], 10000)
#plt.plot(x_interval_for_fit, fit_func_GOE(x_interval_for_fit, *popt)) #, label='GOE')
#plt.plot(x_interval_for_fit, fit_func_GUE(x_interval_for_fit, *popt2), label='GUE')
#plt.plot(x_interval_for_fit, fit_func_GSE(x_interval_for_fit, *popt3), label='GSE')

plt.legend()
plt.xlabel(  str(data+1) + '-' + str(data) )


plt.show()

