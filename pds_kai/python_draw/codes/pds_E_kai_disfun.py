# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 21:14:50 2021

@author: mkroc
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import mathtext
mathtext.FontConstantsBase = mathtext.ComputerModernFontConstants

plt.rcParams.update({'mathtext.default': 'default', 'mathtext.fontset': 'stix'})
plt.rcParams['text.usetex'] = True
plt.rcParams["font.size"] = 50

channel = 2 #0:比較なし, 1:1つと比較, 2:2つと比較
data1 = np.genfromtxt(r"/home/satanka/Documents/fort/pds_kai/pds_E_kai_L=9/pds_E_kai_L=9_2_disfun_s=11_MLT=0.0000000000000000.csv", delimiter=',', unpack=True)
vperp = data1[2][:]
vpara = data1[1][:]
ff = data1[9][:]

if(channel == 1 or channel == 2):
    data2 = np.genfromtxt(r"/home/satanka/Documents/fort/pds_kai/pds_E_kai_L=9/pds_E_kai_L=9_2_disfun_s=1_MLT=0.0000000000000000.csv", delimiter=',', unpack=True)
    vperp2 = data2[2][:]
    vpara2 = data2[1][:]
    ff2 = data2[9][:]
    vperp = np.concatenate([vperp, vperp2], axis=0)
    vpara = np.concatenate([vpara, vpara2], axis=0)
    ff = np.concatenate([ff, ff2], axis=0)

if(channel == 2):
    data3 = np.genfromtxt(r"/home/satanka/Documents/fort/pds_kai/pds_E_kai_L=9/pds_E_kai_L=9_2_disfun_s=6_MLT=0.0000000000000000.csv", delimiter=',', unpack=True)
    vperp3 = data3[2][:]
    vpara3 = data3[1][:]
    ff3 = data3[9][:]
    vperp = np.concatenate([vperp, vperp3], axis=0)
    vpara = np.concatenate([vpara, vpara3], axis=0)
    ff = np.concatenate([ff, ff3], axis=0)

length = len(ff)
maxff = np.nanmax(ff)
for ii in range(length):
    if(np.log10(ff[ii]) < np.log10(maxff)-30.):
        ff[ii] = np.nan
        vperp[ii] = np.nan
        vpara[ii] = np.nan

print(maxff)
print(np.nanmin(np.log10(ff)))
print(np.nanmax(np.log10(ff)))

fig = plt.figure()
ax = fig.add_subplot(111)

cm = plt.cm.get_cmap('turbo')

ax.set_xlabel("$v_{\parallel}$ [m/s] (+ : S→N, - : N→S)")
ax.set_ylabel("$v_{\perp}$ [m/s]")
plt.title("distribution function (scale=log10)")
if(min(np.floor(ff)) != 0.):
    mappable = ax.scatter(vpara, vperp, c=np.log10(ff), vmin=np.floor(np.nanmin(np.log10(ff))), vmax=np.trunc(np.nanmax(np.log10(ff))), cmap=cm, s=700, alpha=0.7)
if(min(np.floor(ff)) == 0.):
    mappable = ax.scatter(vpara, vperp, c=np.log10(ff), vmin=np.floor(np.nanmax(np.log10(ff))-15.), vmax=np.trunc(np.nanmax(np.log10(ff))), cmap=cm, s=700, alpha=0.7)

cbar = fig.colorbar(mappable, ax=ax)
#cbar.ax.tick_params(labelsize=25)

ax.minorticks_on()
ax.grid(which="both")
ax.set_axisbelow(True)

#ax.tick_params(labelsize=25)
#ax.xaxis.offsetText.set_fontsize(25)
#ax.yaxis.offsetText.set_fontsize(25)
plt.subplots_adjust(wspace=0.4, hspace=0.6)



plt.show()
