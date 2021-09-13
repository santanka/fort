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

data1 = np.genfromtxt(r"/home/satanka/Documents/fort/pds_kai/pds_E_kai/pds_E_kai_3_disfun_s=1_MLT=0.46234700760877978.csv", delimiter=',', unpack=True)
vperp = data1[2][:]
vpara = data1[1][:]
ff = data1[9][:]

length = len(ff)
maxff = np.nanmax(ff)
for ii in range(length):
    if(np.log10(ff[ii]) < np.log10(maxff)-20.):
        ff[ii] = np.nan
        vperp[ii] = np.nan
        vpara[ii] = np.nan


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

ax.minorticks_on()
ax.grid(which="both")
ax.set_axisbelow(True)

#ax.tick_params(labelsize=25)
#ax.xaxis.offsetText.set_fontsize(25)
#ax.yaxis.offsetText.set_fontsize(25)
plt.subplots_adjust(wspace=0.4, hspace=0.6)



plt.show()