# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 21:14:50 2021

@author: mkroc
"""

import numpy as np
import matplotlib.pyplot as plt

data1 = np.genfromtxt(r"/home/satanka/Documents/fort/pds_kai/pds_E_kai/pds_E_kai_4_disfun_s=12_MLT=0.46234700760877978.csv", delimiter=',', unpack=True)
vperp = data1[4][:]
vpara = data1[3][:]
ff = data1[8][:]

fig = plt.figure()
ax = fig.add_subplot(111)

cm = plt.cm.get_cmap('rainbow')

ax.set_xlabel("vpara [m/s] (+ : S→N, - : N→S)")
ax.set_ylabel("vperp [m/s]")
plt.title("distribution function (scale=log10)")
if(min(np.floor(ff)) != 0.):
    mappable = ax.scatter(vpara, vperp, c=np.log10(ff), vmin=np.floor(min(np.log10(ff))), vmax=np.trunc(max(np.log10(ff))), cmap=cm)
if(min(np.floor(ff)) == 0.):
    mappable = ax.scatter(vpara, vperp, c=np.log10(ff), vmin=np.floor(max(np.log10(ff))-5.), vmax=np.trunc(max(np.log10(ff))), cmap=cm)

fig.colorbar(mappable, ax=ax)

ax.grid()

plt.subplots_adjust(wspace=0.4, hspace=0.6)
plt.rcParams["font.size"] = 25


plt.show()
