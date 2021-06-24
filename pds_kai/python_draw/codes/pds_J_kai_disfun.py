# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 21:14:50 2021

@author: mkroc
"""

import numpy as np
import matplotlib.pyplot as plt

channel = 2 #0:比較なし, 1:1つと比較, 2:2つと比較
data1 = np.genfromtxt(r"/home/satanka/Documents/fort/pds_kai/pds_J_kai/pds_J_kai_6_cvnminpoint/pds_J_kai_6_cvnminpoint_disfun_s=1_MLAT=-39.072641287423266.csv", delimiter=',', unpack=True)
vperp = data1[2][:]
vpara = data1[1][:]
ff = data1[8][:]

if(channel == 1 or channel == 2):
    data2 = np.genfromtxt(r"/home/satanka/Documents/fort/pds_kai/pds_J_kai/pds_J_kai_6_cvnminpoint/pds_J_kai_6_cvnminpoint_disfun_s=3_MLAT=-39.072641287423266.csv", delimiter=',', unpack=True)
    vperp2 = data2[2][:]
    vpara2 = data2[1][:]
    ff2 = data2[8][:]
    vperp = np.concatenate([vperp, vperp2], axis=0)
    vpara = np.concatenate([vpara, vpara2], axis=0)
    ff = np.concatenate([ff, ff2], axis=0)

if(channel == 2):
    data3 = np.genfromtxt(r"/home/satanka/Documents/fort/pds_kai/pds_J_kai/pds_J_kai_6_cvnminpoint/pds_J_kai_6_cvnminpoint_disfun_s=5_MLAT=-39.072641287423266.csv", delimiter=',', unpack=True)
    vperp3 = data3[2][:]
    vpara3 = data3[1][:]
    ff3 = data3[8][:]
    vperp = np.concatenate([vperp, vperp3], axis=0)
    vpara = np.concatenate([vpara, vpara3], axis=0)
    ff = np.concatenate([ff, ff3], axis=0)

length = len(ff)
maxff = max(ff)
for ii in range(length):
    if(np.log10(ff[ii]) < np.log10(maxff)-20.):
        ff[ii] = np.nan

fig = plt.figure()
ax = fig.add_subplot(111)

cm = plt.cm.get_cmap('turbo')

ax.set_xlabel("vpara [m/s] (+ : S→N, - : N→S)", fontsize=25)
ax.set_ylabel("vperp [m/s]", fontsize=25)
plt.title("distribution function (scale=log10)", fontsize=25)
if(min(np.floor(ff)) != 0.):
    mappable = ax.scatter(vpara, vperp, c=np.log10(ff), vmin=np.floor(min(np.log10(ff))), vmax=np.trunc(max(np.log10(ff))), cmap=cm)
if(min(np.floor(ff)) == 0.):
    mappable = ax.scatter(vpara, vperp, c=np.log10(ff), vmin=np.floor(max(np.log10(ff))-15.), vmax=np.trunc(max(np.log10(ff))), cmap=cm)

cbar = fig.colorbar(mappable, ax=ax)
cbar.ax.tick_params(labelsize=25)

ax.grid()
ax.tick_params(labelsize=25)
ax.xaxis.offsetText.set_fontsize(25)
ax.yaxis.offsetText.set_fontsize(25)
plt.subplots_adjust(wspace=0.4, hspace=0.6)
plt.rcParams["font.size"] = 25


plt.show()
