# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 21:36:01 2021

@author: mkroc
"""

import numpy as np
import matplotlib.pyplot as plt

data1 = np.genfromtxt(r"C:\Users\mkroc\Desktop\fort-backup\pds_E_kai_2_disfun_s=12_mu=4.0794332027624548E-012.csv", delimiter=',', unpack=True)
ss = data1[0][:]
ss = np.rad2deg(ss)
vperp = data1[2][:]
vpara = data1[1][:]
ff = np.log10(data1[8][:])

fig = plt.figure()
ax = fig.add_subplot(111)

cm = plt.cm.get_cmap('rainbow')

ax.set_xlabel("S← MLT [degree] →N")
ax.set_ylabel("vpara [m/s] (+ : S→N, - : N→S)")
plt.title("distribution function (scale=log10)")

if np.trunc(max(ff))-10. > np.floor(min(ff)):
    a = np.trunc(max(ff))-10.
else:
    a = np.floor(min(ff))

mappable = ax.scatter(ss, vpara, c=ff, vmin=a, vmax=np.trunc(max(ff)), cmap=cm)

fig.colorbar(mappable, ax=ax)

ax.grid()

plt.subplots_adjust(wspace=0.4, hspace=0.6)
plt.rcParams["font.size"] = 25


plt.show()

