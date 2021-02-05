# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 21:14:50 2021

@author: mkroc
"""

import numpy as np
import matplotlib.pyplot as plt

data1 = np.genfromtxt(r"C:\Users\mkroc\Desktop\fort-backup\pds_E_kai_2_disfun_s=12_MLT=0.0000000000000000.csv", delimiter=',', unpack=True)
vperp = data1[2][:]
vpara = data1[1][:]
ff = np.log10(data1[8][:])

fig = plt.figure()
ax = fig.add_subplot(111)

cm = plt.cm.get_cmap('rainbow')

ax.set_xlabel("vpara [m/s] (+ : S→N, - : N→S)")
ax.set_ylabel("vperp [m/s]")
plt.title("distribution function (scale=log10)")
mappable = ax.scatter(vpara, vperp, c=ff, vmin=np.floor(min(ff)), vmax=np.trunc(max(ff)), cmap=cm)

fig.colorbar(mappable, ax=ax)

ax.grid()

plt.subplots_adjust(wspace=0.4, hspace=0.6)
plt.rcParams["font.size"] = 25


plt.show()
