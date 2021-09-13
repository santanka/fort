# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 21:36:01 2021

@author: mkroc
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import mathtext
mathtext.FontConstantsBase = mathtext.ComputerModernFontConstants
import matplotlib.ticker as ticker


plt.rcParams.update({'mathtext.default': 'default', 'mathtext.fontset': 'stix'})
plt.rcParams['text.usetex'] = True
plt.rcParams["font.size"] = 50

data1 = np.genfromtxt(r"/home/satanka/Documents/fort/pds_kai/pds_J_kai/pds_J_kai_6_cvnminpoint/pds_J_kai_6_cvnminpoint_disfun_s=5_mu=1.1465648570536209E-013.csv", delimiter=',', unpack=True)
ss = data1[0][:]
ss = np.rad2deg(ss)
vperp = data1[2][:]
vpara = data1[1][:]
ff = data1[9][:]

length = len(ff)
maxff = np.nanmax(ff)
for ii in range(length):
    if(np.log10(ff[ii]) < np.log10(maxff)-15.):
        ff[ii] = np.nan
        vperp[ii] = np.nan
        ss[ii] = np.nan
print(maxff)
print(np.nanmin(np.log10(ff)))
print(np.nanmax(np.log10(ff)))

fig = plt.figure()
ax = fig.add_subplot(111)

cm = plt.cm.get_cmap('turbo')

ax.set_xlabel("S← MLAT [degree] →N")
ax.set_ylabel("$v_{\parallel}$ [m/s] (+ : S→N, - : N→S)")
plt.title("distribution function (scale=log10), $\mu$=1.147E-13")
#ax.yaxis.offsetText.set_fontsize(25)

if np.trunc(max(ff))-10. > np.floor(min(ff)):
    a = np.trunc(np.nanmax(np.log10(ff)))-10.
else:
    a = np.floor(np.nanmin(np.log10(ff)))

mappable = ax.scatter(ss, vpara, c=np.log10(ff), vmin=a, vmax=np.trunc(np.nanmax(np.log10(ff))), cmap=cm, s=700, alpha=1)
ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y,_: f'{y*1E-5:.1f}'+ r'$\times 10^{5}$'))
fig.colorbar(mappable, ax=ax)

ax.minorticks_on()
ax.grid(which="both")
ax.set_axisbelow(True)
ax.tick_params()
plt.subplots_adjust(wspace=0.4, hspace=0.6)

plt.tight_layout()

plt.show()

