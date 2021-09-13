import numpy as np
import matplotlib.pyplot as plt
from matplotlib import mathtext
mathtext.FontConstantsBase = mathtext.ComputerModernFontConstants

data = np.genfromtxt(r'/home/satanka/Documents/fort/DEM/DEM_E_L=4_vsPDS.csv', delimiter=',', unpack=True)

trigger = 0
if (trigger == 1):
    data2 = np.genfromtxt(r'/home/satanka/Documents/fort/pds_kai/pds_E_kai/pds_E_kai_3_all.csv', delimiter=',', unpack=True)


x = data[0, :]
x = np.rad2deg(x)
y = data[4:, :]*1.E-6

fig = plt.figure()
plt.rcParams["font.size"] = 30
plt.rcParams.update({'mathtext.default': 'default', 'mathtext.fontset': 'stix'})
ax = fig.add_subplot(111, xlabel='MLAT [degree]', ylabel='Number Density [$cm^{-3}$]', yscale='log', ylim=(1.E-2, 1.E6)) 

if (trigger == 1):
    x2 = data2[0, 127:]
    x2 = np.rad2deg(x2)
    y2 = data2[4:16, 127:]*1.E-6

    ax.plot(x2, y2[0, :]+y2[5, :], c='dimgrey', label='$H^+$(ionosphere, PDS)', linestyle='solid', linewidth='4', alpha=0.5)
    ax.plot(x2, y2[4, :]+y2[9, :], c='blue', label='$e^-$(ionosphere, PDS)', linestyle='solid', linewidth='4', alpha=0.5)
    ax.plot(x2, y2[1, :]+y2[6, :], c='red', label='$He^+$(ionosphere, PDS)', linestyle='solid', linewidth='4', alpha=0.5)
    ax.plot(x2, y2[2, :]+y2[7, :], c='orange', label='$N^+$(ionosphere, PDS)', linestyle='solid', linewidth='4', alpha=0.5)
    ax.plot(x2, y2[3, :]+y2[8, :], c='green', label='$O^{+}$(ionosphere, PDS)', linestyle='solid', linewidth='4', alpha=0.5)

ax.plot(x, y[2, :], c='dimgrey', label='$H^+$(DEM)', linestyle='dotted', linewidth='4')
ax.plot(x, y[0, :], c='blue', label='$e^-$(DEM)', linestyle='dotted', linewidth='4')
ax.plot(x, y[3, :], c='red', label='$He^+$(DEM)', linestyle='dotted', linewidth='4')
ax.plot(x, y[4, :], c='orange', label='$N^+$(DEM)', linestyle='dotted', linewidth='4')
ax.plot(x, y[1, :], c='green', label='$O^{+}$(DEM)', linestyle='dotted', linewidth='4')

ax.minorticks_on()
ax.grid(which="both")
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0, fontsize=25)

plt.tight_layout()
plt.show()
