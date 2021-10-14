import numpy as np
import matplotlib.pyplot as plt
from matplotlib import mathtext
mathtext.FontConstantsBase = mathtext.ComputerModernFontConstants
from matplotlib import ticker, cm, colors


#CGS Gauss系を使用

ni = np.logspace(-2, 6, 1000, base=10) #イオン数密度[cm^-3]
b0 = np.logspace(-5, -1, 1000, base=10) #背景磁束密度[G]
cc = 299792458E2 #光速[cm s^-1]
ee = 1.60217662E-20 * cc #電気素量[Fr]
qi = ee #イオン電荷[Fr]
pressure_ratio = 1.E1 #イオン圧力/電子圧力
phi = 6E2 * 1E8 / cc #波の静電ポテンシャル[statV] ↔ 200 [V]

NI, B0 = np.meshgrid(ni, b0)

deltaBpara_over_B0 = 4.*np.pi/(B0**2.) * (1. + 1/pressure_ratio) * qi*NI*phi

#plotにはSI系

plt.rcParams["font.size"] = 40
plt.rcParams.update({'mathtext.default': 'default', 'mathtext.fontset': 'stix'})

plt.xlabel(r'ion number density $n_i$ [$cm^{-3}$]')
plt.ylabel(r'$B_0$ [T]')
plt.xscale('log')
plt.yscale('log')
plt.title(r'$\frac{\delta B_{\parallel}}{B_0} = \frac{4 \pi}{B_0^2} (1 + \frac{p_e}{p_i}) q_i n_i \varphi$')

cont = plt.contour(NI, B0*1E-4, np.log10(deltaBpara_over_B0), linewidths=5, colors=['black'], levels=np.linspace(-6, 9, 16))
cont.clabel(fontsize=40)

plt.pcolormesh(NI, B0*1E-4, np.log10(deltaBpara_over_B0), cmap='turbo')
pp = plt.colorbar()
pp.set_label(r'$log_{10} (\frac{\delta B_{\parallel}}{B_0})$')

plt.scatter(1, 0.4111037E-07, marker='D', color='r', s=200)

plt.tight_layout()
plt.show()
