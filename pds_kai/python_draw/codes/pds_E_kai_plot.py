# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 16:18:06 2021

@author: mkroc
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import mathtext
mathtext.FontConstantsBase = mathtext.ComputerModernFontConstants

data = np.genfromtxt(r"/home/satanka/Documents/fort/pds_kai/pds_E_kai_L=9/pds_E_kai_L=9_1_all.csv", delimiter=',', unpack=True)

channel = 12
RE = 6371E3
LL = 9 #8, 10で使用

print(data[58, 254]/data[4, 254]/1.602176634E-19)
print(data[62, 254]/data[8, 254]/1.602176634E-19)
print(data[58, 127]/data[4, 127]/1.602176634E-19)
print(data[62, 127]/data[8, 127]/1.602176634E-19)

#1:静電ポテンシャル, 2:数密度, 3:Alfven速度, 4:圧力, 5:ベータ値, 6:Larmor半径&慣性長, 7:平行圧力, 8:Ozhogin et al. model,
#9:plasma beta, 10:概形図, 11:Bwpara vs. B0 (Schekochihin et al., 2009) [phi = 200 V]

if (channel == 1):
    x = data[0, :]
    x = np.rad2deg(x)
    #x = data[1, :]/RE
    y = data[3, :]
    fig = plt.figure()
    plt.rcParams["font.size"] = 40
    ax = fig.add_subplot(111, xlabel='S← MLAT [degree] →N', ylabel='Electrostatic Potential [V]') 
    #ax = fig.add_subplot(111, xlabel='S← z [RE] →N', ylabel='Electrostatic Potential [V]')   
    ax.plot(x, y, linewidth='4')
    #fig.suptitle('Electrostatic Potential')
    ax.minorticks_on()
    ax.grid(which="both")
    


if (channel == 2):
    x = data[0, :]
    x = np.rad2deg(x)
    #x = data[1, :]/RE
    y = data[4:16, :]*1.E-6
    fig = plt.figure()
    plt.rcParams["font.size"] = 40
    plt.rcParams.update({'mathtext.default': 'default', 'mathtext.fontset': 'stix'})
    ax = fig.add_subplot(111, xlabel='S← MLAT [degree] →N', ylabel='Number Density [$cm^{-3}$]', yscale='log', ylim=(1.E-2, 1.E6))
    #ax = fig.add_subplot(111, xlabel='S← z [RE] →N', ylabel='Number Density [$cm^{-3}$]', yscale='log', ylim=(1.E-2, 1.E6))   
    ax.plot(x, y[0, :]+y[5, :], c='dimgrey', label='$H^+$(ionosphere)', linestyle='solid', linewidth='4')
    ax.plot(x, y[4, :]+y[9, :], c='blue', label='$e^-$(ionosphere)', linestyle='dotted', linewidth='4')
    ax.plot(x, y[1, :]+y[6, :], c='red', label='$He^+$(ionosphere)', linestyle='dotted', linewidth='4')
    ax.plot(x, y[2, :]+y[7, :], c='orange', label='$N^+$(ionosphere)', linestyle='dotted', linewidth='4')
    ax.plot(x, y[3, :]+y[8, :], c='green', label='$O^{+}$(ionosphere)', linestyle='dotted', linewidth='4')
    ax.plot(x, y[10, :], c='purple', label='$H^+$(magnetosphere)', linestyle='solid', linewidth='4')
    ax.plot(x, y[11, :], c='deepskyblue', label='$e^-$(magnetosphere)', linestyle='dotted', linewidth='4')
    #fig.suptitle('Number Density')
    ax.minorticks_on()
    ax.grid(which="both")
    ax.legend()
    #ax.legend(bbox_to_anchor=(1, 0), loc='lower right', borderaxespad=0)
    

if (channel == 3):
    x = data[0, :]
    x = np.rad2deg(x)
    #x = data[1, :]/RE
    y = data[21, :]
    fig = plt.figure()
    plt.rcParams["font.size"] = 40
    ax = fig.add_subplot(111, xlabel='S← MLAT [degree] →N', ylabel='Alfven Speed [/(Light Speed)]')
    #ax = fig.add_subplot(111, xlabel='S← z [RE] →N', ylabel='Alfven Speed [/(Light Speed)]', yscale='log')
    ax.plot(x, y, linewidth='4')
    #fig.suptitle('Alfven Speed')
    ax.minorticks_on()
    ax.grid(which="both")
    

if (channel == 4):
    x = data[0, :]
    x = np.rad2deg(x)
    #x = data[1, :]/RE
    y = data[70:79]*1.E9
    fig = plt.figure()
    plt.rcParams["font.size"] = 40
    plt.rcParams.update({'mathtext.default': 'default', 'mathtext.fontset': 'stix'})
    ax1 = fig.add_subplot(131, title='perpendicular')
    ax1.set_xlabel('S← MLAT [degree] →N')
    #ax1.set_xlabel('S← z [RE] →N')
    ax1.set_ylabel('Pressure [nPa]')
    ax1.set_yscale('log')
    ax1.plot(x, y[0, :], c='purple', label='all', linewidth='4', alpha=0.7)
    ax1.plot(x, y[3, :], c='orange', label='ion', linewidth='4', alpha=0.7)
    ax1.plot(x, y[6, :], c='blue', label='electron', linewidth='4', alpha=0.7)
    ax1.minorticks_on()
    ax1.grid(which="both")
    ax1.legend()
    ax2 = fig.add_subplot(132, title='parallel')
    ax2.set_xlabel('S← MLAT [degree] →N')
    #ax2.set_xlabel('S← z [RE] →N')
    ax2.set_ylabel('Pressure [nPa]')
    ax2.set_yscale('log')
    ax2.plot(x, y[1, :], c='purple', label='all', linewidth='4', alpha=0.7)
    ax2.plot(x, y[4, :], c='orange', label='ion', linewidth='4', alpha=0.7)
    ax2.plot(x, y[7, :], c='blue', label='electron', linewidth='4', alpha=0.7)
    ax2.minorticks_on()
    ax2.grid(which="both")
    ax2.legend()
    ax3 = fig.add_subplot(133, title='total')
    ax3.set_xlabel('S← MLAT [degree] →N')
    #ax3.set_xlabel('S← z [RE] →N')
    ax3.set_ylabel('Pressure [nPa]')
    ax3.set_yscale('log')
    ax3.plot(x, y[2, :], c='purple', label='all', linewidth='4', alpha=0.7)
    ax3.plot(x, y[5, :], c='orange', label='ion', linewidth='4', alpha=0.7)
    ax3.plot(x, y[8, :], c='blue', label='electron', linewidth='4', alpha=0.7)
    ax3.minorticks_on()
    ax3.grid(which="both")
    ax3.legend()
    
    plt.subplots_adjust(wspace=0.4, hspace=0.6)

if (channel == 5):
    x = data[0, :]
    x = np.rad2deg(x)
    y = data[82, :]
    z = data[88, :]
    fig = plt.figure()
    plt.rcParams["font.size"] = 40
    ax = fig.add_subplot(111, xlabel='S← MLAT [degree] →N', ylabel='beta', yscale='log')
    ax.plot(x, y, c='orange', label='beta')
    ax.plot(x, z, c='dimgrey', label='me/mi', linestyle='-.')
    fig.suptitle('beta(ion, perpendicular)')
    ax.grid()
    ax.legend()
    

if (channel == 6):
    x = data[0, :]
    x = np.rad2deg(x)
    y = data[89:93, :]
    BB = data[2][:]
    me = 9.1093837015E-31
    ee = 1.602176634E-19
    Ke = np.array([0.1, 1.0, 10., 100.])
    Ke = Ke*ee*1.E3
    cc = 2.99792458E+08
    length = len(BB)
    elr = np.zeros([4, length])
    for i in range(length):
        for j in range(4):
            elr[j, i] = 1./cc/ee/BB[i]*np.sqrt(Ke[j]*(me*cc*cc+Ke[j]))
    fig = plt.figure()
    plt.rcParams["font.size"] = 40
    ax = fig.add_subplot(111, xlabel='S← MLT [degree] →N', ylabel='length [m]', yscale='log')   
    ax.plot(x, y[0, :], c='orange', label='ion Larmor radius - KAW')
    #ax.plot(x, y[1, :], c='deepskyblue', label='electron Larmor radius')
    #ax.plot(x, y[2, :], c='red', label='ion inertial legth')
    ax.plot(x, y[3, :], c='blue', label='electron inertial length - IAW')
    ax.plot(x, elr[0][:], c='indigo', label='electron Larmor radius - 100eV')
    ax.plot(x, elr[1][:], c='violet', label='electron Larmor radius - 1keV')
    ax.plot(x, elr[2][:], c='red', label='electron Larmor radius - 10keV')
    ax.plot(x, elr[3][:], c='green', label='electron Larmor radius - 100keV')
    fig.suptitle('Larmor radius & inertial length')
    ax.grid()
    ax.legend()

if (channel == 7):
    x = data[0, :]
    x = np.rad2deg(x)
    y = data[46:58, :]*1.E9
    fig = plt.figure()
    plt.rcParams["font.size"] = 40
    ax = fig.add_subplot(111, xlabel='S← MLT [degree] →N', ylabel='parallel pressure [nPa]', yscale='log')   
    ax.plot(x, y[0, :]+y[5, :], c='orange', label='H+(ionosphere)', linestyle='solid', linewidth='2')
    ax.plot(x, y[1, :]+y[6, :], c='red', label='He+', linestyle='dashed', linewidth='2')
    ax.plot(x, y[2, :]+y[7, :], c='purple', label='N+', linestyle='dotted', linewidth='4')
    ax.plot(x, y[3, :]+y[8, :], c='green', label='O+', linestyle='dotted', linewidth='2')
    ax.plot(x, y[4, :]+y[9, :], c='blue', label='e-(ionosphere)', linestyle='dashdot', linewidth='2')
    ax.plot(x, y[10, :], c='goldenrod', label='H+(magnetosphere)', linestyle='solid', linewidth='4')
    ax.plot(x, y[11, :], c='deepskyblue', label='e-(magnetosphere)', linestyle='dashdot', linewidth='4')
    ax.grid()
    ax.legend()
    
    plt.subplots_adjust(wspace=0.4, hspace=0.6)


if (channel == 8):
    x = data[0, :]
    y = data[4:16, :]*1.E-6
    z = np.zeros(len(x))
    nabINV = np.arccos(np.sqrt(1./LL))
    Neq = pow(10., 4.4693-0.4903*LL)
    z = Neq*pow(np.cos(np.pi/2.*1.01*x/nabINV), -0.75)
    x = np.rad2deg(x)
    mlatlimit = np.rad2deg(np.arccos(np.sqrt(8371E3/RE/LL)))
    plt.rcParams["font.size"] = 40
    plt.rcParams.update({'mathtext.default': 'default', 'mathtext.fontset': 'stix'})
    fig = plt.figure()
    ax = fig.add_subplot(111, xlabel='S← MLAT [degree] →N', ylabel='Number Density [$cm^{-3}$]', yscale='log', ylim=(1.E0, 3.E4), xlim=(-mlatlimit, mlatlimit))   
    ax.plot(x, y[0, :]+y[5, :], c='dimgrey', label='$H^+$(ionosphere, PDS)', linestyle='solid', linewidth='4')
    ax.plot(x, y[4, :]+y[9, :], c='blue', label='$e^-$(ionosphere, PDS)', linestyle='dotted', linewidth='4')
    ax.plot(x, y[1, :]+y[6, :], c='red', label='$He^+$(ionosphere, PDS)', linestyle='dotted', linewidth='4')
    ax.plot(x, z, c='orange', label='$e^-$(RPI model [Ozhogin et al., 2012])', linestyle='solid', linewidth='4')
    #fig.suptitle('Number Density')
    ax.minorticks_on()
    ax.grid(which="both")
    ax.legend()

if (channel == 9):
    x = data[0, :]
    x = np.rad2deg(x)
    #x = data[1, :]/RE
    y = data[79:88]
    fig = plt.figure()
    plt.rcParams["font.size"] = 40
    plt.rcParams.update({'mathtext.default': 'default', 'mathtext.fontset': 'stix'})
    ax1 = fig.add_subplot(131, title='perpendicular')
    ax1.set_xlabel('S← MLAT [degree] →N')
    #ax1.set_xlabel('S← z [RE] →N')
    ax1.set_ylabel('plasma beta')
    ax1.set_yscale('log')
    ax1.plot(x, y[0, :], c='purple', label='all', linewidth='4', alpha=0.7)
    ax1.plot(x, y[3, :], c='orange', label='ion', linewidth='4', alpha=0.7)
    ax1.plot(x, y[6, :], c='blue', label='electron', linewidth='4', alpha=0.7)
    ax1.minorticks_on()
    ax1.grid(which="both")
    ax1.legend()
    ax2 = fig.add_subplot(132, title='parallel')
    ax2.set_xlabel('S← MLAT [degree] →N')
    #ax2.set_xlabel('S← z [RE] →N')
    ax2.set_ylabel('plasma beta')
    ax2.set_yscale('log')
    ax2.plot(x, y[1, :], c='purple', label='all', linewidth='4', alpha=0.7)
    ax2.plot(x, y[4, :], c='orange', label='ion', linewidth='4', alpha=0.7)
    ax2.plot(x, y[7, :], c='blue', label='electron', linewidth='4', alpha=0.7)
    ax2.minorticks_on()
    ax2.grid(which="both")
    ax2.legend()
    ax3 = fig.add_subplot(133, title='total')
    ax3.set_xlabel('S← MLAT [degree] →N')
    #ax3.set_xlabel('S← z [RE] →N')
    ax3.set_ylabel('plasma beta')
    ax3.set_yscale('log')
    ax3.plot(x, y[2, :], c='purple', label='all', linewidth='4', alpha=0.7)
    ax3.plot(x, y[5, :], c='orange', label='ion', linewidth='4', alpha=0.7)
    ax3.plot(x, y[8, :], c='blue', label='electron', linewidth='4', alpha=0.7)
    ax3.minorticks_on()
    ax3.grid(which="both")
    ax3.legend()


if (channel == 10):
    rad = np.linspace(0, 2*np.pi, 10000)
    x = np.cos(rad)
    y = np.sin(rad)
    xiono = np.cos(rad) * (1+500E3/RE)
    yiono = np.sin(rad) * (1+500E3/RE)
    mlat = data[0, :]
    xfac = LL*np.cos(mlat)**3.
    yfac = LL*np.cos(mlat)**2.*np.sin(mlat)
    fig = plt.figure()
    plt.rcParams["font.size"] = 40
    ax = fig.add_subplot(111, xlabel='r [$R_E$]', ylabel='z [$R_E$]')
    ax.scatter(xfac[0], yfac[0], s=500, color='orange', marker='D', zorder=5)
    ax.scatter(xfac[len(mlat)-1], yfac[len(mlat)-1], s=500, color='orange', marker='D', zorder=5)
    ax.scatter(xfac[(len(mlat)-1)//2], yfac[(len(mlat)-1)//2], s=500, color='orange', marker='D', zorder=5)
    ax.plot(x, y, linewidth='4', color='c', zorder=3)
    ax.plot(xiono, yiono, linewidth='4', color='c', zorder=2)
    ax.plot(xfac, yfac, linewidth='4', color='purple', zorder=1)
    ax.set_aspect('equal')
    ax.grid()
    ax.set_axisbelow(True)

if (channel == 11):
    x = data[0, :]
    x = np.rad2deg(x)
    beta_i = data[84, :]
    pressure_i = data[75, :]-data[68, :] #data[75, :] #data[68, :]
    pressure_e = data[78, :]-data[69, :] #data[78, :] #data[69, :]
    #n_i = data[4, :] + data[5, :] + data[6, :] + data[7, :]
    #n_i = n_i + data[9, :] + data[10, :] + data[11, :] + data[12, :] + data[14, :]
    n_i = data[4, :] + data[5, :] + data[6, :] + data[7, :]
    n_i = n_i + data[9, :] + data[10, :] + data[11, :] + data[12, :]
    #n_i = data[14, :]
    ee = 1.602176634E-19
    delta_B = beta_i/2. * (1 + pressure_e/pressure_i) * n_i*ee/pressure_i * 200. #200 [V]
    fig = plt.figure()
    plt.rcParams["font.size"] = 40
    fig.suptitle(r'$\frac{\delta B_{\parallel}}{B_0} = \frac{\beta_i}{2} (1 + \frac{Z}{\tau}) \frac{Ze \varphi}{T_i}$')
    ax = fig.add_subplot(111, xlabel='S← MLAT [degree] →N', yscale='log')
    ax.plot(x, delta_B, linewidth='4', label=r'$\frac{\delta B_{\parallel}}{B_0}$')
    #ax.plot(x, 1 + pressure_e/pressure_i, linewidth='4', label=r'$1 + \frac{Z}{\tau}$')
    ax.plot(x, beta_i/2. *n_i*ee*200./pressure_i, linewidth='4', label=r'$\frac{\beta_i}{2} \frac{Ze \varphi}{T_i}$')
    ax.plot(x, beta_i, linewidth='4', label=r'$\frac{\beta_i}{2}$')
    ax.minorticks_on()
    ax.grid(which="both")
    ax.legend()

if (channel == 12):
    x = data[0, :]
    x = np.rad2deg(x)
    pressure_i = data[75, :]
    pressure_i_2 = data[68, :]
    pressure_i_3 = data[75, :]-data[68, :]
    pressure_e = data[78, :]
    n_i = data[4, :] + data[5, :] + data[6, :] + data[7, :]
    n_i = n_i + data[9, :] + data[10, :] + data[11, :] + data[12, :] + data[14, :]
    n_i_2 = data[14, :]
    n_i_3 = n_i - n_i_2
    n_e = data[8, :] + data[13, :] + data[15, :]
    ee = 1.602176634E-19
    T_i = pressure_i / n_i /ee
    T_i_2 = pressure_i_2 / n_i_2 /ee
    T_i_3 = pressure_i_3 / n_i_3 /ee
    T_e = pressure_e / n_e /ee
    fig = plt.figure()
    plt.rcParams["font.size"] = 40
    fig.suptitle(r'Temperature')
    ax = fig.add_subplot(111, xlabel='S← MLAT [degree] →N', ylabel='Temperature [eV]', yscale='log')
    ax.plot(x, T_i, linewidth='4', label=r'$T_i$')
    ax.plot(x, T_i_2, linewidth='4', label=r'$T_i2$')
    ax.plot(x, T_i_3, linewidth='4', label=r'$T_i3$')
    ax.plot(x, T_e, linewidth='4', label=r'$T_e$')
    ax.minorticks_on()
    ax.grid(which="both")
    ax.legend()

plt.tight_layout()
plt.show()
