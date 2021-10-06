# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 16:18:06 2021

@author: mkroc
"""

import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt(r"/home/satanka/Documents/fort/pds_kai/pds_J_kai/pds_J_kai_6_cvnminpoint/pds_J_kai_6_cvnminpoint_all.csv", delimiter=',', unpack=True)

channel = 8
RJ = 7.1492E+07
LL = 5.84760 #11で使用
#1:静電ポテンシャル, 2:数密度, 3:Alfven速度, 4:圧力, 5:ベータ値, 6:Larmor半径&慣性長, 7:平行圧力, 8:垂直圧力, 9:ベータ値(圧力形式)
#10:Su et al. (2006)との温度比較, 11:概形図

if (channel == 1):
    x = data[0, :]
    x = np.rad2deg(x)
    #x = (data[1, :]-data[1, 0])/7.1492E7
    y = data[3, :]
    plt.rcParams["font.size"] = 40
    plt.rcParams.update({'mathtext.default': 'default', 'mathtext.fontset': 'stix'})
    fig, (ax1, ax2) = plt.subplots(2, 1) #, sharex=True
    fig.suptitle('Electrostatic Potential')
    ax1.set_xlabel('S← MLAT [degree] →N')
    ax1.set_ylabel('[kV]')
    ax1.plot(x, y/1E3, linewidth='4')
    ax1.minorticks_on()
    ax1.grid(which="both")
    #ax1.tick_params(labelsize=25)
    ax2.set_xlabel('S← MLAT [degree] →N')
    ax2.set_ylabel('[V]')
    ax2.plot(x, y, linewidth='4')
    ax2.set_ylim(-20, 5)
    ax2.minorticks_on()
    ax2.grid(which="both")
    #ax2.tick_params(labelsize=25)
    #fig.suptitle('Electrostatic Potential', fontsize=25)
    plt.subplots_adjust(wspace=0.4, hspace=0.6)
    
    


if (channel == 2):
    x = data[0, :]
    x = np.rad2deg(x)
    #x = (data[1, :]-data[1, 0])/7.1492E7
    y = data[4:14, :]*1.E-6
    plt.rcParams["font.size"] = 40
    plt.rcParams.update({'mathtext.default': 'default', 'mathtext.fontset': 'stix'})
    fig = plt.figure()
    ax = fig.add_subplot(111, xlabel='S← MLAT [degree] →N', ylabel='Number Density [$cm^{-3}$]', yscale='log', ylim=(1.E-2, 1.E5))   
    ax.plot(x, y[0, :]+y[2, :], c='dimgrey', label='$H^+$(Jupiter)', linestyle='solid', linewidth='4', alpha=0.7)
    ax.plot(x, y[1, :]+y[3, :], c='blue', label='$e^-$(Jupiter)', linestyle='dotted', linewidth='4')
    ax.plot(x, y[4, :], c='purple', label='$H^+$(Io)', linestyle='dotted', linewidth='4')
    ax.plot(x, y[5, :], c='orange', label='$O^+$(Io)', linestyle='dotted', linewidth='4')
    ax.plot(x, y[6, :], c='green', label='$S^+$(Io)', linestyle='dotted', linewidth='4')
    ax.plot(x, y[7, :], c='lime', label='$S^{2+}$(Io)', linestyle='dotted', linewidth='4')
    ax.plot(x, y[8, :], c='deepskyblue', label='cold $e^-$(Io)', linestyle='dotted', linewidth='4')
    ax.plot(x, y[9, :], c='hotpink', label='hot $e^-$(Io)', linestyle='dotted', linewidth='4')
    #ax.tick_params(labelsize=25)
    #fig.suptitle('Number Density')
    ax.minorticks_on()
    ax.grid(which="both")
    ax.legend()
    

if (channel == 3):
    x = data[0, :]
    x = np.rad2deg(x)
    y = data[19, :]
    plt.rcParams["font.size"] = 40
    plt.rcParams.update({'mathtext.default': 'default', 'mathtext.fontset': 'stix'})
    fig = plt.figure()
    ax = fig.add_subplot(111, xlabel='S← MLAT [degree] →N', ylabel='Alfven Speed [/(Light Speed)]', yscale='log')
    ax.plot(x, y, linewidth='4')
    ax.minorticks_on()
    ax.grid(which="both")
    #fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    #ax1.set_ylabel('Alfven Speed [/(Light Speed)]')
    #ax1.plot(x, y)
    #ax1.set_yscale('log')
    #ax1.grid()
    #ax1.tick_params(labelsize=25)
    #ax2.set_xlabel('S← MLAT [degree] →N')
    #ax2.set_ylabel('1 - VA/c')
    #ax2.plot(x, 1-y)
    #ax2.set_yscale('log')
    #ax2.grid()
    #ax2.tick_params(labelsize=25)
    #fig.suptitle('Alfven Speed', fontsize=25)
    #plt.rcParams["font.size"] = 25
    

if (channel == 4):
    x = data[0, :]
    x = np.rad2deg(x)
    y = data[60:69]*1.E9
    plt.rcParams["font.size"] = 40
    plt.rcParams.update({'mathtext.default': 'default', 'mathtext.fontset': 'stix'})
    fig = plt.figure()
    ax1 = fig.add_subplot(131, title='perpendicular')
    ax1.set_xlabel('S← MLAT [degree] →N')
    ax1.set_ylabel('Pressure [nPa]')
    ax1.set_yscale('log')
    ax1.plot(x, y[0, :], c='purple', label='all', linewidth='4', alpha=0.7)
    ax1.plot(x, y[3, :], c='orange', label='ion', linewidth='4', alpha=0.7)
    ax1.plot(x, y[6, :], c='blue', label='electron', linewidth='4', alpha=0.7)
    ax1.minorticks_on()
    ax1.grid(which="both")
    ax1.legend()
    #ax1.tick_params(labelsize=25)
    ax2 = fig.add_subplot(132, title='parallel')
    ax2.set_xlabel('S← MLAT [degree] →N')
    ax2.set_ylabel('Pressure [nPa]')
    ax2.set_yscale('log')
    ax2.plot(x, y[1, :], c='purple', label='all', linewidth='4', alpha=0.7)
    ax2.plot(x, y[4, :], c='orange', label='ion', linewidth='4', alpha=0.7)
    ax2.plot(x, y[7, :], c='blue', label='electron', linewidth='4', alpha=0.7)
    ax2.minorticks_on()
    ax2.grid(which="both")
    ax2.legend()
    #ax2.tick_params(labelsize=25)
    ax3 = fig.add_subplot(133, title='total')
    ax3.set_xlabel('S← MLAT [degree] →N')
    ax3.set_ylabel('Pressure [nPa]')
    ax3.set_yscale('log')
    ax3.plot(x, y[2, :], c='purple', label='all', linewidth='4', alpha=0.7)
    ax3.plot(x, y[5, :], c='orange', label='ion', linewidth='4', alpha=0.7)
    ax3.plot(x, y[8, :], c='blue', label='electron', linewidth='4', alpha=0.7)
    ax3.minorticks_on()
    ax3.grid(which="both")
    ax3.legend()
    #ax3.tick_params(labelsize=25)
    
    plt.subplots_adjust(wspace=0.4, hspace=0.6)

if (channel == 5):
    x = data[0, :]
    x = np.rad2deg(x)
    y = data[72, :]
    z = data[78, :]
    fig = plt.figure()
    plt.rcParams["font.size"] = 40
    ax = fig.add_subplot(111, xlabel='S← MLAT [degree] →N', ylabel='beta', yscale='log')
    ax.plot(x, y, c='orange', label='beta')
    ax.plot(x, z, c='dimgrey', label='me/mi', linestyle='-.')
    fig.suptitle('beta(ion, perpendicular)')
    ax.grid()
    ax.legend()
    ax.tick_params(labelsize=25)
    

if (channel == 6):
    x = data[0, :]
    x = np.rad2deg(x)
    y = data[79:83, :]
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
    ax = fig.add_subplot(111, xlabel='S← MLAT [degree] →N', ylabel='length [m]', yscale='log')   
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
    ax.tick_params(labelsize=25)


if (channel == 7):
    x = data[0, :]
    x = np.rad2deg(x)
    y = data[40:50, :]*1.E9
    fig = plt.figure()
    plt.rcParams["font.size"] = 50
    ax = fig.add_subplot(111, xlabel='S← MLAT [degree] →N', ylabel='parallel pressure [nPa]', yscale='log', ylim=(1.E-4, 1.E2)) 
    ax.plot(x, y[0, :]+y[2, :], c='dimgrey', label=r'$\rm{H^+}$(Jupiter)', linestyle='solid', linewidth='4', alpha=0.7)
    ax.plot(x, y[1, :]+y[3, :], c='blue', label=r'$\rm{e^-}$(Jupiter)', linestyle='dotted', linewidth='4')
    ax.plot(x, y[4, :], c='purple', label=r'$\rm{H^+}$(Io)', linestyle='dotted', linewidth='4')
    ax.plot(x, y[5, :], c='orange', label=r'$\rm{O^+}$(Io)', linestyle='dotted', linewidth='4')
    ax.plot(x, y[6, :], c='green', label=r'$\rm{S^+}$(Io)', linestyle='dotted', linewidth='4')
    ax.plot(x, y[7, :], c='lime', label=r'$\rm{S^{2+}}$(Io)', linestyle='dotted', linewidth='4')
    ax.plot(x, y[8, :], c='deepskyblue', label=r'cold $\rm{e^-}$(Io)', linestyle='dotted', linewidth='4')
    ax.plot(x, y[9, :], c='hotpink', label=r'hot $\rm{e^-}$(Io)', linestyle='dotted', linewidth='4')
    ax.minorticks_on()
    ax.grid(which="both")
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0, fontsize=30)
    #fig.suptitle('Parallel Pressure')
    #ax.tick_params(labelsize=50)
    plt.subplots_adjust(wspace=0.4, hspace=0.6)


    


if (channel == 8):
    x = data[0, :]
    x = np.rad2deg(x)
    y = data[30:40, :]*1.E9
    fig = plt.figure()
    plt.rcParams["font.size"] = 50
    ax = fig.add_subplot(111, xlabel='S← MLAT [degree] →N', ylabel='perpendicular pressure [nPa]', yscale='log', ylim=(1.E-4, 1.E2)) 
    ax.plot(x, y[0, :]+y[2, :], c='dimgrey', label=r'$\rm{H^+}$(Jupiter)', linestyle='solid', linewidth='4', alpha=0.7)
    ax.plot(x, y[1, :]+y[3, :], c='blue', label=r'$\rm{e^-}$(Jupiter)', linestyle='dotted', linewidth='4')
    ax.plot(x, y[4, :], c='purple', label=r'$\rm{H^+}$(Io)', linestyle='dotted', linewidth='4')
    ax.plot(x, y[5, :], c='orange', label=r'$\rm{O^+}$(Io)', linestyle='dotted', linewidth='4')
    ax.plot(x, y[6, :], c='green', label=r'$\rm{S^+}$(Io)', linestyle='dotted', linewidth='4')
    ax.plot(x, y[7, :], c='lime', label=r'$\rm{S^{2+}}$(Io)', linestyle='dotted', linewidth='4')
    ax.plot(x, y[8, :], c='deepskyblue', label=r'cold $\rm{e^-}$(Io)', linestyle='dotted', linewidth='4')
    ax.plot(x, y[9, :], c='hotpink', label=r'hot $\rm{e^-}$(Io)', linestyle='dotted', linewidth='4')
    ax.minorticks_on()
    ax.grid(which="both")
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0, fontsize=30)
    #fig.suptitle('Parallel Pressure')
    #ax.tick_params(labelsize=50)
    plt.subplots_adjust(wspace=0.4, hspace=0.6)

if (channel == 9):
    x = data[0, :]
    x = np.rad2deg(x)
    y = data[69:78]
    plt.rcParams["font.size"] = 40
    plt.rcParams.update({'mathtext.default': 'default', 'mathtext.fontset': 'stix'})
    fig = plt.figure()
    ax1 = fig.add_subplot(131, title='perpendicular')
    ax1.set_xlabel('S← MLAT [degree] →N')
    ax1.set_ylabel('plasma beta')
    ax1.set_yscale('log')
    ax1.plot(x, y[0, :], c='purple', label='all', linewidth='4', alpha=0.7)
    ax1.plot(x, y[3, :], c='orange', label='ion', linewidth='4', alpha=0.7)
    ax1.plot(x, y[6, :], c='blue', label='electron', linewidth='4', alpha=0.7)
    ax1.minorticks_on()
    ax1.grid(which="both")
    ax1.legend()
    #ax1.tick_params(labelsize=25)
    ax2 = fig.add_subplot(132, title='parallel')
    ax2.set_xlabel('S← MLAT [degree] →N')
    ax2.set_ylabel('plasma beta')
    ax2.set_yscale('log')
    ax2.plot(x, y[1, :], c='purple', label='all', linewidth='4', alpha=0.7)
    ax2.plot(x, y[4, :], c='orange', label='ion', linewidth='4', alpha=0.7)
    ax2.plot(x, y[7, :], c='blue', label='electron', linewidth='4', alpha=0.7)
    ax2.minorticks_on()
    ax2.grid(which="both")
    ax2.legend()
    #ax2.tick_params(labelsize=25)
    ax3 = fig.add_subplot(133, title='total')
    ax3.set_xlabel('S← MLAT [degree] →N')
    ax3.set_ylabel('plasma beta')
    ax3.set_yscale('log')
    ax3.plot(x, y[2, :], c='purple', label='all', linewidth='4', alpha=0.7)
    ax3.plot(x, y[5, :], c='orange', label='ion', linewidth='4', alpha=0.7)
    ax3.plot(x, y[8, :], c='blue', label='electron', linewidth='4', alpha=0.7)
    ax3.minorticks_on()
    ax3.grid(which="both", axis='both')
    ax3.legend()
    #ax3.tick_params(labelsize=25)
    plt.subplots_adjust(wspace=0.4, hspace=0.6)

if (channel == 10):
    x = data[1, :]/7.1492E7
    latitude = data[0, :]
    latitude = np.rad2deg(latitude)
    y = np.zeros(len(x))
    for ii in range(10):
        y += data[4+ii, :]
    Tem = data[62, :] / y[:]
    Tem = Tem / 1.6021766208E-19
    Su = (0.1 + 49.9*np.tanh(abs(x[0])-abs(x)))
    plt.rcParams["font.size"] = 40
    plt.rcParams.update({'mathtext.default': 'default', 'mathtext.fontset': 'stix'})
    fig = plt.figure()
    ax = fig.add_subplot(111, xlabel='S← MLAT [degree] →N', ylabel='temperature [eV]', yscale='log')
    ax.plot(latitude, Tem, c='b', label='Case 1 (PDS)', linewidth='4', alpha=0.5)
    ax.plot(latitude, Su, c='r', label='Su et al. (2006)', linewidth='4', alpha=0.5)
    #fig.suptitle('Comparison of temperature')
    ax.minorticks_on()
    ax.grid(which="both")
    ax.legend()

if (channel == 11):
    rad = np.linspace(0, 2*np.pi, 10000)
    x = np.cos(rad)
    y = np.sin(rad)
    xiono = np.cos(rad) * (1+2500E3/RJ)
    yiono = np.sin(rad) * (1+2500E3/RJ)
    xIo = np.cos(rad) * (1821.6E3/RJ) + 5.89856
    yIo = np.sin(rad) * (1821.6E3/RJ)
    mlat = data[0, :]
    xfac = LL*np.cos(mlat)**3.
    yfac = LL*np.cos(mlat)**2.*np.sin(mlat)
    fig = plt.figure()
    plt.rcParams["font.size"] = 40
    ax = fig.add_subplot(111, xlabel='r [$R_J$]', ylabel='z [$R_J$]')
    ax.scatter(xfac[0], yfac[0], s=300, color='orange', marker='D', zorder=5)
    ax.scatter(xfac[len(mlat)-1], yfac[len(mlat)-1], s=300, color='orange', marker='D', zorder=5)
    ax.scatter(xfac[(len(mlat)-1)//2], yfac[(len(mlat)-1)//2], s=300, color='orange', marker='D', zorder=5)
    ax.plot(x, y, linewidth='4', color='c', zorder=3)
    ax.plot(xIo, yIo, linewidth='4', color='c', zorder=4)
    ax.plot(xiono, yiono, linewidth='4', color='c', zorder=2)
    ax.plot(xfac, yfac, linewidth='4', color='purple', zorder=1)
    ax.set_aspect('equal')
    ax.grid()
    ax.set_axisbelow(True)


plt.tight_layout()
plt.show()
