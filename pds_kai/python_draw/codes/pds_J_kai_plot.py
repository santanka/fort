# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 16:18:06 2021

@author: mkroc
"""

import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt(r"/home/satanka/Documents/fort/pds_kai/pds_J_kai/pds_J_kai_5_cvnminpoint/pds_J_kai_5_cvnminpoint_all.csv", delimiter=',', unpack=True)

channel = 7
#1:静電ポテンシャル, 2:数密度, 3:Alfven速度, 4:圧力, 5:ベータ値, 6:Larmor半径&慣性長, 7:平行圧力, 8:垂直圧力

if (channel == 1):
    x = data[0, :]
    x = np.rad2deg(x)
    y = data[3, :]
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    ax1.set_ylabel('Electrostatic Potential [kV]', fontsize=25)
    ax1.plot(x, y/1E3)
    ax1.grid()
    ax1.tick_params(labelsize=25)
    ax2.set_xlabel('S← MLT [degree] →N', fontsize=25)
    ax2.set_ylabel('Electrostatic Potential [V]', fontsize=25)
    ax2.plot(x, y)
    ax2.set_ylim(-15, 50)
    ax2.grid()
    ax2.tick_params(labelsize=25)
    fig.suptitle('Electrostatic Potential', fontsize=25)
    plt.rcParams["font.size"] = 25


if (channel == 2):
    x = data[0, :]
    x = np.rad2deg(x)
    y = data[4:14, :]*1.E-6
    fig = plt.figure()
    plt.rcParams["font.size"] = 25
    ax = fig.add_subplot(111, xlabel='S← MLT [degree] →N', ylabel='Number Density [cm^-3]', yscale='log', ylim=(1.E-2, 1.E5))   
    ax.plot(x, y[0, :]+y[2, :], c='orange', label='H+(Jupiter)', linestyle='solid', linewidth='3')
    ax.plot(x, y[1, :]+y[3, :], c='purple', label='e-(Jupiter)', linestyle='dashdot', linewidth='3')
    ax.plot(x, y[4, :], c='red', label='H+(Io)', linestyle='solid', linewidth='3')
    ax.plot(x, y[5, :], c='deepskyblue', label='O+(Io)', linestyle='dotted', linewidth='3')
    ax.plot(x, y[6, :], c='green', label='S+(Io)', linestyle='dotted', linewidth='3')
    ax.plot(x, y[7, :], c='lime', label='S2+(Io)', linestyle='solid', linewidth='3')
    ax.plot(x, y[8, :], c='blue', label='cold e-(Io)', linestyle='dashdot', linewidth='3')
    ax.plot(x, y[9, :], c='hotpink', label='hot e-(Io)', linestyle='dashdot', linewidth='3')
    ax.tick_params(labelsize=25)
    fig.suptitle('Number Density')
    ax.grid()
    ax.legend()
    

if (channel == 3):
    x = data[0, :]
    x = np.rad2deg(x)
    y = data[19, :]
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    ax1.set_ylabel('Alfven Speed (VA) [/(Light Speed (c))]', fontsize=25)
    ax1.plot(x, y)
    ax1.set_yscale('log')
    ax1.grid()
    ax1.tick_params(labelsize=25)
    ax2.set_xlabel('S← MLT [degree] →N', fontsize=25)
    ax2.set_ylabel('1 - VA/c', fontsize=25)
    ax2.plot(x, 1-y)
    ax2.set_yscale('log')
    ax2.grid()
    ax2.tick_params(labelsize=25)
    fig.suptitle('Alfven Speed', fontsize=25)
    plt.rcParams["font.size"] = 25
    

if (channel == 4):
    x = data[0, :]
    x = np.rad2deg(x)
    y = data[60:69]*1.E9
    fig = plt.figure()
    plt.rcParams["font.size"] = 25
    ax1 = fig.add_subplot(131, title='perpendicular')
    ax1.set_xlabel('S← MLT [degree] →N')
    ax1.set_ylabel('Pressure [nPa]')
    ax1.set_yscale('log')
    ax1.plot(x, y[0, :], c='purple', label='all')
    ax1.plot(x, y[3, :], c='orange', label='ion')
    ax1.plot(x, y[6, :], c='blue', label='electron')
    ax1.grid()
    ax1.legend()
    ax1.tick_params(labelsize=25)
    ax2 = fig.add_subplot(132, title='parallel')
    ax2.set_xlabel('S← MLT [degree] →N')
    ax2.set_ylabel('Pressure [nPa]')
    ax2.set_yscale('log')
    ax2.plot(x, y[1, :], c='purple', label='all')
    ax2.plot(x, y[4, :], c='orange', label='ion')
    ax2.plot(x, y[7, :], c='blue', label='electron')
    ax2.grid()
    ax2.legend()
    ax2.tick_params(labelsize=25)
    ax3 = fig.add_subplot(133, title='all')
    ax3.set_xlabel('S← MLT [degree] →N')
    ax3.set_ylabel('Pressure [nPa]')
    ax3.set_yscale('log')
    ax3.plot(x, y[2, :], c='purple', label='all')
    ax3.plot(x, y[5, :], c='orange', label='ion')
    ax3.plot(x, y[8, :], c='blue', label='electron')
    ax3.grid()
    ax3.legend()
    ax3.tick_params(labelsize=25)
    
    plt.subplots_adjust(wspace=0.4, hspace=0.6)

if (channel == 5):
    x = data[0, :]
    x = np.rad2deg(x)
    y = data[72, :]
    z = data[78, :]
    fig = plt.figure()
    plt.rcParams["font.size"] = 25
    ax = fig.add_subplot(111, xlabel='S← MLT [degree] →N', ylabel='beta', yscale='log')
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
    plt.rcParams["font.size"] = 25
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
    ax.tick_params(labelsize=25)


if (channel == 7):
    x = data[0, :]
    x = np.rad2deg(x)
    y = data[40:50, :]*1.E9
    fig = plt.figure()
    plt.rcParams["font.size"] = 25
    ax = fig.add_subplot(111, xlabel='S← MLT [degree] →N', ylabel='parallel pressure [nPa]', yscale='log', ylim=(1.E-6, 1.E2)) 
    ax.plot(x, y[0, :]+y[2, :], c='orange', label='H+(Jupiter)', linestyle='solid', linewidth='3')
    ax.plot(x, y[1, :]+y[3, :], c='purple', label='e-(Jupiter)', linestyle='dashdot', linewidth='3')
    ax.plot(x, y[4, :], c='red', label='H+(Io)', linestyle='solid', linewidth='3')
    ax.plot(x, y[5, :], c='deepskyblue', label='O+(Io)', linestyle='dotted', linewidth='3')
    ax.plot(x, y[6, :], c='green', label='S+(Io)', linestyle='dotted', linewidth='3')
    ax.plot(x, y[7, :], c='lime', label='S2+(Io)', linestyle='solid', linewidth='3')
    ax.plot(x, y[8, :], c='blue', label='cold e-(Io)', linestyle='dashdot', linewidth='3')
    ax.plot(x, y[9, :], c='hotpink', label='hot e-(Io)', linestyle='dashdot', linewidth='3')
    ax.grid()
    ax.legend()
    fig.suptitle('Parallel Pressure')
    ax.tick_params(labelsize=25)
    plt.subplots_adjust(wspace=0.4, hspace=0.6)


if (channel == 8):
    x = data[0, :]
    x = np.rad2deg(x)
    y = data[30:40, :]*1.E9
    fig = plt.figure()
    plt.rcParams["font.size"] = 25
    ax = fig.add_subplot(111, xlabel='S← MLT [degree] →N', ylabel='perpendicular pressure [nPa]', yscale='log', ylim=(1.E-6, 1.E2)) 
    ax.plot(x, y[0, :]+y[2, :], c='orange', label='H+(Jupiter)', linestyle='solid', linewidth='3')
    ax.plot(x, y[1, :]+y[3, :], c='purple', label='e-(Jupiter)', linestyle='dashdot', linewidth='3')
    ax.plot(x, y[4, :], c='red', label='H+(Io)', linestyle='solid', linewidth='3')
    ax.plot(x, y[5, :], c='deepskyblue', label='O+(Io)', linestyle='dotted', linewidth='3')
    ax.plot(x, y[6, :], c='green', label='S+(Io)', linestyle='dotted', linewidth='3')
    ax.plot(x, y[7, :], c='lime', label='S2+(Io)', linestyle='solid', linewidth='3')
    ax.plot(x, y[8, :], c='blue', label='cold e-(Io)', linestyle='dashdot', linewidth='3')
    ax.plot(x, y[9, :], c='hotpink', label='hot e-(Io)', linestyle='dashdot', linewidth='3')
    ax.grid()
    ax.legend()
    fig.suptitle('Perpendicular Pressure')
    ax.tick_params(labelsize=25)
    plt.subplots_adjust(wspace=0.4, hspace=0.6)

plt.tight_layout()
plt.show()
