# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 16:18:06 2021

@author: mkroc
"""

import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt(r"pds_kai/pds_E_kai/pds_E_kai_5_all.csv", delimiter=',', unpack=True)

channel = 1
#1:静電ポテンシャル, 2:数密度, 3:Alfven速度, 4:圧力, 5:ベータ値, 6:Larmor半径&慣性長

if (channel == 1):
    x = data[0, :]
    x = np.rad2deg(x)
    y = data[3, :]
    fig = plt.figure()
    ax = fig.add_subplot(111, xlabel='S← MLT [degree] →N', ylabel='Electrostatic Potential [V]')    
    ax.plot(x, y)
    fig.suptitle('Electrostatic Potential')
    ax.grid()
    plt.rcParams["font.size"] = 25


if (channel == 2):
    x = data[0, :]
    x = np.rad2deg(x)
    y = data[4:16, :]*1.E-6
    fig = plt.figure()
    ax = fig.add_subplot(111, xlabel='S← MLT [degree] →N', ylabel='Number Density [cm^-3]', yscale='log', ylim=(1.E-2, 1.E6))   
    ax.plot(x, y[0, :]+y[5, :], c='orange', label='H+(ionosphere)')
    ax.plot(x, y[1, :]+y[6, :], c='red', label='He+')
    ax.plot(x, y[2, :]+y[7, :], c='purple', label='N+')
    ax.plot(x, y[3, :]+y[8, :], c='green', label='O+')
    ax.plot(x, y[4, :]+y[9, :], c='blue', label='e-(ionosphere)')
    ax.plot(x, y[10, :], c='goldenrod', label='H+(magnetosphere)')
    ax.plot(x, y[11, :], c='deepskyblue', label='e-(magnetosphere)')
    fig.suptitle('Number Density')
    ax.grid()
    ax.legend()
    plt.rcParams["font.size"] = 25

if (channel == 3):
    x = data[0, :]
    x = np.rad2deg(x)
    y = data[21, :]
    fig = plt.figure()
    ax = fig.add_subplot(111, xlabel='S← MLT [degree] →N', ylabel='Alfven Speed [/(Light Speed)]')
    ax.plot(x, y)
    fig.suptitle('Alfven Speed')
    ax.grid()
    plt.rcParams["font.size"] = 25

if (channel == 4):
    x = data[0, :]
    x = np.rad2deg(x)
    y = data[70:79]*1.E9
    fig = plt.figure()
    ax1 = fig.add_subplot(131, title='perpendicular')
    ax1.set_xlabel('S← MLT [degree] →N')
    ax1.set_ylabel('Pressure [nPa]')
    ax1.set_yscale('log')
    ax1.plot(x, y[0, :], c='purple', label='all')
    ax1.plot(x, y[3, :], c='orange', label='ion')
    ax1.plot(x, y[6, :], c='blue', label='electron')
    ax1.grid()
    ax1.legend()
    ax2 = fig.add_subplot(132, title='parallel')
    ax2.set_xlabel('S← MLT [degree] →N')
    ax2.set_ylabel('Pressure [nPa]')
    ax2.set_yscale('log')
    ax2.plot(x, y[1, :], c='purple', label='all')
    ax2.plot(x, y[4, :], c='orange', label='ion')
    ax2.plot(x, y[7, :], c='blue', label='electron')
    ax2.grid()
    ax2.legend()
    ax3 = fig.add_subplot(133, title='all')
    ax3.set_xlabel('S← MLT [degree] →N')
    ax3.set_ylabel('Pressure [nPa]')
    ax3.set_yscale('log')
    ax3.plot(x, y[2, :], c='purple', label='all')
    ax3.plot(x, y[5, :], c='orange', label='ion')
    ax3.plot(x, y[8, :], c='blue', label='electron')
    ax3.grid()
    ax3.legend()
    plt.rcParams["font.size"] = 25
    plt.subplots_adjust(wspace=0.4, hspace=0.6)

if (channel == 5):
    x = data[0, :]
    x = np.rad2deg(x)
    y = data[82, :]
    z = data[88, :]
    fig = plt.figure()
    ax = fig.add_subplot(111, xlabel='S← MLT [degree] →N', ylabel='beta', yscale='log')
    ax.plot(x, y, c='orange', label='beta')
    ax.plot(x, z, c='dimgrey', label='me/mi', linestyle='-.')
    fig.suptitle('beta(ion, perpendicular)')
    ax.grid()
    ax.legend()
    plt.rcParams["font.size"] = 25

if (channel == 6):
    x = data[0, :]
    x = np.rad2deg(x)
    y = data[89:93, :]
    BB = data[2][:]
    elr1 = np.sqrt(2.*0.1*9.1093837015E-31/1.602176634E-19)/BB
    elr2 = np.sqrt(2.*1.*9.1093837015E-31/1.602176634E-19)/BB
    elr3 = np.sqrt(2.*10.*9.1093837015E-31/1.602176634E-19)/BB
    elr4 = np.sqrt(2.*100.*9.1093837015E-31/1.602176634E-19)/BB
    fig = plt.figure()
    ax = fig.add_subplot(111, xlabel='S← MLT [degree] →N', ylabel='length [m]', yscale='log')   
    ax.plot(x, y[0, :], c='orange', label='ion Larmor radius - KAW')
    #ax.plot(x, y[1, :], c='deepskyblue', label='electron Larmor radius')
    #ax.plot(x, y[2, :], c='red', label='ion inertial legth')
    ax.plot(x, y[3, :], c='blue', label='electron inertial length - IAW')
    ax.plot(x, elr1, c='indigo', label='electron Larmor radius - 100eV')
    ax.plot(x, elr2, c='violet', label='electron Larmor radius - 1keV')
    ax.plot(x, elr3, c='red', label='electron Larmor radius - 10keV')
    ax.plot(x, elr4, c='green', label='electron Larmor radius - 100keV')
    fig.suptitle('Larmor radius & inertial length')
    ax.grid()
    ax.legend()
    plt.rcParams["font.size"] = 25

plt.show()
