import numpy as np
import matplotlib.pyplot as plt
from matplotlib import mathtext
mathtext.FontConstantsBase = mathtext.ComputerModernFontConstants

data1 = np.genfromtxt(r"/home/satanka/Documents/fort/pds_kai/pds_J_kai/SVC_J/SVC_penum_39_50_75.csv", delimiter=',', unpack=True)

data2 = np.genfromtxt(r"/home/satanka/Documents/fort/pds_kai/pds_J_kai/SVC_J/SVC_IC_J_30kV_75_50.csv", delimiter=',', unpack=True)

data3 = np.genfromtxt(r"/home/satanka/Documents/fort/pds_kai/pds_J_kai/pds_J_kai_initialMatsuda/pds_J_kai_initialMatsuda_cvnminpoint_65-2_all.csv", delimiter=',', unpack=True)

channel = 5
#1:静電ポテンシャル、2:数密度、3:MLATと木星電離圏端からの距離の対応関係、4:pdsとの比較(数密度)、5:pdsとの比較(静電ポテンシャル)

if (channel == 1):
    x = data2[2, 1:]

    size = len(x)
    MLAT = np.zeros(size)
    req = 7.1492E+07*5.84760
    
    for jj in range(size):
        MLAT0 = 1.
        for ii in range(1000000):
            if (ii == 1000000):
                print("Error!: solution is not found. z_position = " + str(x(jj)))
                
            ff = req * ((1. / 2.) * np.sin(MLAT0) * np.sqrt(3. * np.sin(MLAT0)**2. + 1.) \
                + (1. / (2. * np.sqrt(3.))) * np.arcsinh(np.sqrt(3.) * np.sin(MLAT0))) \
                - x[jj]
            gg = req * np.cos(MLAT0) * np.sqrt(3. * np.sin(MLAT0)**2. + 1.)

            MLAT1 = float(MLAT0 - ff / gg)
            
            if (abs(MLAT1 - MLAT0) <= 1E-5):
                break

            MLAT0 = MLAT1
        
        MLAT[jj] = MLAT1
        print(x[jj], MLAT[jj])

    MLAT = np.rad2deg(MLAT)
    y = data1[0, :]
    plt.rcParams["font.size"] = 30
    plt.rcParams.update({'mathtext.default': 'default', 'mathtext.fontset': 'stix'})
    fig, (ax1, ax2) = plt.subplots(2, 1) #, sharex=True
    
    ax1.set_xlabel('Io← MLAT [degree] →Jupiter')
    ax1.set_ylabel('Electrostatic Potential [kV]')
    ax1.plot(MLAT, y/1E3, linewidth='4')
    ax1.minorticks_on()
    ax1.grid(which="both")
    #ax1.tick_params(labelsize=25)
    ax2.set_xlabel('Io← MLAT [degree] →Jupiter')
    ax2.set_ylabel('Electrostatic Potential [V]')
    ax2.plot(MLAT, y, linewidth='4')
    ax2.set_ylim(-20, 5)
    ax2.minorticks_on()
    ax2.grid(which="both")


if (channel == 2):
    x = data2[2, 1:]

    size = len(x)
    MLAT = np.zeros(size)
    req = 7.1492E+07*5.84760
    
    for jj in range(size):
        MLAT0 = 1.
        for ii in range(1000000):
            if (ii == 1000000):
                print("Error!: solution is not found. z_position = " + str(x(jj)))
                
            ff = req * ((1. / 2.) * np.sin(MLAT0) * np.sqrt(3. * np.sin(MLAT0)**2. + 1.) \
                + (1. / (2. * np.sqrt(3.))) * np.arcsinh(np.sqrt(3.) * np.sin(MLAT0))) \
                - x[jj]
            gg = req * np.cos(MLAT0) * np.sqrt(3. * np.sin(MLAT0)**2. + 1.)

            MLAT1 = float(MLAT0 - ff / gg)
            
            if (abs(MLAT1 - MLAT0) <= 1E-5):
                break

            MLAT0 = MLAT1
        
        MLAT[jj] = MLAT1
        print(x[jj], MLAT[jj])

    MLAT = np.rad2deg(MLAT)
    y = data1[1:9, :]*1.E-6
    plt.rcParams["font.size"] = 30
    plt.rcParams.update({'mathtext.default': 'default', 'mathtext.fontset': 'stix'})
    fig = plt.figure()
    ax = fig.add_subplot(111, xlabel='Io← MLAT [degree] →Jupiter', ylabel='Number Density [$cm^{-3}$]', yscale='log', ylim=(1.E-2, 1.E6))   
    ax.plot(MLAT, y[0, :], c='dimgrey', label='$H^+$(Jupiter)', linestyle='solid', linewidth='4')
    ax.plot(MLAT, y[1, :], c='blue', label='$e^-$(Jupiter)', linestyle='dotted', linewidth='4')
    ax.plot(MLAT, y[5, :], c='purple', label='$H^+$(Io)', linestyle='dotted', linewidth='4')
    ax.plot(MLAT, y[2, :], c='orange', label='$O^+$(Io)', linestyle='dotted', linewidth='4')
    ax.plot(MLAT, y[3, :], c='green', label='$S^+$(Io)', linestyle='dotted', linewidth='4')
    ax.plot(MLAT, y[4, :], c='lime', label='$S^{2+}$(Io)', linestyle='dotted', linewidth='4')
    ax.plot(MLAT, y[6, :], c='deepskyblue', label='cold $e^-$(Io)', linestyle='dotted', linewidth='4')
    ax.plot(MLAT, y[7, :], c='hotpink', label='hot $e^-$(Io)', linestyle='dotted', linewidth='4')
    ax.tick_params(labelsize=25)
    #fig.suptitle('Number Density')
    ax.minorticks_on()
    ax.grid(which="both")
    ax.legend()
    

if (channel == 3):
    x = data2[2, 1:]

    size = len(x)
    MLAT = np.zeros(size)
    req = 7.1492E+07*5.84760
    
    for jj in range(size):
        MLAT0 = 1.
        for ii in range(1000000):
            if (ii == 1000000):
                print("Error!: solution is not found. z_position = " + str(x(jj)))
                
            ff = req * ((1. / 2.) * np.sin(MLAT0) * np.sqrt(3. * np.sin(MLAT0)**2. + 1.) \
                + (1. / (2. * np.sqrt(3.))) * np.arcsinh(np.sqrt(3.) * np.sin(MLAT0))) \
                - x[jj]
            gg = req * np.cos(MLAT0) * np.sqrt(3. * np.sin(MLAT0)**2. + 1.)

            MLAT1 = float(MLAT0 - ff / gg)
            
            if (abs(MLAT1 - MLAT0) <= 1E-5):
                break

            MLAT0 = MLAT1
        
        MLAT[jj] = MLAT1
        print(x[jj], MLAT[jj])

    MLAT = np.rad2deg(MLAT)
    plt.rcParams["font.size"] = 30
    plt.rcParams.update({'mathtext.default': 'default', 'mathtext.fontset': 'stix'})
    fig = plt.figure()
    ax = fig.add_subplot(111, xlabel='Io← MLAT [degree] →Jupiter', ylabel=' [$R_J$] ')
    ax.plot(MLAT, (x[len(x)-1]-x)/7.149200000000E+07, linewidth='4')
    ax.plot(MLAT, 0.036*np.ones(len(x)))
    ax.plot(MLAT, 0.910*np.ones(len(x)))
    ax.minorticks_on()
    ax.grid(which="both")


if (channel == 4):
    x = data2[2, 1:]

    size = len(x)
    MLAT = np.zeros(size)
    req = 7.1492E+07*5.84760
    
    for jj in range(size):
        MLAT0 = 1.
        for ii in range(1000000):
            if (ii == 1000000):
                print("Error!: solution is not found. z_position = " + str(x(jj)))
                
            ff = req * ((1. / 2.) * np.sin(MLAT0) * np.sqrt(3. * np.sin(MLAT0)**2. + 1.) \
                + (1. / (2. * np.sqrt(3.))) * np.arcsinh(np.sqrt(3.) * np.sin(MLAT0))) \
                - x[jj]
            gg = req * np.cos(MLAT0) * np.sqrt(3. * np.sin(MLAT0)**2. + 1.)

            MLAT1 = float(MLAT0 - ff / gg)
            
            if (abs(MLAT1 - MLAT0) <= 1E-5):
                break

            MLAT0 = MLAT1
        
        MLAT[jj] = MLAT1
        print(x[jj], MLAT[jj])

    MLAT = np.rad2deg(MLAT)
    y = data1[1:9, :]*1.E-6
    x3 = data3[0, 116:]
    x3 = np.rad2deg(x3)
    y3 = data3[4:14, 116:]*1.E-6
    plt.rcParams["font.size"] = 40
    plt.rcParams.update({'mathtext.default': 'default', 'mathtext.fontset': 'stix'})
    fig = plt.figure()
    ax = fig.add_subplot(111, xlabel='Io← MLAT [degree] →Jupiter', ylabel='Number Density [$cm^{-3}$]', yscale='log', ylim=(1.E-2, 1.E5)) 
    ax.plot(x3, y3[0, :]+y3[2, :]+y3[4, :], c='red', label='$H^+$(PDS)', linestyle='solid', linewidth='4', alpha=0.5)
    ax.plot(x3, y3[1, :]+y3[3, :]+y3[8, :]+y3[9, :], c='blue', label='$e^-$(PDS)', linestyle='solid', linewidth='4', alpha=0.5)
    ax.plot(x3, y3[5, :], c='orange', label='$O^+$(PDS)', linestyle='solid', linewidth='4', alpha=0.5)

    ax.plot(MLAT, y[0, :]+y[5, :], c='red', label='$H^+$(SVC)', linestyle='dotted', linewidth='4')
    ax.plot(MLAT, y[1, :]+y[6, :]+y[7, :], c='blue', label='$e^-$(SVC)', linestyle='dotted', linewidth='4')
    ax.plot(MLAT, y[2, :], c='orange', label='$O^+$(SVC)', linestyle='dotted', linewidth='4')
    ax.minorticks_on()
    ax.grid(which="both")
    ax.legend()

    fig.tight_layout()

if (channel == 5):
    x = data2[2, 1:]

    size = len(x)
    MLAT = np.zeros(size)
    req = 7.1492E+07*5.84760
    
    for jj in range(size):
        MLAT0 = 1.
        for ii in range(1000000):
            if (ii == 1000000):
                print("Error!: solution is not found. z_position = " + str(x(jj)))
                
            ff = req * ((1. / 2.) * np.sin(MLAT0) * np.sqrt(3. * np.sin(MLAT0)**2. + 1.) \
                + (1. / (2. * np.sqrt(3.))) * np.arcsinh(np.sqrt(3.) * np.sin(MLAT0))) \
                - x[jj]
            gg = req * np.cos(MLAT0) * np.sqrt(3. * np.sin(MLAT0)**2. + 1.)

            MLAT1 = float(MLAT0 - ff / gg)
            
            if (abs(MLAT1 - MLAT0) <= 1E-5):
                break

            MLAT0 = MLAT1
        
        MLAT[jj] = MLAT1
        print(x[jj], MLAT[jj])

    MLAT = np.rad2deg(MLAT)
    y = data1[0, :]
    x3 = data3[0, 116:]
    x3 = np.rad2deg(x3)
    y3 = data3[3, 116:]
    plt.rcParams["font.size"] = 40
    plt.rcParams.update({'mathtext.default': 'default', 'mathtext.fontset': 'stix'})
    fig, (ax1, ax2) = plt.subplots(2, 1) #, sharex=True
    
    fig.suptitle('Electrostatic Potential')
    ax1.set_xlabel('Io← MLAT [degree] →Jupiter')
    ax1.set_ylabel('[kV]')
    ax1.plot(x3, y3/1E3, linewidth='4', linestyle='solid', label='PDS', alpha=0.5)
    ax1.plot(MLAT, y/1E3, linewidth='4', linestyle='dotted', label='SVC')
    ax1.minorticks_on()
    ax1.grid(which="both")
    ax1.legend()
    #ax1.tick_params(labelsize=25)
    ax2.set_xlabel('Io← MLAT [degree] →Jupiter')
    ax2.set_ylabel('[V]')
    ax2.plot(x3, y3, linewidth='4', linestyle='solid', label='PDS', alpha=0.5)
    ax2.plot(MLAT, y, linewidth='4', linestyle='dotted', label='SVC')
    ax2.set_ylim(-20, 5)
    ax2.minorticks_on()
    ax2.grid(which="both")
    ax2.legend()
    plt.subplots_adjust(wspace=0.4, hspace=0.6)


plt.tight_layout()
plt.show()
