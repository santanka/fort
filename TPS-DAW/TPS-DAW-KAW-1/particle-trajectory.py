import numpy as np
import matplotlib.pyplot as plt
from matplotlib import mathtext
from numpy.lib.twodim_base import tri
mathtext.FontConstantsBase = mathtext.ComputerModernFontConstants

df = np.genfromtxt(r"/home/satanka/Documents/fort/TPS-DAW/TPS-DAW-KAW-1/results_particle/myrank000/particle_trajectory01-102.dat")

#規格化定数
c  = 299792458E0
q  = 1.6021766208E-19
m  = 9.10938356E-31

R_E = 6371E3
L = 9E0
T_ion = 1E3 #[eV]
T_electron = 1E2 #[eV]
moment  = 7.75E22 #the Earth's dipole moment model
mu_0    = 4E0 * np.pi * 1E-7
B0_eq     = (1E-7 * moment) / (L * R_E)**3
Omega0_eq = q * B0_eq / m
z_unit = c / Omega0_eq
t_unit = 1E0 / Omega0_eq
J_unit = m * c**2E0
V_unit = m * c**2E0 / q
#print(t_unit, z_unit)


time = df[:, 1]
z_particle = df[:, 2]/R_E
u_z_particle = df[:, 3]
u_perp_particle = df[:, 4]
u_phase_particle = df[:, 5]
energy_particle = df[:, 6]
pitch_angle_eq = df[:, 7] #deg
v_z_particle = u_z_particle / np.sqrt(1 + (u_z_particle**2 + u_perp_particle**2)/c**2)
v_perp_particle = u_perp_particle / np.sqrt(1 + (u_z_particle**2 + u_perp_particle**2)/c**2)

channel = 1
trigger = 0 #(1: wave_check)

if (trigger == 1):
    df2 = np.genfromtxt(r"/home/satanka/Documents/fort/TPS-DAW/TPS-DAW-KAW-1/results_particle/myrank000/potential_prof.dat")
    z_position = df2[:, 0]
    wave_number_para = df2[:, 1]
    wave_number_perp = df2[:, 2]
    wave_frequency = df2[:, 3]
    V_resonant = df2[:, 4]
    electrostatic_potential = df2[:, 5]
    EE_wave_para = df2[:, 6]
    EE_wave_para_nonphase =  df2[:, 1] * df2[:, 5] * (T_electron/T_ion)
    EE_wave_perp_perp = df2[:, 7]
    EE_wave_perp_phi = df2[:, 8]
    BB_wave_para = df2[:, 9]
    BB_wave_perp = df2[:, 10]
    Alfven_velocity = df2[:, 11]
    V_resonant_wide_plus = V_resonant + np.sqrt(np.abs(q * EE_wave_para_nonphase / m / wave_number_para))/c
    V_resonant_wide_minus = V_resonant - np.sqrt(np.abs(q * EE_wave_para_nonphase / m / wave_number_para))/c
    ion_Larmor_radius = df2[:, 12]

if (channel == 1):
    length = len(z_particle)
    fig = plt.figure()
    plt.rcParams["font.size"] = 40
    plt.rcParams.update({'mathtext.default': 'default', 'mathtext.fontset': 'stix'})
    ax = fig.add_subplot(111, xlabel='z [RE]', ylabel='$v_{\parallel}/c$')    
    ax.plot(z_particle, v_z_particle/c, zorder=1, color='b')
    ax.scatter(z_particle[0], v_z_particle[0]/c, marker='o', color='r', label='start', zorder=3, s=200)
    ax.scatter(z_particle[length-1], v_z_particle[length-1]/c, marker='D', color='r', label='goal', zorder=3, s=200)  #[len(z_particle)-1]
    if (trigger == 1):
        ax.plot(z_position, V_resonant/c, linestyle='-.', color='red', linewidth='4')
        ax.plot(z_position, Alfven_velocity/c, linestyle='-.', color='orange', linewidth='4')
        ax.plot(z_position, V_resonant_wide_plus/c, linestyle='-.', color='green', linewidth='4')
        ax.plot(z_position, V_resonant_wide_minus/c, linestyle='-.', color='green', linewidth='4')
    #fig.suptitle('particle trajectory')
    ax.minorticks_on()
    ax.grid(which="both")
    ax.set_axisbelow(True)
    ax.legend()

if (channel == 2):
    fig = plt.figure()
    plt.rcParams["font.size"] = 40
    ax = fig.add_subplot(111, xlabel='time [s]', ylabel='energy [eV]')    
    ax.plot(time, energy_particle)
    #fig.suptitle('Evolution of particle energy')
    ax.minorticks_on()
    ax.grid(which="both")

if (channel == 3):
    fig = plt.figure()
    plt.rcParams["font.size"] = 40
    ax = fig.add_subplot(111, xlabel='time [s]', ylabel='equatorial pitch angle [degree]')    
    ax.plot(time, pitch_angle_eq)
    #fig.suptitle('Evolution of equatorial pitch angle')
    ax.minorticks_on()
    ax.grid(which="both")

if (channel == 23):
    fig = plt.figure()
    plt.rcParams["font.size"] = 40
    ax1 = fig.add_subplot(121)
    ax1.plot(time[0:900000], energy_particle[0:900000])
    ax1.set_xlabel('time [s]')
    ax1.set_ylabel('energy [eV]')
    ax1.minorticks_on()
    ax1.grid(which="both")
    ax2 = fig.add_subplot(122)
    ax2.plot(time[0:900000], pitch_angle_eq[0:900000])
    ax2.set_xlabel('time [s]')
    ax2.set_ylabel('equatorial pitch angle [degree]')
    ax2.minorticks_on()
    ax2.grid(which="both")

if (channel == 4 and trigger == 1):
    fig = plt.figure()
    plt.rcParams["font.size"] = 40
    ax = fig.add_subplot(111, xlabel='z [RE]', ylabel='wave frequency [Hz]', yscale='log')    
    ax.plot(z_position, wave_frequency/2/np.pi, linewidth='4')
    #fig.suptitle('wave frequency [Hz]')
    ax.minorticks_on()
    ax.grid(which="both")

if (channel == 5 and trigger == 1):
    plt.rcParams["font.size"] = 40
    plt.rcParams.update({'mathtext.default': 'default', 'mathtext.fontset': 'stix'})
    plt.rcParams['text.usetex'] = True
    fig = plt.figure()
    ax = fig.add_subplot(111, xlabel='z [RE]', ylabel='length [km]', yscale='log')    
    ax.plot(z_position, np.abs(2*np.pi/wave_number_para/10**3), label=r'$\lambda_{\parallel}$', linewidth='4')
    ax.plot(z_position, np.abs(2*np.pi/wave_number_perp/10**3), label=r'$\lambda_{\perp}$', linewidth='4')
    ax.plot(z_position, np.abs(ion_Larmor_radius/ 10**3), label=r'$\rho_i$', linewidth='4')
    #ax.plot(z_position, ion_Larmor_radius*wave_number_perp)
    #fig.suptitle('wavelength [km]')
    ax.minorticks_on()
    ax.grid(which="both")
    ax.legend()

if (channel == 6 and trigger == 1):
    fig = plt.figure()
    plt.rcParams["font.size"] = 40
    plt.rcParams.update({'mathtext.default': 'default', 'mathtext.fontset': 'stix'})
    ax = fig.add_subplot(111, xlabel='z [$R_E$]', ylabel='$v_A / c$', yscale='log')    
    ax.plot(z_position, Alfven_velocity/c, linewidth='4')
    #fig.suptitle('Alfven velocity [km/s]')
    ax.minorticks_on()
    ax.grid(which="both")
    

if (channel == 7 and trigger == 1):
    fig = plt.figure()
    plt.rcParams["font.size"] = 40
    plt.rcParams.update({'mathtext.default': 'default', 'mathtext.fontset': 'stix'})
    ax = fig.add_subplot(111, xlabel='z [RE]', ylabel='$E_{\parallel} [mV/m]$')    
    ax.plot(z_position, EE_wave_para * 1E3, linewidth='4')
    #fig.suptitle('Ewpara [mV/m]')
    ax.minorticks_on()
    ax.grid(which="both")
    

if (channel == 8 and trigger == 1):
    fig = plt.figure()
    plt.rcParams["font.size"] = 40
    ax = fig.add_subplot(111, xlabel='z [RE]', ylabel='electrostatic potential [V]')    
    ax.plot(z_position, electrostatic_potential, linewidth='4')
    #fig.suptitle('electrostatic potential [V]')
    ax.minorticks_on()
    ax.grid(which="both")
    



if (channel == 9 and trigger == 1):
    size = len(z_position)
    MLAT = np.zeros(size)
    
    for jj in range(size):
        MLAT0 = 1.
        for ii in range(1000000):
            if (ii == 1000000):
                print("Error!: solution is not found. z_position = " + str(z_position(jj)))
                
            ff = R_E*L * ((1. / 2.) * np.sin(MLAT0) * np.sqrt(3. * np.sin(MLAT0)**2. + 1.) \
                + (1. / (2. * np.sqrt(3.))) * np.log(np.sqrt(3.) * np.sin(MLAT0) + np.sqrt(3. * np.sin(MLAT0)**2. + 1.))) \
                - z_position[jj]*R_E
            gg = R_E*L * np.cos(MLAT0) * np.sqrt(3. * np.sin(MLAT0)**2. + 1.)

            MLAT1 = float(MLAT0 - ff / gg)
            
            if (abs(MLAT1 - MLAT0) <= 1E-5):
                break

            MLAT0 = MLAT1
        
        MLAT[jj] = MLAT1
        print(z_position[jj], MLAT[jj])
    
    Ke_100 = 100 * q
    Ke_1000 = 1000 * q
    cyclotron_radius_100 = 1/q/c/B0_eq * np.cos(MLAT)**6 / np.sqrt(1+3*np.sin(MLAT)**2) * np.sqrt(Ke_100*(Ke_100+m*c*c))
    cyclotron_radius_1000 = 1/q/c/B0_eq * np.cos(MLAT)**6 / np.sqrt(1+3*np.sin(MLAT)**2) * np.sqrt(Ke_1000*(Ke_1000+m*c*c))

    fig = plt.figure()
    plt.rcParams["font.size"] = 40
    ax = fig.add_subplot(111, xlabel='z [RE]', ylabel='length [km]', yscale='log')    
    #ax.plot(z_position, np.abs(2*np.pi/wave_number_para/10**3), label='para')
    ax.plot(z_position, np.abs(2*np.pi/wave_number_perp/10**3), label='wave_perp', linewidth='4')
    ax.plot(z_position, cyclotron_radius_100/10**3, label='100 eV', linewidth='4')
    ax.plot(z_position, cyclotron_radius_1000/10**3, label='1 keV', linewidth='4')
    ax.plot(z_position, np.abs(ion_Larmor_radius/ 10**3), label='ion_acoustic_gyroradius', linewidth='4')
    #ax.plot(z_position, ion_acoustic_gyroradius*wave_number_perp)
    fig.suptitle('length [km]')
    ax.grid()
    ax.legend()


if (channel == 10):
    size = len(z_particle)
    MLAT = np.zeros(size)
    
    for jj in range(size):
        MLAT0 = 1.
        for ii in range(1000000):
            if (ii == 1000000):
                print("Error!: solution is not found. z_position = " + str(z_particle(jj)))
                
            ff = R_E*L * ((1. / 2.) * np.sin(MLAT0) * np.sqrt(3. * np.sin(MLAT0)**2. + 1.) \
                + (1. / (2. * np.sqrt(3.))) * np.log(np.sqrt(3.) * np.sin(MLAT0) + np.sqrt(3. * np.sin(MLAT0)**2. + 1.))) \
                - z_particle[jj]*R_E
            gg = R_E*L * np.cos(MLAT0) * np.sqrt(3. * np.sin(MLAT0)**2. + 1.)

            MLAT1 = float(MLAT0 - ff / gg)
            
            if (abs(MLAT1 - MLAT0) <= 1E-7):
                break

            MLAT0 = MLAT1
        
        MLAT[jj] = MLAT1
    
    mu = m * (v_perp_particle*c)**2. / 2. / B0_eq * np.cos(MLAT)**6 / np.sqrt(1+3*np.sin(MLAT)**2)

    fig = plt.figure()
    plt.rcParams["font.size"] = 40
    plt.rcParams.update({'mathtext.default': 'default', 'mathtext.fontset': 'stix'})
    ax = fig.add_subplot(111, xlabel='time [$s$]', ylabel='1st adiabatic invariant [$Am^2$]')
    ax.plot(time[0:900000], mu[0:900000])
    #fig.suptitle('1st adiabatic invariant [Am^2]')
    ax.minorticks_on()
    ax.grid(which="both")


if (channel == 11 and trigger == 1):
    wave_phase = np.zeros(len(z_position))
    for jj in range(len(z_position)-1):
        wave_phase[jj+1] = wave_phase[jj] + (wave_number_para[jj]+wave_number_para[jj+1])/2E0 * (z_position[jj+1]-z_position[jj])
    fig = plt.figure()
    plt.rcParams["font.size"] = 40
    plt.rcParams.update({'mathtext.default': 'default', 'mathtext.fontset': 'stix'})
    ax = fig.add_subplot(111, xlabel='z [$R_E$]', ylabel='wave phase (t=0) [rad]')
    ax.plot(z_position, wave_phase)
    ax.minorticks_on()
    ax.grid(which="both")


plt.tight_layout()
plt.show()


