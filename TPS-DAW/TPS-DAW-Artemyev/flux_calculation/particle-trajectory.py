import numpy as np
import matplotlib.pyplot as plt

df = np.genfromtxt(r"/home/satanka/Documents/fort/TPS-DAW/TPS-DAW-Artemyev/flux_calculation/myrank000/particle_trajectory00-3.dat")

#規格化定数
c  = 299792458E0
q  = 1.6021766208E-19
m  = 9.10938356E-31

R_E = 6371E3
L = 9E0
moment  = 7.75E22 #the Earth's dipole moment model
mu_0    = 4E0 * np.pi * 1E-7
B0_eq     = (1E-7 * moment) / (L * R_E)**3
Omega0_eq = q * B0_eq / m
z_unit = c / Omega0_eq
t_unit = 1E0 / Omega0_eq
J_unit = m * c**2E0
V_unit = m * c**2E0 / q
print(t_unit, z_unit)


time = df[:, 1] * t_unit
z_particle = df[:, 2] * z_unit / R_E
u_z_particle = df[:, 3]
u_perp_particle = df[:, 4]
u_phase_particle = df[:, 5] % (2*np.pi) #rad
energy_particle = df[:, 6] * J_unit / q
pitch_angle_eq = df[:, 7] #deg
v_z_particle = u_z_particle / np.sqrt(1 + u_z_particle**2 + u_perp_particle**2)

channel = 7
trigger = 1 #(1: wave_check)

if (trigger == 1):
    df2 = np.genfromtxt(r"/home/satanka/Documents/fort/TPS-DAW/TPS-DAW-Artemyev/flux_calculation/myrank000/potential_prof.dat")
    z_position = df2[:, 0] * z_unit / R_E
    wave_number_para = df2[:, 1] / z_unit
    wave_number_perp = df2[:, 2] / z_unit
    wave_frequency = df2[:, 3] / t_unit
    V_resonant = df2[:, 4]
    electrostatic_potential = df2[:, 5] * V_unit
    EE_wave_para = df2[:, 6] * V_unit / z_unit
    EE_wave_para_nonphase =  - df2[:, 1] * df2[:, 5] * V_unit / z_unit
    Alfven_velocity = df2[:, 7]
    V_resonant_wide_plus = V_resonant + np.sqrt(np.abs(q * EE_wave_para_nonphase / m / wave_number_para))/c
    V_resonant_wide_minus = V_resonant - np.sqrt(np.abs(q * EE_wave_para_nonphase / m / wave_number_para))/c
    ion_acoustic_gyroradius = df2[:, 8] * z_unit

if (channel == 1):
    fig = plt.figure()
    plt.rcParams["font.size"] = 20
    ax = fig.add_subplot(111, xlabel='z [RE]', ylabel='v_z [/c]')    
    ax.plot(z_particle, v_z_particle)
    if (trigger == 1):
        ax.plot(z_position, V_resonant, linestyle='-.', color='red')
        #ax.plot(z_position, Alfven_velocity, linestyle='-.', color='orange')
        ax.plot(z_position, V_resonant_wide_plus, linestyle='-.', color='green')
        ax.plot(z_position, V_resonant_wide_minus, linestyle='-.', color='green')
    fig.suptitle('particle trajectory')
    ax.grid()

if (channel == 2):
    fig = plt.figure()
    plt.rcParams["font.size"] = 20
    ax = fig.add_subplot(111, xlabel='time [s]', ylabel='energy [eV]')    
    ax.plot(time, energy_particle)
    fig.suptitle('Evolution of particle energy')
    ax.grid()

if (channel == 3):
    fig = plt.figure()
    plt.rcParams["font.size"] = 20
    ax = fig.add_subplot(111, xlabel='time [s]', ylabel='equatorial pitch angle [degree]')    
    ax.plot(time, pitch_angle_eq)
    fig.suptitle('Evolution of equatorial pitch angle')
    ax.grid()

if (channel == 4 and trigger == 1):
    fig = plt.figure()
    plt.rcParams["font.size"] = 20
    ax = fig.add_subplot(111, xlabel='z [RE]', ylabel='wave frequency [Hz]', yscale='log')    
    ax.plot(z_position, wave_frequency/2/np.pi)
    fig.suptitle('wave frequency [Hz]')
    ax.grid()

if (channel == 5 and trigger == 1):
    fig = plt.figure()
    plt.rcParams["font.size"] = 20
    ax = fig.add_subplot(111, xlabel='z [RE]', ylabel='wavelength [km]', yscale='log')    
    ax.plot(z_position, np.abs(2*np.pi/wave_number_para/10**3), label='para')
    ax.plot(z_position, np.abs(2*np.pi/wave_number_perp/10**3), label='perp')
    ax.plot(z_position, np.abs(ion_acoustic_gyroradius/ 10**3), label='ion_acoustic_gyroradius')
    #ax.plot(z_position, ion_acoustic_gyroradius*wave_number_perp)
    fig.suptitle('wavelength [km]')
    ax.grid()
    ax.legend()

if (channel == 6 and trigger == 1):
    fig = plt.figure()
    plt.rcParams["font.size"] = 20
    ax = fig.add_subplot(111, xlabel='z [RE]', ylabel='Alfven velocity [km/s]', yscale='log')    
    ax.plot(z_position, Alfven_velocity/10**3)
    fig.suptitle('wave frequency [km/s]')
    ax.grid()

if (channel == 7 and trigger == 1):
    fig = plt.figure()
    plt.rcParams["font.size"] = 20
    ax = fig.add_subplot(111, xlabel='z [RE]', ylabel='Ewpara [mV/m]')    
    ax.plot(z_position, EE_wave_para * 1E3)
    fig.suptitle('Ewpara [V/m]')
    ax.grid()

if (channel == 8 and trigger == 1):
    fig = plt.figure()
    plt.rcParams["font.size"] = 20
    ax = fig.add_subplot(111, xlabel='z [RE]', ylabel='electrostatic potential [V]')    
    ax.plot(z_position, electrostatic_potential)
    fig.suptitle('electrostatic potential [V]')
    ax.grid()

plt.tight_layout()
plt.show()


