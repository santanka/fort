import numpy as np
import matplotlib.pyplot as plt

df = np.genfromtxt(r"/home/satanka/Documents/fort/TPS-DAW/TPS-DAW-Artemyev/flux_calculation/myrank000-wave_frequency * vA_eq/particle_trajectory00-48.dat")

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
print(t_unit, z_unit)


time = df[:, 1] * t_unit
z_particle = df[:, 2] * z_unit / R_E
u_z_particle = df[:, 3]
u_perp_particle = df[:, 4]
u_phase_particle = df[:, 5] % (2*np.pi) #rad
energy_particle = df[:, 6] * J_unit / q
pitch_angle_eq = df[:, 7] #deg

channel = 3

if (channel == 1):
    fig = plt.figure()
    plt.rcParams["font.size"] = 20
    ax = fig.add_subplot(111, xlabel='z [RE]', ylabel='u_z [/c]')    
    ax.plot(z_particle, u_z_particle)
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


plt.tight_layout()
plt.show()


