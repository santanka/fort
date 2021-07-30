import numpy as np
import matplotlib.pyplot as plt

df = np.genfromtxt(r"/home/satanka/Documents/fort/TPS-DAW/TPS-DAW-Artemyev/flux_calculation/particle_trajectory00-1.dat")

time = df[:, 1]
z_particle = df[:, 2]
u_z_particle = df[:, 3]
u_perp_particle = df[:, 4]
u_phase_particle = df[:, 5] % (2*np.pi)
energy_particle = df[:, 6]
pitch_angle_eq = df[:, 7]

channel = 3

if (channel == 1):
    fig = plt.figure()
    plt.rcParams["font.size"] = 20
    ax = fig.add_subplot(111, xlabel='z', ylabel='u_z')    
    ax.plot(z_particle, u_z_particle)
    fig.suptitle('particle trajectory')
    ax.grid()

if (channel == 2):
    fig = plt.figure()
    plt.rcParams["font.size"] = 20
    ax = fig.add_subplot(111, xlabel='time', ylabel='energy')    
    ax.plot(time, energy_particle)
    fig.suptitle('Evolution of particle energy')
    ax.grid()

if (channel == 3):
    fig = plt.figure()
    plt.rcParams["font.size"] = 20
    ax = fig.add_subplot(111, xlabel='z', ylabel='u_z')    
    ax.plot(time, pitch_angle_eq)
    fig.suptitle('Evolution of equatorial pitch angle')
    ax.grid()


plt.tight_layout()
plt.show()


