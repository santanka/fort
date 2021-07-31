import numpy as np


input_file = '/home/satanka/Documents/fort/TPS-DAW/TPS-DAW-Artemyev/flux_calculation/myrank000/particle_trajectory01'
df = np.genfromtxt(input_file + '.dat', filling_values=np.nan)

line_count = 0
for line in df:
    line_count += 1

particle_number = list(set(df[:, 0]))
particle_number_len = len(particle_number)
    
time = line_count / len(particle_number)
time = int(time)


file_name = []
data_name = []
for i in range(particle_number_len):
    file_name.append(input_file + '-' + str(int(particle_number[i])) + '.dat')

xx = np.zeros((time, 8))
for k in range(particle_number_len):
    count = 0
    for l in range(line_count):
        if (df[l, 0] == particle_number[k]):
            xx[count, :] = df[l, :]
            count += 1
        
    np.savetxt(file_name[k], xx)
    xx = np.zeros((time, 8))


