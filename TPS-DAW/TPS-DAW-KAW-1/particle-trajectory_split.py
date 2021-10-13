import numpy as np


input_file = '/home/satanka/Documents/fort/TPS-DAW/TPS-DAW-KAW-1/results_particle/myrank000/particle_trajectory00'
df = np.genfromtxt(input_file + '.dat', filling_values=np.nan)

line_count = 0
for line in df:
    line_count += 1
print(line_count)

particle_number = list(set(df[:, 0]))
particle_number_len = len(particle_number)
print(particle_number_len)

    
time = line_count / particle_number_len
time = int(time)


file_name = []
data_name = []
for i in range(particle_number_len):
    file_name.append(input_file + '-' + str(int(particle_number[i])) + '.dat')

xx = np.zeros((time, 8))
for k in range(particle_number_len):
    print(file_name[k])
    count = 0
    for l in range(line_count):
        if (l % particle_number_len == k):
            xx[count, :] = df[l, :]
            count += 1
        
    np.savetxt(file_name[k], xx)
    xx = np.zeros((time, 8))


