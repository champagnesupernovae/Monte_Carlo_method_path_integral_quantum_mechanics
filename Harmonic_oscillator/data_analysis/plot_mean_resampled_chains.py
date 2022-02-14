import numpy as np
import matplotlib.pyplot as plt

input_data = open("mean_res_chains_file.txt","r")
data = np.loadtxt(input_data, unpack = True)
input_data.close()

plt.plot(np.arange(0,len(data),1), data)
Y = np.mean(data)*np.ones(len(data))
plt.plot(np.arange(0,len(data),1), Y)
plt.show()