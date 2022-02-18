import numpy as np
import matplotlib.pyplot as plt

file_name = "bootstrap_old_mean_res_chains_file.txt"
data = open(file_name,"r")
mean = np.loadtxt(data, unpack = True)
data.close()


plt.plot(np.arange(0,len(mean),1), mean, '.')
plt.plot(np.linspace(0,100,1000), np.ones(1000)*np.mean(mean))
plt.show()