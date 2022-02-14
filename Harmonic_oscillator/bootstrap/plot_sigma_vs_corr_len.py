import numpy as np
import matplotlib.pyplot as plt

file_name = "./bootstrap_mean_sigma_len_file.txt"
data = open(file_name,"r")
x, dx, len_corr = np.loadtxt(data, unpack = True)
data.close()

print(len_corr)

plt.plot(len_corr, dx, '.')
plt.show()