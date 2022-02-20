import numpy as np
import matplotlib.pyplot as plt

k_max = 9
tau_max = 10

file_name = "./bootstrap_matrix_mean_sigma_len_file.txt"
data = open(file_name,"r")
x, dx, len_corr = np.loadtxt(data, unpack = True)
data.close()

len_corr = len_corr[0:k_max]


for i in range(0, tau_max, 1):
	sigma = dx[i*k_max:i*k_max+k_max]
	plt.title('tau = %d' %i)
	plt.plot(len_corr, sigma, '.')
	plt.show()