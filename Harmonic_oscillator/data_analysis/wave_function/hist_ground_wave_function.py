import numpy as np
import matplotlib.pyplot as plt

N_list = [50,60,70,80,90]

for N in N_list:
	file_name = "./field_0_out_file_N_" + str(N) + ".txt"
	data = open(file_name,"r")
	field_0 = np.loadtxt(data, unpack=True)
	data.close()

	plt.hist(field_0, bins=100, density=True)


plt.show()