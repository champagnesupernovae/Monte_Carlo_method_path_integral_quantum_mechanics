import numpy as np
import matplotlib.pyplot as plt

A = list()
N_list = [30,70]

for N in N_list:
	bh = 5
	omega = 4
	eta = 1. / N
	a = bh / N


	file_name = "../results/field/bh_5_omega_4/field_out_file_N_" + str(N) + ".txt"
	data = open(file_name,"r")
	field = np.loadtxt(data, unpack = True)
	data.close()

	A.append(field)
	x = np.arange(0, bh, a)
	#plt.plot(x, field, linewidth=0.0001)
	plt.plot(x, field, marker='.')

#plt.ylim(-2,10)
plt.legend()
plt.show()