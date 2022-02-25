import numpy as np
import matplotlib.pyplot as plt


file_name = "../../results/field_0/bh_100_omega_1/field_0_out_file_N_100.txt"
data = open(file_name,"r")
field0_10 = np.loadtxt(data, unpack=True)
data.close()

file_name = "../../results/field_0/bh_100_omega_1/field_0_out_file_N_200.txt"
data = open(file_name,"r")
field0_100 = np.loadtxt(data, unpack=True)
data.close()

file_name = "../../results/field_0/bh_100_omega_1/field_0_out_file_N_500.txt"
data = open(file_name,"r")
field0_1000 = np.loadtxt(data, unpack=True)
data.close()


def f(x, C0):
	return C0*np.exp(-x**2)


C0 = np.sqrt(1./np.pi)
plt.hist(field0_10, bins=100, density=True, alpha=0.5, label='N=100', color='red')
plt.hist(field0_100, bins=100, density=True, alpha=0.5, label='N=200', color='blue')
plt.hist(field0_1000, bins=100, density=True, alpha=0.5, label='N=500', color='orange')
xdata = np.linspace(-5, 5, 1000)
plt.legend()
plt.plot(xdata, f(xdata,C0))
plt.show()