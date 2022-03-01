import numpy as np
import matplotlib.pyplot as plt


file_name = "../../results/field_0/eta_0.01_omega_1/field_0_out_file_bh_100.txt"
data = open(file_name,"r")
data1 = np.loadtxt(data, unpack=True)
data.close()

file_name = "../../results/field_0/eta_0.01_omega_1/field_0_out_file_bh_130.txt"
data = open(file_name,"r")
data2 = np.loadtxt(data, unpack=True)
data.close()

file_name = "../../results/field_0/eta_0.01_omega_1/field_0_out_file_bh_150.txt"
data = open(file_name,"r")
data3 = np.loadtxt(data, unpack=True)
data.close()

def f(x, C0):
	return C0*np.exp(-x**2)


C0 = np.sqrt(1./np.pi)
plt.hist(data1, bins=21, density=True, alpha=0.5, label="bh=100")
plt.hist(data2, bins=21, density=True, alpha=0.5, label="bh=130")
plt.hist(data3, bins=21, density=True, alpha=0.5, label='bh=150')
xdata = np.linspace(-5, 5, 1000)
plt.legend()
plt.plot(xdata, f(xdata,C0))
plt.show()