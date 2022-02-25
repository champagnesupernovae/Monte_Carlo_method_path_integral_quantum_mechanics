import numpy as np
import matplotlib.pyplot as plt


file_name = "../../results/field_0/eta_0.01_omega_1/field_0_out_file_bh_10.txt"
data = open(file_name,"r")
field0_10 = np.loadtxt(data, unpack=True)
data.close()

file_name = "../../results/field_100/eta_0.01_omega_1/field_100_out_file_bh_10.txt"
data = open(file_name,"r")
field0_100 = np.loadtxt(data, unpack=True)
data.close()



def f(x, C0):
	return C0*np.exp(-x**2)


C0 = np.sqrt(1./np.pi)
plt.hist(field0_10, bins=100, density=True, alpha=0.5, label='bh=10, field0', color='red')
plt.hist(field0_100, bins=100, density=True, alpha=0.5, label='bh=10, field100', color='blue')
xdata = np.linspace(-5, 5, 1000)
plt.legend()
plt.plot(xdata, f(xdata,C0))
plt.show()