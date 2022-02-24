import numpy as np
import matplotlib.pyplot as plt


file_name = "energy_eta_0.01_omega_1.txt"
data = open(file_name,"r")
U, dU, bh = np.loadtxt(data, unpack = True)
data.close()

eta = 0.001
omega = 1
N = bh*omega / eta;

bh = bh * omega

def f(x):
	return 0.5 + 1/(np.exp(1/x)-1)

xdata = np.linspace(0.001, 12, 1000)

plt.errorbar(1/bh, U, dU, fmt='.', color='red')
plt.plot(xdata, f(xdata), color='black')
plt.xlabel(r'$\beta$E')
plt.ylabel('U')
plt.ylim(0,12)
plt.xlim(0,12)
plt.rcParams.update({'font.size': 15})
plt.show()