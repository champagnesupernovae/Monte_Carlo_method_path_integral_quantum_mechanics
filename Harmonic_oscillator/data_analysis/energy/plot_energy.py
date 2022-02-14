import numpy as np
import matplotlib.pyplot as plt

file_name = "energy_eta_1_omega_4.txt"
data = open(file_name,"r")
U1, dU1, N1 = np.loadtxt(data, unpack = True)
data.close()

file_name = "energy_eta_0.2_omega_4.txt"
data = open(file_name,"r")
U2, dU2, N2 = np.loadtxt(data, unpack = True)
data.close()

file_name = "energy_eta_0.05_omega_4.txt"
data = open(file_name,"r")
U3, dU3, N3 = np.loadtxt(data, unpack = True)
data.close()


eta = 0.2
bh_omega1 = N1*eta
bh_omega2 = N2*eta
bh_omega3 = N3*eta

def f(x):
	return 0.5 + 1/(np.exp(x)-1)


plt.errorbar(bh_omega1, U1, dU1, fmt='.', label='$\eta$=1')
plt.errorbar(bh_omega2, U2, dU2, fmt='.', label='$\eta$=0.2')
plt.errorbar(bh_omega3, U3, dU3, fmt='.', label='$\eta$=0.05')

xdata = np.linspace(0, max(bh_omega2), 1000)

plt.plot(xdata, f(xdata), color='black')
plt.xlabel(r'$\beta$E')
plt.ylabel('U')
plt.ylim(0,10)
plt.xlim(0,15)
plt.legend()
plt.rcParams.update({'font.size': 15})
plt.show()