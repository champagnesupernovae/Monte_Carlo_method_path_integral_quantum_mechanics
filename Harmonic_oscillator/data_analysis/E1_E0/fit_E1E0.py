import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit

data = open("../results/output/C2/bh_5_omega_4/C2_N_20.txt","r")
C2, tau = np.loadtxt(data, unpack=True)
data.close()

L = int(len(tau)/2)
print(L)
C2 = C2[0:L]
tau = tau[0:L]

dC2 = 0.0003*np.ones(len(C2))


N = 40
bh = 5
omega = 4
eta = bh*omega/N


def f(x, A, E):
	return A*np.exp(-E*x)

xdata = np.linspace(0, bh/2, 1000)
init = [0.5, 1]


popt, pcov = curve_fit(f, tau, C2, init, dC2)
A_fit, E_fit = popt
dA_fit, dE_fit = np.sqrt(pcov.diagonal())

chi2 = (((C2 - f(tau, *popt))/dC2)**2).sum()
ndof = len(tau) - len(init)

print('A = %f +- %f ' % (A_fit, dA_fit))
print('E1-E0 = %f +- %f ' % (E_fit/omega,  dE_fit/omega))

print('chi2_scaled = %f +/- %f' %(chi2/ndof, np.sqrt(2*ndof)/ndof))

plt.rcParams.update({'font.size': 13})
plt.plot(xdata, f(xdata, *popt), color='red')
plt.errorbar(tau, C2, dC2, fmt='.', color='black') 
plt.xlabel(r'$\tau$') 
plt.ylabel('C2')


plt.show()