import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit

N = 250

print("N =",N)
bh = 100
omega = 1
file_name = f"../../results/output/C2/bh_{bh}_omega_{omega}_mean/C2_N_{N}.txt"
data = open(file_name,"r")
C2 = np.loadtxt(data, unpack=True)
data.close()

tau_max = 15

eta = bh*omega/N
print(eta)
a = bh/N

tau = a*np.linspace(0, N, N)

mask = tau < tau_max
tau = tau[mask]
C2 = C2[mask]

dC2 = 0.0003*np.ones(len(C2))


'''
plt.plot(tau, C2, '.')
plt.show()
'''

def f(x, A, E, phi):
	return A*np.exp(-E*x) + phi

xdata = a*np.linspace(0, len(tau), 1000)
init = [1, 1, -0.5]

popt, pcov = curve_fit(f, tau, C2, init, dC2)
A_fit, E_fit, phi_fit = popt
dA_fit, dE_fit, dphi_fit = np.sqrt(pcov.diagonal())

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