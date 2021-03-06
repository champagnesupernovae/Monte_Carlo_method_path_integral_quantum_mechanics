import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit

N = 500
print("N =",N)
bh = 100
omega = 1

file_name = f"../../results/output/C4/bh_{bh}_omega_{omega}_mean/C4_N_{N}.txt"
data = open(file_name,"r")
C4 = np.loadtxt(data, unpack=True)
data.close()

tau_max = 8

eta = bh*omega/N
a = bh/N

tau = a*np.linspace(0, N, N)

mask = tau < tau_max
tau = tau[mask]
C4 = C4[mask]

dC4 = 0.0001*np.ones(len(C4))

'''
plt.plot(tau, C4, '.')
plt.show()
'''

def f(x, A, E, phi):
	return A*np.exp(-E*x) + phi

xdata = a*np.linspace(0, len(tau), 1000)
init = [1, 2, -0.5]

popt, pcov = curve_fit(f, tau, C4, init, dC4)
A_fit, E_fit, phi_fit = popt
dA_fit, dE_fit, dphi_fit = np.sqrt(pcov.diagonal())

chi2 = (((C4 - f(tau, *popt))/dC4)**2).sum()
ndof = len(tau) - len(init)

print('A = %f +- %f ' % (A_fit, dA_fit))
print('E2-E0 = %f +- %f ' % (E_fit/omega,  dE_fit/omega))

print('chi2_scaled = %f +/- %f' %(chi2/ndof, np.sqrt(2*ndof)/ndof))

plt.rcParams.update({'font.size': 13})
plt.plot(xdata, f(xdata, *popt), color='red')
plt.errorbar(tau, C4, dC4, fmt='.', color='black') 
plt.xlabel(r'$\tau$') 
plt.ylabel('C4')
plt.show()