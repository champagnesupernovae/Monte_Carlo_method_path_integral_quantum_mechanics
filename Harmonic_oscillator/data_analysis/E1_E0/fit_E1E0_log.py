import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit

N = 500
bh = 100
omega = 1

data = open("../../Results/output/C2/bh_100_omega_1_mean_fast/C2_N_500.txt","r")
C2 = np.loadtxt(data, unpack=True)
data.close()

div = 15
L = int(N/div)
C2 = C2[0:L]
eta = bh*omega/N
a = bh/N

dC2 = 0.0005*np.ones(len(C2))/C2
C2 = np.log(C2[0:L])

tau = a*np.linspace(0, L, L)

'''
plt.plot(tau, C2, '.')
plt.show()
'''

def f(x, E, phi):
	return -E*x + phi

xdata = np.linspace(0, bh/div, 1000)
init = [1, 0]

popt, pcov = curve_fit(f, tau, C2, init, dC2)
E_fit, phi_fit = popt
dE_fit, dphi_fit = np.sqrt(pcov.diagonal())

chi2 = (((C2 - f(tau, *popt))/dC2)**2).sum()
ndof = len(tau) - len(init)

#print('A = %f +- %f ' % (A_fit, dA_fit))
print('E1-E0 = %f +- %f ' % (E_fit/omega,  dE_fit/omega))

print('chi2_scaled = %f +/- %f' %(chi2/ndof, np.sqrt(2*ndof)/ndof))

plt.rcParams.update({'font.size': 13})
plt.plot(xdata, f(xdata, *popt), color='red')
plt.errorbar(tau, C2, dC2, fmt='.', color='black') 
plt.xlabel(r'$\tau$') 
plt.ylabel('ln(C2)')
plt.show()