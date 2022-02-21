import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit

N = 120
bh = 100
omega = 1

file_name = f"../../results/output/C4/bh_{bh}_omega_{omega}_mean_fast/C4_N_{N}.txt"
data = open(file_name,"r")
C4 = np.loadtxt(data, unpack=True)
data.close()

div = 12
L = int(N/div)
C4 = C4[0:L]
eta = bh*omega/N
a = bh/N

dC4 = 0.0005*np.ones(len(C4))/C4
C4 = np.log(C4[0:L])

tau = a*np.linspace(0, L, L)

'''
plt.plot(tau, C4, '.')
plt.show()
'''

def f(x, E, phi):
	return -E*x + phi

xdata = np.linspace(0, bh/div, 1000)
init = [1, 0]

popt, pcov = curve_fit(f, tau, C4, init, dC4)
E_fit, phi_fit = popt
dE_fit, dphi_fit = np.sqrt(pcov.diagonal())

chi2 = (((C4 - f(tau, *popt))/dC4)**2).sum()
ndof = len(tau) - len(init)

#print('A = %f +- %f ' % (A_fit, dA_fit))
print('E2-E0 = %f +- %f ' % (E_fit/omega,  dE_fit/omega))

print('chi2_scaled = %f +/- %f' %(chi2/ndof, np.sqrt(2*ndof)/ndof))

plt.rcParams.update({'font.size': 13})
plt.plot(xdata, f(xdata, *popt), color='red')
plt.errorbar(tau, C4, dC4, fmt='.', color='black') 
plt.xlabel(r'$\tau$') 
plt.ylabel('ln(C4)')
plt.show()