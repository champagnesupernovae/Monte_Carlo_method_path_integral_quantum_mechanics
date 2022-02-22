import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit

bh = 100
omega = 1

data = open("risultati_fit.txt","r")
DE, dDE, N = np.loadtxt(data, unpack=True)
data.close()

eta2 = (bh*omega/N)**2


###### FIT
def f(x, m, q):
	return m*x+q

init = [0.17, 1]

popt, pcov = curve_fit(f, eta2, DE, init, dDE)
m_fit, q_fit = popt
dm_fit, dq_fit = np.sqrt(pcov.diagonal())

chi2 = (((DE - f(eta2, *popt))/dDE)**2).sum()
ndof = len(eta2) - len(init)

print('m_fit = %f +- %f ' % (m_fit, dm_fit))
print('q_fit = %f +- %f ' % (q_fit, dq_fit))
print('chi2_scaled = %f +/- %f' %(chi2/ndof, np.sqrt(2*ndof)/ndof))


###### PLOT
xdata = np.linspace(0, max(eta2), 1000)
plt.rcParams.update({'font.size': 13})
plt.errorbar(eta2, DE, dDE, fmt='.', color='black')
plt.plot(xdata, f(xdata, *popt), color='red')
plt.xlabel(r"$\eta^2$")
plt.ylabel(r"$(E_1-E_0)/(\hbar\omega)$")
plt.show()