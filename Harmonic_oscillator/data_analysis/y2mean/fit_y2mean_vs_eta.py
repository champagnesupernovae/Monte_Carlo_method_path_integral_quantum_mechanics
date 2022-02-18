import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit

data = open("../../bootstrap/y2mean_std_eta_bh3_omega1.txt","r")
y2, dy2, N = np.loadtxt(data, unpack=True)
data.close()

bh = 3
omega = 1
eta = bh*omega/N


def f(x, A, y_exp):
	return A*x**2 + y_exp

xdata = np.linspace(min(eta), max(eta), 1000)
init = [-0.05, 0.5]


popt, pcov = curve_fit(f, eta, y2, init, dy2)
A_fit, y_exp_fit = popt
dA_fit, dy_exp_fit = np.sqrt(pcov.diagonal())

chi2 = (((y2 - f(eta, *popt))/dy2)**2).sum()
ndof = len(eta) - len(init)

print('A = %f +- %f ' % (A_fit, dA_fit))
print('y_exp = %f +- %f ' % (y_exp_fit,  dy_exp_fit))

print('chi2_scaled = %f +/- %f' %(chi2/ndof, np.sqrt(2*ndof)/ndof))

plt.rcParams.update({'font.size': 15})
plt.plot(xdata, f(xdata, *popt), color='red')
plt.errorbar(eta, y2, dy2, fmt='.', color='black') 
plt.xlabel('$\eta$') 
plt.ylabel(r'$<y^2>$')


plt.show()
