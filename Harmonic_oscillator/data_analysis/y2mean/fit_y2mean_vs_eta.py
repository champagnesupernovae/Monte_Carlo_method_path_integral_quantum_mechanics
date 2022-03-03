import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit

data = open("y2mean_std_eta_bh100_omega1.txt","r")
y2, dy2, eta = np.loadtxt(data, unpack=True)
data.close()

bh = 100
omega = 1


def f(x, A, y_exp):
	return A*x**2 + y_exp

xdata = np.linspace(min(eta), max(eta), 1000)
init = [-0.05, 0.5]


popt, pcov = curve_fit(f, eta, y2, init, dy2)
A_fit, y_exp_fit = popt
dA_fit, dy_exp_fit = np.sqrt(pcov.diagonal())

chi2 = (((y2 - f(eta, *popt))/dy2)**2).sum()
ndof = len(eta) - len(init)

y_teo = 0.5 + 1/(np.exp(bh*omega)-1)

print('y_teorico = %lf' %y_teo)
print('A = %f +- %f ' % (A_fit, dA_fit))
print('y_exp = %f +- %f ' % (y_exp_fit,  dy_exp_fit))

print('chi2_scaled = %f +/- %f' %(chi2/ndof, np.sqrt(2*ndof)/ndof))

plt.rcParams.update({'font.size': 15})
plt.plot(xdata, f(xdata, *popt), color='red')
plt.errorbar(eta, y2, dy2, fmt='.', color='black') 
plt.xlabel('$\eta$') 
plt.ylabel(r'$\langle y^2\rangle$')


plt.show()
