import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit


file_name = "y_mean_vs_bh.txt"
data = open(file_name,"r")
ymean, dymean, bh = np.loadtxt(data, unpack = True)
data.close()

def f(x, m, q):
	return m*x+q

xdata = np.linspace(min(bh)-10, max(bh)+10, 1000)
init = [0, 0]

popt, pcov = curve_fit(f, bh, ymean, init, dymean)
m_fit, q_fit = popt
dm_fit, dq_fit = np.sqrt(pcov.diagonal())

chi2 = (((ymean - f(bh, *popt))/dymean)**2).sum()
ndof = len(bh) - len(init)

print('m = %f +- %f ' % (m_fit, dm_fit))
print('q = %f +- %f ' % (q_fit,  dq_fit))

print('chi2_scaled = %f +/- %f' %(chi2/ndof, np.sqrt(2*ndof)/ndof))

plt.rcParams.update({'font.size': 15})
plt.plot(xdata, f(xdata, *popt), color='red')
plt.errorbar(bh, ymean, dymean, fmt='.', color='black')
plt.xlabel(r"$\beta\hbar\omega$")
plt.ylabel(r"$\langle y\rangle$")
plt.ylim(-0.02,0.02)
plt.yticks(np.arange(-0.02,0.021,0.01))
plt.xlim(min(bh)-10, max(bh)+10)
plt.show()