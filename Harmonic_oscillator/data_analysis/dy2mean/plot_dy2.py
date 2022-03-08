import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit

data = open("Dy2_eta_bh_3.txt","r")
Dy2_bh_3, dDy2_bh_3, eta = np.loadtxt(data, unpack=True)
data.close()

data = open("Dy2_eta_bh_100.txt","r")
Dy2_bh_100, dDy2_bh_100, eta = np.loadtxt(data, unpack=True)
data.close()


plt.errorbar(eta, -0.5*Dy2_bh_3/eta**2, 0.5*dDy2_bh_3/eta**2, fmt='.', color='red', label=r'$\beta\hbar\omega = 3$')
plt.errorbar(eta, -0.5*Dy2_bh_100/eta**2, 0.5*dDy2_bh_100/eta**2, fmt='.', color='green', label=r'$\beta\hbar\omega = 100$')
plt.legend()
plt.title(r'$-\langle\Delta y^2\rangle/2$')
plt.show()