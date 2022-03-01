import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit

data = open("Dy2_N.txt","r")
Dy2, dDy2, N = np.loadtxt(data, unpack=True)
data.close()

bh = 3
omega = 1
eta = bh*omega/N


def f(x):
	return 1/(2*x)

print(f(eta))

yy = f(eta) - 0.5*Dy2/eta**2

plt.errorbar(eta, yy, 0.5*dDy2/eta**2, fmt='.')
plt.plot(np.linspace(min(eta),max(eta),1000), np.zeros(1000), color='red')
plt.show()