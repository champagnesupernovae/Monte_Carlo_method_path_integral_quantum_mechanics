import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit

data = open("Dy2_eta_fixed.txt","r")
Dy2, dDy2, N = np.loadtxt(data, unpack=True)
data.close()

eta = 0.1
omega = 1
#eta = bh*omega / N;

def f(x, A, C):
	return A*1/x + C

diff = eta*np.ones(len(Dy2)) - Dy2

xx = np.linspace(min(N),500, 1000)

plt.plot(N, diff, '.')
plt.plot(xx, f(xx, 0.04, 0.005))
plt.xlabel("N")
plt.ylim(0,0.02)
#plt.ylabel(r"$\eta - \Deltay^2$")
plt.show()