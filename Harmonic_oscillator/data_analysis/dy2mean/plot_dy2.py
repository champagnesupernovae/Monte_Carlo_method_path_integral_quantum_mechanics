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

plt.errorbar(eta, -0.5*Dy2/eta**2, 0.5*dDy2/eta**2, fmt='.')
plt.show()

print(0.5*dDy2/eta**2)