import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit

bh = 100
omega = 1

data = open("risultati_fit.txt","r")
DE, dDE, N = np.loadtxt(data, unpack=True)
data.close()

eta = bh*omega/N

plt.rcParams.update({'font.size': 13})
plt.errorbar(eta**2, DE, dDE, fmt='.', color='black')
plt.xlabel(r"$\eta^2$")
plt.ylabel(r"$(E_1-E_0)/(\hbar\omega)$")
plt.show()