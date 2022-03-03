import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit


data = open("kin_energy_bh_1.txt","r")
U_k1, dU_k1, eta1 = np.loadtxt(data, unpack=True)
data.close()


data = open("kin_energy_bh_100.txt","r")
U_k2, dU_k2, eta2 = np.loadtxt(data, unpack=True)
data.close()


plt.rcParams.update({'font.size': 15})
plt.errorbar(eta1, U_k1, dU_k1, fmt='.', label=r"$\beta\hbar\omega=1$")
plt.errorbar(eta2, U_k2, dU_k2, fmt='.', label=r"$\beta\hbar\omega=100$")
plt.xlabel(r"$\eta$")
plt.ylabel(r"$U_{kin}$")
plt.legend()
plt.show()