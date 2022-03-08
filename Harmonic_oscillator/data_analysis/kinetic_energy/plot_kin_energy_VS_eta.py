import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit


data = open("kin_energy_bh_3.txt","r")
U_k1, dU_k1, eta1 = np.loadtxt(data, unpack=True)
data.close()
U_kin_1_ext = 0.5 * (0.5 + 1/(np.exp(3)-1))
print("energia cinetica attesa per bh=3: \n", U_kin_1_ext)


data = open("kin_energy_bh_100.txt","r")
U_k2, dU_k2, eta2 = np.loadtxt(data, unpack=True)
data.close()
U_kin_2_ext = 0.5 * (0.5 + 1/(np.exp(100)-1))
print("energia cinetica attesa per bh=100: \n", U_kin_2_ext)


plt.rcParams.update({'font.size': 15})
plt.errorbar(eta1, U_k1, dU_k1, fmt='.', label=r"$\beta\hbar\omega=3$")
plt.plot(np.linspace(0, 0.8, 1000), U_kin_1_ext*np.ones(1000), color='red', linewidth=0.8)
plt.errorbar(eta2, U_k2, dU_k2, fmt='.', label=r"$\beta\hbar\omega=100$")
plt.plot(np.linspace(0, 0.8, 1000), U_kin_2_ext*np.ones(1000), color='red', linewidth=0.8)
plt.xlabel(r"$\eta$")
plt.ylabel(r"$U_{kin}$")
plt.xlim(0, 0.8)
plt.legend()
plt.show()