import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit

omega=1

data = open("kin_energy_eta_0.1.txt","r")
U_k1, dU_k1, N1 = np.loadtxt(data, unpack=True)
data.close()
eta1=0.1
bh1 = N1*eta1/omega

data = open("kin_energy_eta_0.01.txt","r")
U_k2, dU_k2, N2 = np.loadtxt(data, unpack=True)
data.close()
eta2=0.01
bh2 = N2*eta2/omega

plt.rcParams.update({'font.size': 15})
plt.errorbar(1/bh1, U_k1, dU_k1, fmt='.', label=r"$\eta=0.1$")
plt.errorbar(1/bh2, U_k2, dU_k2, fmt='.', label=r"$\eta=0.01$")
plt.xlabel(r"$1/\beta\hbar$")
plt.ylabel(r"$U_{kin}$")
plt.legend()
plt.show()