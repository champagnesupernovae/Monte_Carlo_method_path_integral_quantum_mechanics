import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit

div = 10
bh = 100
omega = 1

################## N = 120
N = 120
file_name = f"../../results/output/C2/bh_{bh}_omega_{omega}_mean/C2_N_{N}.txt"
data = open(file_name,"r")
C2 = np.loadtxt(data, unpack=True)
data.close()

L = int(N/div)
C2 = C2[0:L]

eta = bh*omega/N
a = bh/N

tau = a*np.linspace(0, L, L)

plt.plot(tau, C2, '.', label="$\eta=0.83$")


################## N = 170
N = 170
file_name = f"../../results/output/C2/bh_{bh}_omega_{omega}_mean/C2_N_{N}.txt"
data = open(file_name,"r")
C2 = np.loadtxt(data, unpack=True)
data.close()

L = int(N/div)
C2 = C2[0:L]

eta = bh*omega/N
a = bh/N

tau = a*np.linspace(0, L, L)

plt.plot(tau, C2, '.', label="$\eta=0.59$")

'''
################## N = 250
N = 250
file_name = f"../../results/output/C2/bh_{bh}_omega_{omega}_mean/C2_N_{N}.txt"
data = open(file_name,"r")
C2 = np.loadtxt(data, unpack=True)
data.close()

L = int(N/div)
C2 = C2[0:L]

eta = bh*omega/N
a = bh/N

tau = a*np.linspace(0, L, L)

plt.plot(tau, C2, '.', label="$\eta=0.4$")


################## N = 500
N = 500
file_name = f"../../results/output/C2/bh_{bh}_omega_{omega}_mean/C2_N_{N}.txt"
data = open(file_name,"r")
C2 = np.loadtxt(data, unpack=True)
data.close()

L = int(N/div)
C2 = C2[0:L]

eta = bh*omega/N
a = bh/N

tau = a*np.linspace(0, L, L)

plt.plot(tau, C2, '.', label="$\eta=0.2$")

'''

plt.xlabel(r"$\tau$")
plt.ylabel("C2")
plt.legend()
plt.rcParams.update({'font.size': 13})
plt.show()