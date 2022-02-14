import numpy as np
import matplotlib.pyplot as plt

file_name = "../bootstrap/y2mean_std_eta_bh3_omega3.txt"
data = open(file_name,"r")
y2mean, dy2mean, N = np.loadtxt(data, unpack = True)
data.close()


bh=3
omega=1

eta = bh*omega/N

print(y2mean)

plt.errorbar(eta, y2mean, dy2mean, fmt='.', color='red')
#plt.ylim(0.49,0.58)
plt.show()