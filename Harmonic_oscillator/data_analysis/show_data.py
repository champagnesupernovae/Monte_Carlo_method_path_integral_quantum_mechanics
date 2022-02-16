import numpy as np
import matplotlib.pyplot as plt

data = open("../results/output/mean_y2/bh_1_omega_1/mean_y2_N_10.txt","r")
y2 = np.loadtxt(data, unpack = True)
data.close()

mean = np.mean(y2)
std = np.std(y2)

print("mean = ", mean, "+/-", std)

plt.figure(1)
plt.plot(np.arange(0,len(y2)), y2, '.')
plt.plot(np.arange(0,len(y2)), mean*np.ones(len(y2)), color='red')
plt.ylabel("$<y^2>$")
plt.show()