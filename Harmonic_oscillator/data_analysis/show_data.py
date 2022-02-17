import numpy as np
import matplotlib.pyplot as plt

data = open("../results/output/mean_y2/bh_1_omega_1/mean_y2_N_10.txt","r")
y2 = np.loadtxt(data, unpack = True)
data.close()

L = len(y2)
l = 1000

means = np.ones(l)

j=0
for i in range(0, L, int(L/l)):
	array = y2[i:i+int(L/l)]
	means[j] = np.mean(array)
	j += 1

M = np.mean(means)
std = np.std(means)/np.sqrt(l-1)

print("mean = ", M, "+/-", std)

plt.figure(1)
plt.plot(np.arange(0,1000), means, '.')
plt.plot(np.arange(0,1000), M*np.ones(1000), color='red')
plt.ylabel("$<y^2>$")
plt.show()