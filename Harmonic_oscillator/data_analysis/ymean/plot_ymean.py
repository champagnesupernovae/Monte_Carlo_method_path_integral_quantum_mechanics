import numpy as np
import matplotlib.pyplot as plt



file_name = ".txt"
data = open(file_name,"r")
ymean = np.loadtxt(data, unpack = True)
data.close()


bh = 50
omega = 1
eta = bh*omega/N

