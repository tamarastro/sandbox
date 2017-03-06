# Code to plot one coordinate system as another.

import numpy as np
import matplotlib.pyplot as plt

c = 3.0e5 #km/s
H0kmsMpc= 70. #km/s/Mpc
nx=11
nR=101
x = np.arange(nx)/(nx-1)*c/H0kmsMpc
R = np.arange(nR)/(nR-1)*2

Dt = np.zeros([nR,nx])

for j in range(nx):
    Dt[:,j] = R*x[j]

plt.plot(range(nR),Dt[:,0])
for j in range(nx):
    plt.plot(range(nR),Dt[:,j])
plt.show()