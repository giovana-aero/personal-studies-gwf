import numpy as np
import matplotlib.pyplot as plt

coo = np.loadtxt("coordinates.dat",delimiter=' ')

plt.figure(1)
plt.plot(coo[:,0],coo[:,1])
plt.grid(True)
plt.axis("Equal")
plt.show()
