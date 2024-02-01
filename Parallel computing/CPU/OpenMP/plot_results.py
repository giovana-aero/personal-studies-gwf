import matplotlib.pyplot as plt
import numpy as np

data = np.genfromtxt("results.txt",delimiter=' ')
x = np.linspace(0,data[0][0],data.shape[1])
y = np.linspace(0,data[0][1],data.shape[1])
X,Y = np.meshgrid(x,y)

fig = plt.contourf(X,np.flip(Y,0),data[1:,:])
bar = plt.colorbar()
bar.ax.set_ylabel("T [degrees celsius]")
plt.show()

# plt.contour()

