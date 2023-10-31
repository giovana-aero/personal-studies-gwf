import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("naca4_results.txt",delimiter=' ')
x = np.arange(1,data.shape[0]+1)

plt.figure(1)
plt.plot(x,data[:,0])
plt.grid(True)
plt.title("L/D")

plt.figure(2)
plt.title("Gradient")
if data.shape[1] == 2:
    plt.plot(x,data[:,1])
    plt.grid(True)
    

else:
    plt.plot(x,data[:,1],x,data[:,2],x,data[:,3])
    plt.grid(True)

foil1 = np.loadtxt("original_coordinates.dat")
foil2 = np.loadtxt("coordinates.dat")

plt.figure(3)
plt.plot(foil1[:,0],foil1[:,1],label="Original")
plt.plot(foil2[:,0],foil2[:,1],label="Optimized")
plt.axis("equal")
plt.grid(True)
plt.legend()
plt.show()
