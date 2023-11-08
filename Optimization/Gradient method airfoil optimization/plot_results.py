import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("naca4_results.txt",delimiter=' ')
x = np.arange(1,data.shape[0]+1)

# normalize gradient to better view its development
data[:,1] = abs(data[:,1])/max(abs(data[:,1]))
data[:,2] = abs(data[:,2])/max(abs(data[:,2]))
data[:,3] = abs(data[:,3])/max(abs(data[:,3]))

plt.figure(1)
plt.plot(x,data[:,0])
plt.grid(True)
plt.title("L/D")

plt.figure(2)
plt.title("Gradient")
if data.shape[1] == 2:
    plt.plot(x,data[:,1],label="m")
    plt.legend()
    plt.yscale("log")
    plt.grid(True)
    

else:
    plt.plot(x,data[:,1],label="m")
    plt.plot(x,data[:,2],label="p")
    plt.plot(x,data[:,3],label="t")
    plt.legend()
    plt.yscale("log")
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
