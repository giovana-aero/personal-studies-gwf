import numpy as np
import matplotlib.pyplot as plt

def remove_boundaries(mat,Nx,Ny):
    for i in range(Nx):
        mat[0,i] = 0
        mat[-1,i] = 0

    for j in range(Ny):
        mat[j,0] = 0
        mat[j,-1] = 0

    return mat

# Normalizar vetores no gráfico de magnitude de velocidade?
normalize_vectors = 1


data = np.genfromtxt("cavity_flow_p.txt",delimiter = ' ')
lx = data[0,0]
ly = data[0,1]
pmat = data[1:,:]
umat = np.genfromtxt("cavity_flow_u.txt",delimiter = ' ')
vmat = np.genfromtxt("cavity_flow_v.txt",delimiter = ' ')
Nx = pmat.shape[1]
Ny = pmat.shape[0]

x = np.linspace(0,lx,Nx)
y = np.linspace(0,ly,Ny)
X,Y = np.meshgrid(x,y)
V = np.zeros((Ny,Nx))
for j in range(Ny):
    for i in range(Nx):
        V[j,i] = np.sqrt(umat[j,i]**2 + vmat[j,i]**2)

pmat = remove_boundaries(pmat,Nx,Ny)
umat = remove_boundaries(umat,Nx,Ny)
vmat = remove_boundaries(vmat,Nx,Ny)

plt.figure()
fig1 = plt.contourf(X,Y,pmat)
bar = plt.colorbar()
bar.ax.set_ylabel("Normalized pressure")
# plt.show()

plt.figure()
fig2 = plt.contourf(X,Y,umat)
bar = plt.colorbar()
bar.ax.set_ylabel("Horizontal velocity [m/s]")
# plt.show()

plt.figure()
fig3 = plt.contourf(X,Y,vmat)
bar = plt.colorbar()
bar.ax.set_ylabel("Vertical velocity [m/s]")
# plt.show()

if normalize_vectors:
    for j in range(Ny):
        for i in range(Nx):
            if V[j,i] == 0: V[j,i] = 1
    umat = np.divide(umat,V)
    vmat = np.divide(vmat,V)

remove_boundaries(V,Nx,Ny)

plt.figure()
fig4 = plt.contourf(X,Y,V)
bar = plt.colorbar()
bar.ax.set_ylabel("Velocity magnitude [m/s]")
plt.quiver(X,Y,umat,vmat)
plt.show()
