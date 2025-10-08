#ME22B174
import matplotlib.pyplot as plt
import numpy as np

#Define constants
L = 1
grid_sizex = 3
grid_sizey = 3
gamma = 3
rho = 2
lx = 3
ly = 3
del_x = lx/grid_sizex
del_y = ly/grid_sizey
phi = 0.5

#initialiazing matrices

a_p = np.zeros((grid_sizey,grid_sizex))
a_e = np.zeros((grid_sizey,grid_sizex))
a_w = np.zeros((grid_sizey,grid_sizex))
a_n = np.zeros((grid_sizey,grid_sizex))
a_s = np.zeros((grid_sizey,grid_sizex))
b = np.zeros((grid_sizey,grid_sizex))

Fe = np.zeros((grid_sizey,grid_sizex))
Fw = np.zeros((grid_sizey,grid_sizex))
Fn = np.zeros((grid_sizey,grid_sizex))
Fs = np.zeros((grid_sizey,grid_sizex))

x_grid = np.linspace(0+del_x/2,lx-del_x/2,grid_sizex)
y_grid = np.linspace(0+del_y/2,ly-del_y/2,grid_sizey)
# python - A[y,x]
#filling the matrices
for j in range(grid_sizey):
    for i in range(grid_sizex):
        #center = (grid_sizex[i], grid_sizey[j])
        #F = rho*u*A
        Fe[j,i] = ((x_grid[i] + del_x/2)**2 + 1)*rho*del_y*L
        Fw[j,i] = ((x_grid[i] - del_x/2)**2 + 1)*rho*del_y*L
        Fn[j,i] = ((y_grid[j] + del_y/2)**2 + 1)*rho*del_y*L
        Fs[j,i] = ((y_grid[j] - del_y/2)**2 + 1)*rho*del_y*L

Dx = 0
Dy = 0

a_e[:, :-1] = Dx
a_n[:-1, :] = Dy

a_s = Fs+ Dy*np.ones(np.shape(Fs))
a_w = Fw + Dx*np.ones(np.shape(Fw))
b[0,:] = 1*Fs[0,:]

for j in range(grid_sizey):
    for i in range(grid_sizex):
        #aw = Dw + Fw
        b[j,i] += del_x*del_y*( 2* x_grid[i] + 2*y_grid[j])

#ap = aw + as + b
a_p = Fe + Fn

phi = np.ones((grid_sizey, grid_sizex))*100
it = 0
it_max = 100000
tol = 0.000001
delta = np.inf
conv = False
trunc = False

# Create a history list
phi_history = [phi.copy()]

def gauss_seidel(phi_old):
    phi_new = phi_old.copy()
    for j in range(grid_sizey):
        for i in range(grid_sizex):
            west  = phi_new[j, i-1] if i-1 >= 0 else 0
            east  = phi_new[j, i+1] if i+1 < grid_sizex else 0
            north = phi_new[j+1, i] if j+1 < grid_sizey else 0
            south = phi_new[j-1, i] if j-1 >= 0 else 0
            phi_new[j, i] = (a_w[j, i]*west + a_e[j, i]*east + a_n[j, i]*north + a_s[j, i]*south + b[j, i]) / a_p[j, i]
    return phi_new

while not conv and not trunc:
    print(it,"iteration")
    phi_new = gauss_seidel(phi)
    print(phi)
    delta = np.max(np.abs(phi_new - phi))
    phi = np.array(phi_new)
    phi_history.append(phi.copy())  
    it += 1
    if delta < tol:
        conv = True
    elif it >= it_max:
        trunc = True

#PLOTTING

plt.figure(figsize=(8, 6))
plt.imshow(phi, origin='lower', cmap='hot', interpolation='nearest')


plt.colorbar(label='Phi')
plt.title('Heatmap of $\phi$')
plt.xlabel('i (x-direction)')
plt.ylabel('j (y-direction)')


nrows, ncols = phi.shape
for j in range(nrows):
    for i in range(ncols):
        plt.text(i, j, f"{phi[j, i]:.2f}", 
                 ha='center', va='center', color='white', fontsize=8)

plt.tight_layout()
plt.show()
x_grid = np.linspace(0 + del_x/2, lx - del_x/2, grid_sizex)


center_j = phi.shape[0] // 2  
phi_centerline = phi[center_j, :]

plt.figure(figsize=(8, 5))
plt.plot(x_grid, phi_centerline, marker='o', markersize=5, linewidth=2, label=f"Upwind Scheme ({grid_sizex}x{phi.shape[0]} grid)")

plt.xlabel(r"$x$ coordinate", fontsize=14)
plt.ylabel(r"$\phi$ value", fontsize=14)
plt.title(r"Variation of $\phi$ along horizontal center-line", fontsize=16)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend()
plt.tight_layout()
plt.show()

