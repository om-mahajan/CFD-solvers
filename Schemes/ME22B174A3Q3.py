#ME22B174
import matplotlib.pyplot as plt
import numpy as np


#Define constants
L = 1
grid_sizex = 50
grid_sizey = 50
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
Dx = gamma*del_x/del_x
Dy = gamma*del_y/del_x

a_e[:, :-1] = Dx
a_n[:-1, :] = Dy

a_s = Fs+ Dy*np.ones(np.shape(Fs))
a_w = Fw + Dx*np.ones(np.shape(Fw))

a_w[:, 0] = 0
a_s[0, :] = 0


a_p = Fe + Fn + 2*(Dx+Dy)
a_p[0,:] += Dy
a_p[:,0] += Dx
a_p[:,grid_sizex-1] -= Dx
a_p[grid_sizey-1,:] -= Dy

phi = np.zeros((grid_sizey, grid_sizex))
it = 0
it_max = 100000
tol = 0.000001
delta = np.inf
conv = False
trunc = False

print("as",a_s)
print("ae",a_e)
print("ap",a_p)
print("an",a_n)
print("aw",a_w)
print("b",b)

phi_history = [phi.copy()]

def refreshb():
    b = np.zeros((grid_sizey,grid_sizex))
    b[0,:] = 1*Fs[0,:] + 2*1*Dy
    for j in range(grid_sizey):
        for i in range(grid_sizex):
            b[j,i] += del_x*del_y*( 2* x_grid[i] + 2*y_grid[j])
    return b

def gauss_seidel(phi_old):
    phi_new = phi_old.copy()
    b = refreshb()
    for j in range(grid_sizey):
        for i in range(grid_sizex): 
            current = phi_new[j,i]
            west  = phi_new[j, i-1] if i-1 >= 0 else 0
            east  = phi_new[j, i+1] if i+1 < grid_sizex else 0
            north = phi_new[j+1, i] if j+1 < grid_sizey else 0
            south = phi_new[j-1, i] if j-1 >= 0 else 0
            west_west = phi_new[j,i-2] if i-2 >= 0 else 0
            south_south = phi_new[j-2, i] if j-2 >= 0 else 0
            #East quick upwinding
            if 0 < i < grid_sizex-1:
                b[j,i] -= Fe[j,i] * (-2/8*(current) + 3/8*(east) -1/8*(west))
            #West quick upwinding
            if i > 1:
                b[j,i] += Fw[j,i] *(-2/8*(west) + 3/8*(current) -1/8*(west_west))
            #South quick upwinding
            if j > 1:
                b[j,i] += Fs[j,i] *(-2/8*(south) + 3/8*(current) -1/8*(south_south))               
            #North quick upwinding
            if 0<j<grid_sizey-1:
                b[j,i] -= Fn[j,i] * (-2/8*(current) + 3/8*(north) -1/8*(south))

            phi_new[j, i] = (a_w[j, i]*west + a_e[j, i]*east + a_n[j, i]*north + a_s[j, i]*south + b[j, i]) / a_p[j, i]

    return phi_new
#np.flipufd]
#driver
while not conv and not trunc:
    print(it,"iteration")
    phi_new = gauss_seidel(phi)
    delta = np.max(np.abs(phi_new - phi))
    phi = np.array(phi_new)
    phi_history.append(phi.copy())  # Save current state
    it += 1
    if delta < tol:
        conv = True
    elif it >= it_max:
        trunc = True
print(phi)

plt.figure(figsize=(8, 6))
plt.imshow(phi, origin='lower', cmap='hot', interpolation='nearest')
plt.colorbar(label='Phi')
plt.title('Heatmap of $\phi$ QUICK 20x20')
plt.xlabel('i (x-direction)')
plt.ylabel('j (y-direction)')
plt.tight_layout()
plt.show()

x_grid = np.linspace(0 + del_x/2, lx - del_x/2, grid_sizex)


center_j = phi.shape[0] // 2  
phi_centerline = phi[center_j, :]

plt.figure(figsize=(8, 5))
plt.plot(x_grid, phi_centerline, marker='o', markersize=5, linewidth=2, label=f"QUICK Scheme ({grid_sizex}x{phi.shape[0]} grid)")

plt.xlabel(r"$x$ coordinate", fontsize=14)
plt.ylabel(r"$\phi$ value", fontsize=14)
plt.title(r"Variation of $\phi$ along horizontal center-line", fontsize=16)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend()
plt.tight_layout()
plt.show()


