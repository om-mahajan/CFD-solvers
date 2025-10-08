import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd
import sys
import math

path = sys.path[0] + "\\u.txt"
uMatrix = pd.read_csv(path, delimiter="\t", on_bad_lines='skip', header=None).to_numpy()

path = sys.path[0] + "\\v.txt"
vMatrix = pd.read_csv(path, delimiter="\t", on_bad_lines='skip', header=None).to_numpy()

path = sys.path[0] + "\\p.txt"
pMatrix = pd.read_csv(path, delimiter="\t", on_bad_lines='skip', header=None).to_numpy()

path = sys.path[0] + "\\speed.txt"
speedMatrix = pd.read_csv(path, delimiter="\t", on_bad_lines='skip', header=None).to_numpy()
path = sys.path[0] + "\\BC.txt"
BCs = pd.read_csv(path, delimiter="\t", on_bad_lines='skip', header=None).to_numpy()

f1 = plt.figure(1)

plt.imshow(uMatrix, cmap='plasma', interpolation='none')#, aspect='auto')
ax = plt.gca()

cbar = plt.colorbar(fraction=0.046, pad=0.04, shrink=0.8)
cbar.set_label('Flow speed (m/s) u')
#plt.clim(-0.5,2)
plt.gca().invert_yaxis()
ax = plt.gca().set_xlim(left=1)
plt.ylabel("y")     
plt.xlabel("x")
plt.title("FLOW OVER A CYLINDER - u") 

f2 = plt.figure(2)

plt.imshow(vMatrix, cmap='plasma', interpolation='none')#, aspect='auto')
ax = plt.gca()

cbar = plt.colorbar(fraction=0.046, pad=0.04, shrink=0.8)
cbar.set_label('Flow speed (m/s) v')
#plt.clim(-0.5,2)
plt.gca().invert_yaxis()
ax = plt.gca().set_xlim(left=1)
plt.ylabel("y") 
plt.xlabel("x")
plt.title("FLOW OVER A CYLINDER - v") 

f3 = plt.figure(3)

plt.imshow(pMatrix, cmap='plasma', interpolation='none')#, aspect='auto')
ax = plt.gca()
cbar = plt.colorbar(fraction=0.046, pad=0.04, shrink=0.8)
cbar.set_label('Pressure (bar)')
plt.gca().invert_yaxis()
ax = plt.gca().set_xlim(left=1)
plt.ylabel("y") 
plt.xlabel("x")

plt.title("FLOW OVER A CYLINDER - p") 

#plt.savefig("C:\\Users\\Niels\\Google Drive\\Cpp\\SIMPLE\\pressure.png",bbox_inches='tight')

f4 = plt.figure(4)

plt.imshow(BCs, cmap='plasma', interpolation='none')#, aspect='auto')
ax = plt.gca()
cbar = plt.colorbar()
plt.gca().invert_yaxis()
ax = plt.gca().set_xlim(left=1)
plt.ylabel("y") 
plt.xlabel("x")

plt.title("FLOW OVER A CYLINDER - BCs") 

#plt.savefig("C:\\Users\\Niels\\Google Drive\\Cpp\\SIMPLE\\pressure.png",bbox_inches='tight')

f5 = plt.figure(5)

ny, nx = uMatrix.shape  # Get matrix dimensions

# Create meshgrid for quiver positions
ny, nx = uMatrix.shape
x = np.linspace(0, nx - 1, nx)
y = np.linspace(0, ny - 1, ny)
X, Y = np.meshgrid(x, y)

# Compute speed magnitude
speed = np.sqrt(uMatrix**2 + vMatrix**2)

# Downsample for clarity
step = 5
Xq = X[::step, ::step]
Yq = Y[::step, ::step]
Uq = uMatrix[::step, ::step]
Vq = vMatrix[::step, ::step]
speedQ = speed[::step, ::step]

# Plot quiver
plt.figure(figsize=(10, 6))
quiv = plt.quiver(Xq, Yq, Uq, Vq, speedQ, cmap='plasma', scale=5)
plt.colorbar(quiv, label='Speed (m/s)',fraction=0.046, pad=0.04)

# Optional: overlay BCs
plt.imshow(BCs, cmap='plasma', alpha=0.3, interpolation='none')

plt.gca().invert_yaxis()
plt.title("FLOW OVER A CYLINDER - quiver plot")
plt.xlabel("x")
plt.ylabel("y")
plt.tight_layout()
plt.show()

#x = np.arange(0, uMatrix.shape[0], 10)
#y = np.arange(0, vMatrix.shape[0], 10)

#print(np.add.outer(x,y))

#plt.quiver(np.add.outer(x,y), np.add.outer(x,y), uMatrix[x,y], vMatrix[x,y], scale = 10)

#plt.quiver(uMatrix, vMatrix, scale = 100)
#plt.contour(speedMatrix, cmap='plasma', levels=50, linewidths=1.5)
f6 = plt.figure(6)
plt.imshow(speedMatrix, cmap='plasma', interpolation='bicubic')#, aspect='auto')
ax = plt.gca()
rect = patches.Rectangle((31, -2), 10, 20, linewidth=0, edgecolor='b', facecolor='w')
#circ = patches.Circle((101, 198), 5, facecolor='w')
#ax.add_patch(rect)
#ax.add_patch(circ)



cbar = plt.colorbar(fraction=0.046, pad=0.04, shrink=0.8)
cbar.set_label('Flow speed (m/s)')
#plt.clim(-0.5,2)
plt.gca().invert_yaxis()
ax = plt.gca().set_xlim(left=1)
plt.ylabel("y") 
plt.xlabel("x")
plt.title("FLOW OVER A CYLINDER - v magnitude") 

x = np.linspace(0, nx - 1, nx)
y = np.linspace(0, ny - 1, ny)
X, Y = np.meshgrid(x, y)

# Compute speed magnitude
speed = np.sqrt(uMatrix**2 + vMatrix**2)

plt.figure(figsize=(10, 6))

# Plot streamlines with color based on speed
stream = plt.streamplot(X, Y, uMatrix, vMatrix, color=speed, cmap='plasma', density=1.2, linewidth=1)
plt.colorbar(stream.lines, label='Speed (m/s)',fraction=0.046, pad=0.04)

# Optional: overlay BCs if needed
# plt.imshow(BCs, cmap='gray', alpha=0.2, interpolation='none')

plt.gca().invert_yaxis()
plt.title("FLOW OVER A CYLINDER - streamlines")
plt.xlabel("x")
plt.ylabel("y")
plt.tight_layout()
plt.show()
#plt.xlim([2,100])

#plt.savefig("C:\\Users\\Niels\\Google Drive\\Cpp\\SIMPLE\\velocity.png",bbox_inches='tight')
#plt.legend(['c$_1$', 'c$_2$'])
#plt.legend(['CH$_4$', 'O$_2$', 'H$_2$', 'H$_2$O', 'CO', 'CO$_2$', 'C'])
#plt.legend(['DMC', 'MeCN'])
#plt.legend(['$\phi = 0.5$', '$\phi = 0.8$', '$\phi = 1$', '$\phi = 1.2$', '$\phi = 1.4$', '$\phi = 1.6$'])
#plt.legend(['Full rate', '1/10', '1/20', '1/40'], loc=1)
#plt.legend(['Steel', 'Aluminium'])
#plt.xticks([0,1,2,3,4,5,6,7,8,9])

#plt.axvline(x=1450, ls = '--', label='axvline - full height')
#plt.axvline(x=0.08, ls = '--', label='axvline - full height')

plt.show()

    
    
