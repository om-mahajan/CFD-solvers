import numpy as np
import matplotlib.pyplot as plt

# Parameters
M, N = 300, 100  # Grid size
cx, cy = M / 2, N / 2  # Cylinder center
r = min(M, N) / 6  # Cylinder radius

# Create grid
x = np.arange(M+1)
y = np.arange(N+1)
X, Y = np.meshgrid(x, y, indexing='ij')

# Initialize velocity boundary condition arrays
u_bc = np.zeros((M+1, N+1))  # u velocity boundary condition
v_bc = np.zeros((M+1, N+1))  # v velocity boundary condition

# Set velocity boundary conditions (uIn)
u_bc[0, :] = 0        # Left wall u=0
u_bc[M, :] = 1        # Right wall u=1
u_bc[:, 0] = 1        # Bottom u=1
u_bc[:, N] = u_bc[:, N-1]  # Top free slip (approximate)

# Set velocity boundary conditions (vIn)
v_bc[0, :] = 0        # Left wall v=0
v_bc[M-1, :] = 0      # Right wall v=0
v_bc[:, 0] = 0        # Bottom v=0
v_bc[:, N] = 0        # Top v=0

# Mark cylinder obstacle cells
obstacle = (X - cx)**2 + (Y - cy)**2 <= r**2

# Plotting
fig, ax = plt.subplots(figsize=(8, 8))

# Plot velocity boundary conditions as arrows
skip = 3  # skip grid points for clarity
ax.quiver(X[::skip, ::skip], Y[::skip, ::skip], u_bc[::skip, ::skip], v_bc[::skip, ::skip], color='blue', label='Velocity BC')

# Plot cylinder obstacle
circle = plt.Circle((cx, cy), r, color='red', alpha=0.5, label='Cylinder Obstacle')
ax.add_patch(circle)

# Plot pressure boundary condition as dashed box (Neumann BC)
ax.plot([0, M, M, 0, 0], [0, 0, N, N, 0], 'k--', label='Pressure BC (Neumann)')

ax.set_xlim(-1, M+1)
ax.set_ylim(-1, N+1)
ax.set_aspect('equal')
ax.set_title('CFD Boundary Conditions Visualization')

ax.set_xlabel('i (grid index)')
ax.set_ylabel('j (grid index)')
plt.gca().invert_yaxis()  # Invert y-axis to match matrix indexing
plt.show()