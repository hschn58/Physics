import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags, eye
from scipy.sparse.linalg import spsolve

def build_Dx_robin(Nx, dx, alpha, beta):
    """
    Build 1D second‑derivative D_x with Robin BC:
        alpha * u ± beta * u_x = 0
    """
    # central diagonals
    main = -2*np.ones(Nx, dtype=complex)
    off  =  np.ones(Nx-1, dtype=complex)
    D = diags([off, main, off], [-1,0,1], shape=(Nx,Nx), format='lil')

    # modified diagonal entry at boundaries
    star = (alpha/dx)/(alpha/dx + beta) - 2
    D[0,0]   = star
    D[-1,-1] = star
    return D.tocsc()

# 1) PARAMETERS
wavelength = 550e-9       # m
k          = 2*np.pi/wavelength # wavenumber (m^-1)
w0         = 1e-6         # beam waist (m)
dist       = 10           # range in units of w0

x_num = 2000               # transverse points
z_num = 2000               # propagation steps

# physical grids
x_vals = np.linspace(-dist*w0, dist*w0, x_num)
z_vals = np.linspace(0, dist*w0, z_num)
dx = x_vals[1] - x_vals[0]
dz = z_vals[1] - z_vals[0]

# Robin BC constants (Sommerfeld: α = -i k, β = 1)
alpha = -1j * k 
beta  = +1.0

# 2) BUILD SECOND‑DIFFERENCE & CN MATRICES
D = build_Dx_robin(x_num, dx, alpha, beta)

# Crank–Nicolson factor:
fac = 1j * dz / (4 * k * dx**2)

I = eye(x_num, format='csc') # 'Compressed Sparse Column', efficient for large sparse matrices
L = (I - fac * D).tocsc()   # left‑hand side matrix
R = (I + fac * D).tocsc()   # right‑hand side matrix

# 3) INITIALIZE & MARCH
u = np.exp(-x_vals**2 / w0**2).astype(complex)
soln = np.zeros((x_num, z_num), dtype=complex)
soln[:,0] = u

for m in range(z_num-1):
    rhs    = R.dot(u)
    u_next = spsolve(L, rhs) # "Sparse Solve" is used because this solves the system for sparse matrices much faster
    soln[:, m+1] = u_next
    u = u_next

# 4) PLOT REAL PART

for zs in range(z_num):
    for xs in range(x_num):
        soln[xs, zs] = np.real((soln[xs, zs]*np.exp(1j*k*z_vals[zs]) )) 

plt.figure(figsize=(6,5))
plt.imshow(
    np.real(soln),
    extent=[0, dist*w0*1e6, -dist*w0*1e6, dist*w0*1e6],
    origin='lower',
    aspect='auto',
    cmap='binary',
)
plt.colorbar(label='Relative Amplitude')
plt.xlabel('z (µm)')
plt.ylabel('x (µm)')
plt.title('Crank–Nicolson: Electric Field Propagation')
plt.savefig('/Users/henryschnieders/Desktop/Crank_Nicolson.png', dpi=1200)
plt.show()
