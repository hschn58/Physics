import numpy as np
import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import diags, lil_matrix, csc_matrix

def build_Dx_robin(Nx, dx, alpha, beta):
    """
    Build the 1D second‑derivative matrix D_x with a Robin BC at both ends:
        alpha * u +/- beta * u_x = 0.
    Boundary‑row diagonal entries (“star”) come from your formula:
        star = (alpha/dx)/(alpha/dx + beta) - 2.
    """
    # 1) compute the modified diagonal entry at x=±L
    star = (alpha/dx) / (alpha/dx + beta) - 2

    # 2) create the three central diagonals [1, -2, 1]
    main_diag = -2 * np.ones(Nx, dtype=complex)
    off_diag  =  1 * np.ones(Nx-1, dtype=complex)

    # 3) assemble into a sparse LIL matrix
    Dx = diags(
        diagonals=[off_diag, main_diag, off_diag],
        offsets=[-1, 0, 1],
        shape=(Nx, Nx),
        format='lil'
    )

    # 4) overwrite the two boundary rows with the "star" on the diagonal
    Dx[0,   0] = star
    Dx[-1, -1] = star
    # (the off‑diagonal +1 entries at [0,1] and [-1,-2] stay as they are)

    # 5) convert to CSC for fast solves / mat‑vecs
    return Dx


x_num = 450
z_num = 450
wavelength = 550e-9
w0 = 1e-6  # beam waist
dist = 10
k = 2*np.pi/wavelength
dx = 2*dist*w0 / x_num
dz = dist*w0 / z_num
x_vals = np.linspace(-dist*w0, dist*w0, x_num)
z_vals = np.linspace(0, dist*w0, z_num)


alpha = -1j * 2*np.pi/wavelength   # e.g. Sommerfeld: α = -i k
beta  = +1.0                    #         β = 1
Dx    = build_Dx_robin(x_num, dx, alpha, beta)


soln_map = np.zeros((x_num, z_num), dtype=complex)

#center the initial profile in the middle of the image 
initial_data = np.exp(-(x_vals)**2/w0**2) 

soln_map[:, 0] = initial_data[:]

soln_map[:, 1] = 2*soln_map[:, 0] + (3*dz*1j / (k*dx**2)) * (Dx @ soln_map[:, 0]) #first step in z


for zmarch in range(2,z_num): #iterate from x_1 to x_N
    soln_map[:, zmarch] = soln_map[:, zmarch-2] + ((dz*1j) / (2*k*dx**2)) * (Dx @ soln_map[:, zmarch-1]) #march in z

for zs in range(z_num):
    for xs in range(x_num):
        soln_map[xs, zs] = np.real(soln_map[xs, zs]*np.exp(1j*k*z_vals[zs])) #apply the phase factor to the solution

X, Z = np.meshgrid(x_vals, z_vals)

plt.figure(figsize=(6,5))
plt.imshow(np.real(soln_map), extent=[0, dist*w0*1e6, -dist*w0*1e6, dist*w0*1e6], aspect='auto', cmap='binary')
plt.colorbar(label='Relative Amplitude')
plt.title('Leapfrog method: Electric field propagation')
plt.ylabel('x (µm)')
plt.xlabel('z (µm)')

plt.savefig('/Users/henryschnieders/desktop/leapfrog_method.png', dpi=1200, bbox_inches='tight')
plt.show()


