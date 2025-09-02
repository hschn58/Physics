import sympy as sp
import mpmath as mp
import numpy as np
import matplotlib.pyplot as plt

xpnum = 150
zpnum = 150
dist = 10
# 1) Symbolic definitions
w, x, z, W0, lam = sp.symbols('w x z W0 lam', real=True, positive=True)
k = 2*sp.pi/lam

# distance from source point w to field point (x,z)
r = sp.sqrt(z**2 + (x - w)**2)

# Corrected Gaussian profile and kernel (dropping constant prefactors and obliquity)
E0  = sp.exp(-w**2/W0**2)
kernel = E0 * z * sp.exp(sp.I * k * r) / r**2

# 2) Lambdify to mpmath for fast numerics
integrand = sp.lambdify((w, x, z, W0, lam), kernel, 'mpmath')

def E_numeric(x_val, z_val, W0_val, lam_val=377e-9):
    """
    Numerically computes E(x_val, z_val) via quadrature of the
    kernel over w ∈ [-5 W0, 5 W0].
    """
    L = 5 * W0_val
    f = lambda w_val: integrand(w_val, x_val, z_val, W0_val, lam_val)
    return mp.quad(f, [-L, L])

# 3) Parameters & sampling
W0_val = 1e-6        # 1 μm beam waist
z_val  = 1e-3        # evaluate at z = 1 mm
lam    = 377e-9      # 377 nm

xs     = np.linspace(-dist*W0_val, dist*W0_val, xpnum)
zs     = np.linspace(0, dist*W0_val, zpnum)

soln_map = np.zeros((len(xs), len(zs)), dtype=np.float64)

# 4) Compute field and magnitude

for ziter in range(len(zs)):
    zi = zs[ziter]

    print('ziter', ziter, 'of', len(zs), 'at z =', zi)
    for xiter in range(len(zs)):

        xi = xs[xiter]
        soln_map[xiter, ziter] = np.real(E_numeric(xi, zi, W0_val, lam))


plt.imshow(soln_map, extent=[0, dist*W0_val*1e6, -dist*W0_val*1e6, dist*W0_val*1e6], aspect='auto', cmap='binary')
plt.show()

# E_vals = [E_numeric(xi, z_val, W0_val, lam) for xi in xs]
# E_mag  = np.array([abs(ev) for ev in E_vals])
