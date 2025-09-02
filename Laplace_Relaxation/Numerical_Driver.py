import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

ROWS = 100
COLS = 500
ITERS = 90000

# set the relative length scales
Lx = 1.0
Ly = 2.0

#hx, hy
dx = Lx / (COLS - 1)
dy = Ly / (ROWS - 1)


u_0 = np.zeros((ROWS,COLS))

# set the boundary conditions ()
u_0[0, :] = 1
u_0[ROWS-1, :] = 0 
u_0[:, 0] = 0
u_0[:, COLS-1] = 0 


# define the analytic solution
def analytic_result(x, y, nmax = 1000):
    
    total = 0

    # Equation is correct, I checked it. 

    for n in range(1, nmax, 2):

        coeff = 4 / (n * np.pi)

        x_comp = np.sin( (n*np.pi*x) / Lx)
        y_comp = np.exp( ( - n*np.pi*y) / Lx)

        total += coeff * x_comp * y_comp

    return total

x = np.linspace(0, Lx, COLS)
y = np.linspace(0, Ly, ROWS)

X, Y = np.meshgrid(x, y)


Z = np.zeros_like(u_0)

# Calculate the analytic solution
for row in range(ROWS):
    for col in range(COLS):
        Z[row, col] = analytic_result(X[row, col], Y[row, col])


denom = (2*(dx**2 + dy**2))

u_k = np.zeros_like(u_0)
u_k[0, :] = 1
u_k[ROWS-1, :] = 0 
u_k[:, 0] = 0
u_k[:, COLS-1] = 0 


# Relaxation approximation block
for iter in range(ITERS):

    if iter % 100 == 0:

        if iter % 1000 == 0:
            print(f"iteration {iter}")
        
        err = Z - u_0
        rms_err = np.sqrt(np.mean(err**2))
        rms_exact = np.sqrt(np.mean(Z**2))
        percent_accuracy = 100 * (1 - rms_err / rms_exact)

        with open('/Users/henryschnieders/desktop/proj_output.txt', 'a') as f:
            f.write(f"Iteration {iter}\n")
            f.write(f"RMS‑based % accuracy: {percent_accuracy:.2f}%\n")

    for row in range(1, ROWS-1):
        for col in range(1, COLS-1):
            
            
            u_k[row, col] = ((dy**2)*(u_0[row, col+1] + u_0[row, col-1]) + 
                            (dx**2)*(u_0[row+1, col] + u_0[row-1, col])) / denom

    u_0 = u_k


x = np.linspace(0, Lx, COLS)
y = np.linspace(0, Ly, ROWS)

X, Y = np.meshgrid(x, y)

fig = plt.figure(figsize=(12, 4))
ax = fig.add_subplot(121, projection='3d')
ax.plot_surface(X, Y, u_0, linewidth=0, antialiased=True)

# Plot the relaxation approximation
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Function Value')
ax.set_title(f'Relaxation approximation after {iter} iterations')

# Plot the analytic solution
ax = fig.add_subplot(122, projection='3d')
ax.plot_surface(X, Y, Z, linewidth=0, antialiased=True)

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Function Value')
ax.set_title(f'Analytic solution')

plt.show()