import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation



"""
This script creates a 2D wave interference animation using the superposition principle.
The script creates a grid of points and calculates the electric field at each point due 
oscillating charge sources (of randomly assigned total charge) set at random points on the grid.
The electric field magnitude dissipates with distance from the source as 1/r for energy conservation.
"""




save_loc = 'path/efield.mp4'



class Ewave:

    def __init__(self,x0,y0,k,omega,E0,phi):
        self.x0=x0
        self.y0=y0
        self.k=k
        self.omega=omega
        self.E0=E0
        self.phi=phi


sources=100

max_x=20
max_y=20
min_y=0
min_x=0

max_omega=np.pi
max_k=3*np.pi
max_E0=5
max_phi=np.pi/2



initials={}
for source in range(sources):
    initials[f'{source}']=Ewave(x0=(max_x-3)*np.random.rand(),
                                y0=(max_y-3)*np.random.rand(),
                                k=max_k*np.random.rand(),
                                omega=max_omega*np.random.rand(),
                                E0=max_E0*np.random.rand(),
                                phi=max_phi*np.random.rand()
                               )

grid_pnum=800
x_span=np.linspace(0,max_x,grid_pnum)
y_span=np.linspace(0,max_y,grid_pnum)

X,Y=np.meshgrid(x_span,y_span)


def eqn(X,x0,Y,y0,t,k,phi,omega,E0,grid_pnum):

  # Create an array to store the results
    dist=np.sqrt((Y - y0)**2 + (X - x0)**2)
    
    processed_array = np.ones_like(dist)  # Initialize with 1

    mask = dist > 1
    processed_array[mask] = 1 / dist[mask]
    
    return E0*processed_array * np.cos(dist * k - omega * t+phi)




def total_sum(initials,grid_pnum,X,Y,t):

    frame_sum=np.zeros((grid_pnum,grid_pnum))
    
    for source in range(len(initials)):

        cvar=initials[f'{source}']
        
        frame_sum+=eqn(X,cvar.x0,Y,cvar.y0,t,cvar.k,cvar.phi,cvar.omega,cvar.E0,grid_pnum)

    return frame_sum


fig, ax = plt.subplots()
cax = ax.imshow(total_sum(initials,grid_pnum,X,Y,0), cmap='RdBu', extent=(min_x, max_x, min_y, max_y))
ax.axis('off')

def animate(t):
    cax.set_data(total_sum(initials,grid_pnum,X,Y,t))
    return cax,


ani = FuncAnimation(fig, animate, frames=np.linspace(0, 20 * np.pi, 1440), interval=15, blit=True)

ani.save('interference2_conserved.mp4')
        
