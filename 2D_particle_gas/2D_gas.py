import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


"""
This version has no trail lines.
This python3 program simulates the movement of a 2D gas based on momentum conservation
Coloumbic interactions are not considered

Datasets:

positions: np.zeros((num_particles, 2, num_steps))    x,y positions for each particle for each step number

col_mat=np.zeros((num_particles, num_particles), dtype=bool)  collision matrix, checks if particle collisions have been accounted for in each step

particles: dictionary of all Particle class objects
"""

save_loc = "path/2d_gas.mp4"    

def are_particles_approaching(rel_pos, vel1, vel2):
   
    rel_vel = np.array(vel2) - np.array(vel1)
    dot_product = np.dot(rel_pos, rel_vel)
    
    return dot_product > 0

class Particle:
    
    def __init__(self, mass, position, velocity, i):
        #define all class attributes
        self.mass = mass
        self.position = np.array(position)
        self.velocity = np.array(velocity)
        self.i=i

    def __delta_check(self, position, velocity, time_step, p_rad): 

        pad_size=3*p_rad
        
        #if particle is in some range between the boundary, check to see if it will fall out
        for i in range(2):

            delta=self.velocity[i]*time_step
            
            if self.position[i]<pad_size:
                if delta<-3*p_rad:
                    
                    self.position[i]=delta-self.position[i]
                    self.position[1-i]+=self.velocity[1-i]*time_step
                    
                    self.velocity[i]=-self.velocity[i]

                    return True

            if self.position[i]>1-pad_size:
                if delta>3*p_rad:

                    self.position[i]=1-(delta-self.position[i])
                    self.position[1-i]+=self.velocity[1-i]*time_step
                    
                    self.velocity[i]=-self.velocity[i]

                    return True
                    
        
    def update_position(self, position,velocity,time_step):
        
        p_rad=0.01470588

        if self.__delta_check(position, velocity, time_step, p_rad)==True:
            return
        
        
        self.position += self.velocity * time_step
        
        # Check for collisions with the grid edges
        for i in range(2):  # Check x and y components
            
            if self.position[i] < 0+2*p_rad:
                if self.velocity[i] <0:
                    self.velocity[i] = -self.velocity[i]  # Reflect the velocity
                    return

            if self.position[i] > 1-2*p_rad:
                if self.velocity[i] >0:
                    self.velocity[i] = -self.velocity[i]
                    return


    def apply_force(self, force, time_step):
        # Update velocity based on force (F = ma)
        acceleration = force / self.mass
        self.velocity += acceleration * time_step
        # Update position after applying force

    def check_collision(self,other,col_mat,i,k):
        
        col_axis=self.position-other.position
        magnitude=np.linalg.norm(col_axis)
        p_rad=0.01470588
        
        if (magnitude <= 1.3*p_rad):
            if are_particles_approaching(rel_pos=col_axis,vel1=self.velocity,vel2=other.velocity)==True:
        
                ucol_axis=col_axis/magnitude
    
                #self component
                #smag+omag=(smag-omag)/v_2,1 +v_2,1
                
                #mv_1,0 + mv_2,0 = mv_1,1 + mv_2,1
                #coefficient of restitution set to 1
                #thus, smag-vmag=-(v_1,1-v_2,1)
        
                smag=np.dot(ucol_axis,self.velocity)
                omag=np.dot(ucol_axis,other.velocity)
                
                omag_out=smag
                smag_out=omag
        
                self.velocity, other.velocity = (smag_out*(ucol_axis)+(self.velocity-smag*ucol_axis)), (omag_out*(ucol_axis)+(other.velocity-omag*ucol_axis))
                # Mark the collision as handled
                col_mat[i, k] = True
                col_mat[k, i] = True  # S
        return
        
# Parameters
mass = 1.0
num_particles = 700
time_step = 0.1
num_steps = 1000
max_vel=0.2
force = np.array([0, 0])  


# Initialize datasets
particles = {}
for i in range(num_particles):
    particles[f"x{i}"] = Particle(
        mass=mass, 
        position=np.random.rand(2), 
        velocity=max_vel * np.random.rand() * np.random.rand(2),
        i=i
    )

positions = np.zeros((num_particles, 2, num_steps))

# Simulate particle movement
for j in range(num_steps):
    
    
    col_mat=np.zeros((num_particles, num_particles), dtype=bool)
    update_check = np.zeros((num_particles),dtype=bool)
    
    for i in range(num_particles):
        
        cvar = particles[f"x{i}"]
        cvar.apply_force(force,time_step)
        
        # Check for collisions with other particles
        for k in range(i + 1, num_particles):
            if not (col_mat[i, k] or k==i):
                particles[f"x{i}"].check_collision(particles[f"x{k}"], col_mat, i, k)

        cvar.update_position(cvar.position, cvar.velocity, time_step)
        
        positions[i, 0, j] = cvar.position[0]
        positions[i, 1, j] = cvar.position[1]

# Set up the figure and axis
fig, ax = plt.subplots()
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_xticks([])
ax.set_yticks([])

# Create initial plot elements
lines = [ax.plot([], [], 'bo-', markersize=5, zorder=2)[0] for _ in range(num_particles)]  # Particle markers

# Initialization function
def init():
    for line in lines:
        line.set_data([], [])
    return lines

# Animation function
def animate(i):
    for j, line in enumerate(lines):
        x = positions[j, 0, :i+1]
        y = positions[j, 1, :i+1]
        line.set_data([positions[j, 0, i]], [positions[j, 1, i]])  # Update particle position as single-element lists
    return lines 

# Create the animation
ani = FuncAnimation(fig, animate, init_func=init, frames=np.arange(num_steps), interval=25, blit=False)

# Save the animation
ani.save(save_loc)