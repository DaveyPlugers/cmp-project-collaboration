import numpy as np
import random
import math
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D




#Current list of problems and improvements to make:
#4: We still have to tinker with h, Timesteps and delay of animation to get it working properly
#5: (minor)  We have h which gives us how long our timestep is and a constant named Timesteps, these 2 might be a bit confusing if we use them through eachother
#7: We do not have anything set up to have interaction with mirror particles (at x-L for example) but this is for later
#8 Our 3d animation is not running yet (only got the scatter to plot once), found some examples online but got distracted with other problems
#9 We need to add the potential and kinetic energy plot
#10 (minor and negigeble) Code can crash sometimes, this is due to instability causing a zero division in distance
#11 rewrite position and velocity update to leap-frog (But this will probably be next week's milestone)
#12 Still did not derive dimensionless expressions or update the code to it yet



#Completed problems: (for morale boost)
#1: I am unsure if the potential is working correctly, it seems like the interaction might be wrong when looking at close range effects in the animation
#Kinda resolved, one issue comes that we were using int32 data instead of float64, other is numerical instability (try Particles = np.array([[[20,50],[2,-3]],[[30,50],[-1,-3]]]).astype('float64') if you want to see what happens)
#Should get resolved next week when we implement leap-frog, if not get back to this

#2: We do not have a randomisation for the particle set up and running
#Works now

#3: We do not have our animation correctly set up for N arbitrary particles but only 2 at the moment
#Works fine

#6: We do not have anything set up to make the boundary a mirror (L+x -> x)
#We do now, simple modulo operator



h = 0.005   #Timestep
L = 60 #Box size
m=3     #Mass of particle
eps = 119
sigma = 3.405
Particleamount = 10
Particlevelocity = 3
dimension = 2
Timesteps = 200   #Total amount of steps we will take
#Particles = np.array([[[20,50],[2,-2.8]],[[30,50],[-1,-3]]]).astype('float64') #This array will store the positions and velocities at this time, make random later
                                                                #[Timestep0,Timestep1,Timestep2,...] Will be appended like that where
                                                                #Timestep0 = [Particle1,Particle2] where Particle1 = [pos1vector,vel1vector]
Particles = []
for l in range(Particleamount):

    Partic = [[L*random.random() for d in range(dimension)],[2*Particlevelocity*(random.random()-0.5) for d in range(dimension)]]
    Particles.append(Partic)

print(Particles)
Particles = np.array(Particles).astype('float64')
print(Particles)
Allpositions = []    #Make array where we will save all the positions (and velocity) for all the timesteps
Allpositions.append(Particles.copy())   #Easier if this isn't a numpy array


def Distancepoints(Vector1,Vector2):  #Calculate distance between 2 positions
    Distance = math.sqrt(sum([x * x for x in (Vector2 - Vector1)]))
    return Distance



def LJpot(Vector1,Vector2):   #Can be updated to also include the distance already and then an if to see if distance was given
    r = Distancepoints(Vector1,Vector2)
    Constant = 24*eps*sigma**6*(2*sigma**6/r**14 - 1/r**8)  #Neglect - since we have -deltaU for our force
    Force = [Constant*x for x in Vector1-Vector2] #Maybe switch 2 and 1

    return Force


def Force(j):           #This calculates the force vector, will need improvement later for other mirrored particles
    force = np.array([0 for l in range(dimension)], dtype=np.float64)
    for i in range(Particleamount): #Loop over all the different particles
        if i !=j:
            force += LJpot(Allpositions[tstep][j][0], Allpositions[tstep][i][0])    #We use the Allpositions instead of Particles since these contain already updated values

    return force

def minimage(j): #Takes index, checks if corresp. particle outside box and moves it
    Particles[j,0] = Particles[j, 0] % L


for tstep in range(Timesteps):      #Here we let the magic happen, we loop and let the system evolve

    for j in range(Particleamount): #We make sure to go over every particle

        Particles[j,0] = Particles[j,0] + Particles[j,1]*h
        Particles[j,1] = Particles[j,1] + h/m*Force(j)
        minimage(j) #Dunno why I made this a function, it's a one line code don't think we need it anywhere else?



    Allpositions.append(Particles.copy())       #We save the updated values in our list

    print("timestep= " + str(tstep))            #Not needed but it's here anyways


print(Allpositions) #Used for checking mistakes in simulation, remove later





#Attempt at making animation

fig, ax = plt.subplots()

ln1, = plt.plot([], [], 'ro')




def init():
    ax.set_xlim(0, L)
    ax.set_ylim(0, L)

print(Allpositions[0])
print(3)
print(Allpositions[0][0])
print(4)
print(Allpositions[0][:,0])


def update(q):
    ln1.set_data([Allpositions[int(q)][:,0,0]], [Allpositions[int(q)][:,0,1]])   #i is the timestep, 0 is particle 1, 0 is for position not velocity, 0 is for x not y, second one we have 1 in the last for y not x

ani = FuncAnimation(fig, update,frames=Timesteps+1,interval=10, init_func=init)
plt.show()
