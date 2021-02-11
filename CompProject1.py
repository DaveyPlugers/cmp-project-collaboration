import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


#With regards to problem 1, I checked and something is wrong for sure but I will have to take a look at what exactly later


#Current list of problems and improvements to make:
#1: I am unsure if the potential is working correctly, it seems like the interaction might be wrong when looking at close range effects in the animation
#2: We do not have a randomisation for the particle set up and running
#3: We do not have our animation correctly set up for N arbitrary particles but only 2 at the moment
#4: We still have to tinker with h, Timesteps and delay of animation to get it working properly
#5: (minor)  We have h which gives us how long our timestep is and a constant named Timesteps, these 2 might be a bit confusing if we use them through eachother
#6: We do not have anything set up to make the boundary a mirror (L+x -> x)
#7: We do not have anything set up to have interaction with mirror particles (at x-L for example) but this is for later


h = 0.02   #Timestep
L = 100 #Box size
m=3     #Mass of particle
eps = 119
sigma = 3.405
Particleamount = 2
Timesteps = 1000   #Total amount of steps we will take
Particles = np.array([[[20,50],[0.5,0.2]],[[30,50],[-0.2,-0.3]]]) #This array will store the positions and velocities at this time, make random later
                                                                #[Timestep0,Timestep1,Timestep2,...] Will be appended like that where
                                                                #Timestep0 = [Particle1,Particle2] where Particle1 = [pos1vector,vel1vector]

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
    force = np.array([0,0], dtype=np.float64)
    for i in range(Particleamount): #Loop over all the different particles
        if i !=j:
            force += LJpot(Allpositions[tstep][j][0], Allpositions[tstep][i][0])    #We use the Allpositions instead of Particles since these contain already updated values

    return force


for tstep in range(Timesteps):      #Here we let the magic happen, we loop and let the system evolve

    for j in range(Particleamount): #We make sure to go over every particle
        Particles[j,0] = Particles[j,0] + Particles[j,1]*h
        Particles[j,1] = Particles[j,1] + h/m*Force(j)


    Allpositions.append(Particles.copy())       #We save the updated values in our list
    print("timestep= " + str(tstep))            #Not needed but it's here anyways





#Attempt at making animation




fig, ax = plt.subplots()
ln1, = plt.plot([], [], 'ro')
ln2, = plt.plot([], [], 'm*')


def init():
    ax.set_xlim(0, L)
    ax.set_ylim(0, L)


def update(q):
    ln1.set_data([Allpositions[int(q)][0][0][0]], [Allpositions[int(q)][0][0][1]])   #i is the timestep, 0 is particle 1, 0 is for position not velocity, 0 is for x not y, second one we have 1 in the last for y not x
    ln2.set_data([Allpositions[int(q)][1][0][0]], [Allpositions[int(q)][1][0][1]])

ani = FuncAnimation(fig, update,frames=Timesteps+1,interval=10, init_func=init)
plt.show()
