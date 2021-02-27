
#Framework of Allposition array: Allposition[i][j][k] =>
#i gives timestep [Timestep0,Timestep1,Timestep2,...] Will be appended like that where =>
#j gives particle number Timestep0 = [Particle1,Particle2] where =>
#k=0 gives position vector, k=1 gives velocity vector Particlej = [posjvector,veljvector]




import numpy as np
import random
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation





#Current list of problems and improvements to make:
#4: We still have to tinker with h, Timesteps and delay of animation to get it working properly
#5: (minor)  We have h which gives us how long our timestep is and a constant named Timesteps, these 2 might be a bit confusing if we use them through eachother
#6 Still did not derive dimensionless expressions or update the code to it yet

#8 Our 3d animation is not running yet (only got the scatter to plot once), found some examples online but got distracted with other problems

#10 (minor and negigeble) Code can crash sometimes, this is due to instability causing a zero division in distance
#11 rewrite position and velocity update to leap-frog (But this will probably be next week's milestone)



#Done tasks
#7: We do not have anything set up to have interaction with mirror particles (at x-L for example) but this is for later
#Done
#9 We need to add the potential and kinetic energy plot
#Done



Tsteplength = 0.002
BoxSize = 15
mass=6
eps = 119
sigma = 3.405
Particleamount = 8
Particlevelocity = 1
dimension = 2
Timesteps = 1000

#Generating particles and Allpositions array
Particles = []
for l in range(Particleamount):

    Partic = [[BoxSize*random.random() for d in range(dimension)],[2*Particlevelocity*(random.random()-0.5) for d in range(dimension)]]
    Particles.append(Partic)


Particles = np.array(Particles).astype('float64')

Allpositions = [[] for x in range(Timesteps+1)]  #Premaking is faster than appending in Matlab not sure about Python
Allpositions[0] = Particles.copy()


def Distancepoints(Vector1,Vector2):
    #Takes 2 position vectors and calculates distance between them
    Distance = sum([x * x for x in ((Vector1 - Vector2 + BoxSize/2)%BoxSize - BoxSize/2)])**0.5
    return Distance



def LJForce(Vector1,Vector2):
    #Takes 2 position vectors and returns the force on Vector1 due to Vector2
    r = Distancepoints(Vector1,Vector2)
    Constant = 24*(2/r**14 - 1/r**8)  #Neglect - since it's in the force, extra 1/r for normalisation Vector1 - Vector2 in next line
    Force = [Constant*x for x in ((Vector1 - Vector2 + BoxSize/2)%BoxSize - BoxSize/2)] #Maybe switch 2 and 1
    return Force


def Force(j,tstep):
    #This calculates the total force vector experienced by particle j
    #Can be improved to exploit symmetry of the force (LJForce(j,i) = LJForce(i,j))
    force = np.array([0 for l in range(dimension)], dtype=np.float64)
    if tstep ==-1: #We don't save position tstep+1 until we have done the velocities as well, this is a workaround to get the updated force
        for i in range(Particleamount):
            if i !=j:
                force += LJForce(Particles[j,0], Particles[i,0])
    else:
        for i in range(Particleamount):
            if i !=j:
                force += LJForce(Allpositions[tstep][j,0], Allpositions[tstep][i][0])    #We use the Allpositions instead of Particles since these contain already updated values

    return force

Epot = []
Ekin = []
Etot = []


def Energy():
    #Can be called to save the potential,kinetic and total energy of the system
    pot=0
    kin=0
    for j in range(Particleamount):
        pot1 = 0
        kin += (0.5*mass)*(sum([x*x for x in Particles[j,1]]))
        for k in range(Particleamount):
            if k !=j:
                r = Distancepoints(Particles[j,0],Particles[k,0])
                pot1 += 4*((1/r)**12 - (1/r)**6)
        pot += pot1/2  # *0.5 to prevent double counting
    Epot.append(pot)
    Ekin.append(kin)
    Etot.append(pot+kin)









for tstep in range(Timesteps):      #Here we let the magic happen, we loop and let the system evolve

    for j in range(Particleamount): #We make sure to update all the particle positions

        Particles[j,0] = Particles[j,0] + Particles[j,1]*Tsteplength +Tsteplength**2/(2*mass)*Force(j,tstep)
        Particles[j,0] = Particles[j, 0] % BoxSize
    for j in range(Particleamount): #Here we update the particle velocities
        Particles[j,1] = Particles[j,1] + Tsteplength/(2*mass)*(Force(j,-1) + Force(j,tstep))

    if tstep%10==0: #For now I do every 10 steps
        Energy()



    Allpositions[tstep + 1] = (Particles.copy())       #We save the updated values in our list

    print("timestep= " + str(tstep))            #Not needed but it's here anyways


#print(Allpositions) #Used for checking mistakes in simulation, remove later


plot2 = plt.figure(2)
plt.plot(range(0,Timesteps,10),Epot)
plt.plot(range(0,Timesteps,10),Ekin)
plt.plot(range(0,Timesteps,10),Etot)
plt.legend(["Epot","Ekin","Etot"])



#Attempt at making animation

fig, ax = plt.subplots()

ln1, = plt.plot([], [], 'ro')




def init():
    ax.set_xlim(0, BoxSize)
    ax.set_ylim(0, BoxSize)


#The 100 I divided and multiplied with is something I will have to manually change, it's so our animation isn't slow when we have many steps


def update(q):
    ln1.set_data([Allpositions[int(10*q)][:,0,0]], [Allpositions[int(10*q)][:,0,1]])   #i is the timestep, 0 is particle 1, 0 is for position not velocity, 0 is for x not y, second one we have 1 in the last for y not x

ani = FuncAnimation(fig, update,frames=int((Timesteps)/10+1),interval=10, init_func=init)
plt.show()

