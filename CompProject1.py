#Framework of Allposition array: Allposition[i][j][k] =>
#i gives timestep [Timestep0,Timestep1,Timestep2,...] Will be appended like that where =>
#j gives particle number Timestep0 = [Particle1,Particle2] where =>
#k=0 gives position vector, k=1 gives velocity vector Particlej = [posjvector,veljvector]


import numpy as np
import random
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


Temperature = 100  #Make this an input variable later
Density= 0.8 #In units, atoms/sigma**3, make input variable later
mass=1 #Mass in atomic mass units



Tsteplength = 0.001
BoxSize = 3*(4*mass/Density)**(1/3)  #Times 4 since 4 particles per cube
print("Boxsize = " + str(BoxSize))




Timesteps = 600
Random = False
#Generating particles and Allpositions array
Particles = []
if Random:
    Particleamount = 8
    dimension = 2
    Particlevelocity = 0.001
    for l in range(Particleamount):

        Partic = [[BoxSize*random.random() for d in range(dimension)],[2*Particlevelocity*(random.random()-0.5) for d in range(dimension)]]
        Particles.append(Partic)



else:
    dimension=3
    Particleamount = 108
    Number = 0
    var = (Temperature / 100) ** (1 / 2)  # 100 comes from eps/k_b = 100K
    Velocities = np.random.normal(0, var, 324)
    for i in range(6):
        for j in range(6):
            for k in range(3):
                Partic = [[k * BoxSize / 3 + (j + i) % 2 * BoxSize / 6, j * BoxSize / 6, i * BoxSize / 6],
                          [Velocities[3 * Number], Velocities[3 * Number + 1], Velocities[3 * Number + 2]]]
                Number += 1
                Particles.append(Partic)
Particles = np.array(Particles).astype('float64')
Allpositions = [[] for x in range(Timesteps+1)]  #Premaking is faster than appending in Matlab not sure about Python
Allpositions[0] = Particles.copy()


def Distancepoints(Vector1,Vector2):
    '''
    Takes 2 position vectors and calculates distance between them
    :param Vector1: position vector for particle one
    :param Vector2: position vector for particle two
    :return: distance between two particles
    '''
    Distance = sum([x * x for x in ((Vector1 - Vector2 + BoxSize/2)%BoxSize - BoxSize/2)])**0.5
    return Distance



def LJForce(Vector1,Vector2):
    '''
    Takes 2 position vectors and returns the Lennard-Jones force on Vector1 due to Vector2
    :param Vector1: position vector for particle one
    :param Vector2: position vector for particle two
    :return: Lennard-Jones force between two particles
    '''
    r = Distancepoints(Vector1,Vector2)
    Constant = 24*(2/r**14 - 1/r**8)  #Neglect - since it's in the force, extra 1/r for normalisation Vector1 - Vector2 in next line
    Force = [Constant*x for x in ((Vector1 - Vector2 + BoxSize/2)%BoxSize - BoxSize/2)] #Maybe switch 2 and 1
    return Force


def Force(j,tstep):
    '''
    This calculates the total force vector experienced by particle j
    :param j: index of a particle
    :param tstep: time step
    :return:total force
    '''
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
    '''
    Can be called to save the potential,kinetic and total energy of the system
    '''
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


def Rescale():
    '''
    Can be called to rescale the kinetic energy of the system
    :return: Scale factor for the velocities
    '''
    VelSquaredSummed = 0
    for j in range(Particleamount):
        VelSquaredSummed += sum([x*x for x in Particles[j,1]])
    Lambda = (((Particleamount-1)*3*Temperature)/(VelSquaredSummed*100))**(1/2)
    return Lambda

def Histogram(Bins):
    '''
    Can be called to calculate a histogram giving number of particles in certain ranges
    :param Bins: Amount of intervals of the histogram
    :return: An array of Histogram values
    '''
    Histo = np.histogram([BoxSize], bins=[BoxSize / (2 * Bins) * x for x in range(Bins + 1)])[0]
    for i in range(Particleamount):
        Distances = [Distancepoints(Particles[1+i+j,0],Particles[i,0]) for j in range(Particleamount-i-1)]  #Skipping some Distance calculations here but compensating the counting by a factor 2 in the end
        Histo += np.histogram(Distances, bins=[BoxSize/(2*Bins)*x for x in range(Bins+1)])[0]
    Histo = 2*Histo
    return Histo









for tstep in range(Timesteps):      #Here we let the magic happen, we loop and let the system evolve

    for j in range(Particleamount): #We make sure to update all the particle positions

        Particles[j,0] = Particles[j,0] + Particles[j,1]*Tsteplength +Tsteplength**2/(2*mass)*Force(j,tstep)
        Particles[j,0] = Particles[j,0] %BoxSize

    for j in range(Particleamount): #Here we update the particle velocities
        Particles[j,1] = Particles[j,1] + Tsteplength/(2*mass)*(Force(j,-1) + Force(j,tstep))


    if tstep%10==0: #For now I do every 10 steps
        Energy()
    if tstep%50==0 or tstep in [10,20,35,80]:
        ScaleFactor = Rescale()
        print("Scalefactor= " + str(ScaleFactor))
        for k in range(Particleamount):
            Particles[k,1] = ScaleFactor*Particles[k,1]
    if tstep==0:
        Histo = Histogram(20)
        print(Histo)



    Allpositions[tstep + 1] = (Particles.copy())       #We save the updated values in our list

    print("timestep = " + str(tstep))            #Not needed but it's here anyways


#print(Allpositions[0:5]) #Used for checking mistakes in simulation, remove later


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

