import random
import numpy as np
from tqdm import tqdm

Softening_Parameter = 0.1
N = 1000
M = 0.333
R = 10
a = 0.1
G = 1
x = []
y = []
z = []
vx = []
vy = []
vz = []
rs = []
V = []
Particle = []
Add_Zero_For_Plot = False
for i in range(N):

    while True:
        m = np.random.uniform(0, M)

        r = a / ((M / m) ** 0.5 - 1)
        # print(r)
        if r < R:
            break
    rs.append(r)
    theta = random.random()*2*np.pi
    phi = random.random()*np.pi

    Position_Mass = []
    Velocity = []
    Position_Mass.append(r * np.sin(phi) * np.cos(theta))
    Position_Mass.append(r * np.sin(phi) * np.sin(theta))
    Position_Mass.append(r * np.cos(phi))
    Position_Mass.append(M/N)
    Vcirc = (G*m/r) ** 0.5
    V.append(Vcirc)
    b = random.uniform(-1,1)
    c = random.uniform(-1, 1)

    d = - np.tan(phi)*(b*np.cos(theta) + c*np.sin(theta))

    Norm = np.sqrt(b**2+c**2+d**2)
    Velocity.append(Vcirc * b/Norm)
    Velocity.append(Vcirc * c/Norm)
    Velocity.append(Vcirc * d/Norm)

    if Add_Zero_For_Plot:
        Velocity.append(0)
    Particle.append([Position_Mass,Velocity])

Particle.append([[0.0,0.0,0.0,0.00005],[0.0,0.0,0.0]])
print(Particle)

import numpy as np
from tqdm import tqdm
G = 1
Length = 100
Steps = 50
Boxstructure=[[],[],[],[],[],[],[],[],[0,0,0,0]]         #[[],[]] -> [[],[[pos1],[]]]
#Boxstructure=[[[],[[1,1,1,20],[],[1,1,1,20]],[1,1,1,20]]]
Center_Of_Mass = [0,0,0]
#Particle = [[[3500000000.0+150000000.0,1000000000.0,0.1,0.000003],[0.0,30,0.0]],[[3500000000.0,1000000000.0,0.1,1],[0.0,0.0,0.0]],[[3500000000.0,1000000000.0-108000000.0,0.1,0.000000002],[35,0.0,0.0]],[[3500000000.0,1000000000.0-228000000.0,0.1,0.000000002],[24,0.0,0.0]],[[3500000000.0,1000000000.0-5967265000.0,0.1,0.000000002],[-4.72,0.0,0.0]],[[3500000000.0,1000000000.0-1400000000.0,0.1,0.000000002],[-9.7,0.0,0.0]]]
Fillername = Boxstructure
Looping = True
Angle_Criterion = 0.5
Timestep_Size = 0.1


def CoMUpdater(FirstCoM,SecondCoM):
    UpdatedCoM = [0,0,0,0]
    Totalmass = FirstCoM[3] + SecondCoM[3]
    UpdatedCoM[0] = (FirstCoM[3] * FirstCoM[0] + SecondCoM[3] * SecondCoM[0])/Totalmass
    UpdatedCoM[1] = (FirstCoM[3] * FirstCoM[1] + SecondCoM[3] * SecondCoM[1])/Totalmass
    UpdatedCoM[2] = (FirstCoM[3] * FirstCoM[2] + SecondCoM[3] * SecondCoM[2])/Totalmass
    UpdatedCoM[3] = Totalmass
    return UpdatedCoM

def DivisionIndex(Particle,Center_Of_Mass):
    if Particle[0] < Center_Of_Mass[0]:
        if Particle[1] < Center_Of_Mass[1]:
            if Particle[2] < Center_Of_Mass[2]:
                return 0
            else:
                return 4
        else:
            if Particle[2] < Center_Of_Mass[2]:
                return 2
            else:
                return 6
    else:
        if Particle[1] < Center_Of_Mass[1]:
            if Particle[2] < Center_Of_Mass[2]:
                return 1
            else:
                return 5
        else:
            if Particle[2] < Center_Of_Mass[2]:
                return 3
            else:
                return 7

def TemporaryCoM(DivIndex,Center_Of_Mass,Depth):
    if DivIndex ==0:
        return [Center_Of_Mass[0] - Length / 2 ** (Depth), Center_Of_Mass[1] - Length / 2 ** (Depth),Center_Of_Mass[2] - Length / 2 ** (Depth)]
    elif DivIndex ==1:
        return [Center_Of_Mass[0] + Length / 2 ** (Depth), Center_Of_Mass[1] - Length / 2 ** (Depth),Center_Of_Mass[2] - Length / 2 ** (Depth)]
    elif DivIndex ==2:
        return [Center_Of_Mass[0] - Length / 2 ** (Depth), Center_Of_Mass[1] + Length / 2 ** (Depth),Center_Of_Mass[2] - Length / 2 ** (Depth)]
    elif DivIndex ==3:
        return [Center_Of_Mass[0] + Length / 2 ** (Depth), Center_Of_Mass[1] + Length / 2 ** (Depth),Center_Of_Mass[2] - Length / 2 ** (Depth)]
    elif DivIndex ==4:
        return [Center_Of_Mass[0] - Length / 2 ** (Depth), Center_Of_Mass[1] - Length / 2 ** (Depth),Center_Of_Mass[2] + Length / 2 ** (Depth)]
    elif DivIndex ==5:
        return [Center_Of_Mass[0] + Length / 2 ** (Depth), Center_Of_Mass[1] - Length / 2 ** (Depth),Center_Of_Mass[2] + Length / 2 ** (Depth)]
    elif DivIndex ==6:
        return [Center_Of_Mass[0] - Length / 2 ** (Depth), Center_Of_Mass[1] + Length / 2 ** (Depth),Center_Of_Mass[2] + Length / 2 ** (Depth)]
    elif DivIndex ==7:
        return [Center_Of_Mass[0] + Length / 2 ** (Depth), Center_Of_Mass[1] + Length / 2 ** (Depth),Center_Of_Mass[2] + Length / 2 ** (Depth)]

def DeeperArray(Div_Index,Temp_CoM,Fillername):
    Index = DivisionIndex(Fillername[Div_Index],Temp_CoM)
    Array = [[] for i in range(9)]
    Array[8] = Fillername[Div_Index]
    Array[Index] = Fillername[Div_Index]
    return Array


def Array_updater(Daughter_Array):
    Updated_CoM = [0,0,0,0]
    for i in range(8):
        if len(Daughter_Array[i]) > 0:
            if isinstance(Daughter_Array[i][0],float):
                Updated_CoM = CoMUpdater(Updated_CoM,Daughter_Array[i])
            else:
                Daughter_Array[i] = Array_updater(Daughter_Array[i])
                Updated_CoM = CoMUpdater(Updated_CoM,Daughter_Array[i][8])
    Daughter_Array[8] = Updated_CoM
    return Daughter_Array

def DistanceParticles(Vector1, Vector2):
    """
    Takes 2 position vectors and calculates distance between them
    :param Vector1: position vector for particle one
    :param Vector2: position vector for particle two
    :return: distance between two particles
    """
    Distance = sum([x * x for x in (Vector1 - Vector2)]) ** 0.5
    return Distance



def TreeGenerator():
    Boxstructure = [[], [], [], [],[], [], [], [], [0, 0, 0, 0]]
    for i in range(len(Particle)):
        Looping = True
        Depth=0
        Center_Of_Mass = [0,0,0]
        Fillername = Boxstructure

        while Looping:

            Div_Index = DivisionIndex(Particle[i][0],Center_Of_Mass)


            if len(Fillername[Div_Index]) > 1: #This means it is not an empty array [], here we need to update our array correspondingly
                if isinstance(Fillername[Div_Index][0], float):  #If True, these are positions and we need to make new daughter array
                    Temp_CoM = TemporaryCoM(Div_Index,Center_Of_Mass,Depth+1)
                    Array = DeeperArray(Div_Index,Temp_CoM,Fillername)
                    Fillername[Div_Index] = Array
                else: #This is a daughter array so we update our CoM and Array and repeat to see where in this array our particle should go
                    Fillername[8] = CoMUpdater(Fillername[8],Particle[i][0])
                    Fillername = Fillername[Div_Index]
                    Center_Of_Mass = TemporaryCoM(Div_Index, Center_Of_Mass, Depth+1)

                    Depth += 1
            else: #It is empty so we just insert it
                Fillername[8] = CoMUpdater(Fillername[8],Particle[i][0])
                Fillername[Div_Index] = Particle[i][0]
                Looping = False
    return Boxstructure

def ForceParticles(Vector1,Vector2):
    Distance = DistanceParticles(np.array(Vector1),np.array(Vector2[0:3]))

    return Vector2[3]*G*(np.array(Vector2[0:3])-np.array(Vector1))/(Distance**2+Softening_Parameter**2)**1.5

def TotalForce(Fillername,Force,Length,Part_Index,Calculations):
    if (Length / DistanceParticles(np.array(Particle[Part_Index][0][0:3]), np.array(Fillername[8][0:3]))) < Angle_Criterion:
        Force += ForceParticles(Particle[Part_Index][0][0:3],Fillername[8])
        Calculations += 1
        #print(t)
        #print(Part_Index)
        #print(Length)
        #print(DistanceParticles(np.array(Particle[Part_Index][0][0:3]), np.array(Fillername[8][0:3])))
    else:
        for j in range(8):
            if len(Fillername[j]) > 0: #Is empty? we skip
                if isinstance(Fillername[j][0], float): #Is number? We cannot go deeper and do force with this
                    if not Fillername[j] == Particle[Part_Index][0]: #Except if this was our original particle
                        Force += ForceParticles(Particle[Part_Index][0][0:3],Fillername[j])
                        Calculations += 1
                else:   #We can and must go deeper in the TreeStructure
                    AddedForce, Calculations = TotalForce(Fillername[j],0,Length/2,Part_Index,Calculations)
                    Force += AddedForce
    return Force,Calculations

print("before vel update")

Boxstructure = TreeGenerator()

print("after vel update")
def VelocityUpdate(Timestep):
    ForceCalculations = [() for i in range(len(Particle))]
    for i in range(len(Particle)):
        Force,ForceCalculations[i]  = TotalForce(Boxstructure,0,Length,i,0)
        for k in range(3):
            Particle[i][1][k] += Timestep*Force[k]
    return Particle, ForceCalculations


def PositionUpdate(Timestep):
    for i in range(len(Particle)):
        for k in range(3):
            Particle[i][0][k] += Timestep*Particle[i][1][k]
    return Particle


ParticlePosHistory = [[] for k in range(Steps+2)]
ParticlePosHistory[0] = [Particle[k][0][0:3] for k in range(len(Particle))]

Particle = VelocityUpdate(Timestep=Timestep_Size/2)[0]
print("Vel update" + str(Particle))
Particle = PositionUpdate(Timestep=Timestep_Size)
print("pos update" + str(Particle))

ParticlePosHistory[1] = [Particle[k][0][0:3] for k in range(len(Particle))]
bar = tqdm(range(Steps))
ForceCalculations = [() for i in range(Steps)]
for t in bar:
    Particle, ForceCalculations[t] = VelocityUpdate(Timestep=Timestep_Size)
    Particle = PositionUpdate(Timestep=Timestep_Size)
    ParticlePosHistory[t+2] = [Particle[k][0][0:3] for k in range(len(Particle))]
    Boxstructure = Array_updater(Boxstructure)

    if t%2 == 1:
        Boxstructure = TreeGenerator()
        #print(Boxstructure)
ParticlePosHistory = np.array(ParticlePosHistory)
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation

fig, ax = plt.subplots()
ln1, = plt.plot([], [], 'ro')


def init():
    ax.set_xlim(-Length/30, Length/30)
    ax.set_ylim(-Length/30, Length/30)


print(ForceCalculations)
def update(q):
    ln1.set_data([ParticlePosHistory[int(5*q)][:,0]], [ParticlePosHistory[int(5*q)][:,1]])

ani = FuncAnimation(fig, update, frames=int(int(Steps/5)), interval=50, init_func=init)
plt.figure(2)
plt.plot(ParticlePosHistory[:,:,0],ParticlePosHistory[:,:,1])
ForceCalculations = np.array(ForceCalculations)
plt.figure(3)
plt.plot(range(Steps),ForceCalculations[:,:])
plt.show()


