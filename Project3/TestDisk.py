import numpy as np
import random
import os
import datetime
import matplotlib.pyplot as plt

# create folder
nowTime = datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S')
folder = os.getcwd() + '\\' + nowTime + '\\' + 'plot'

if not os.path.exists(folder):
    os.makedirs(folder)


G = 1
Max_R = 15
Particle = []
def DistanceParticles(Vector1, Vector2):
    """
    Takes 2 position vectors and calculates distance between them
    :param Vector1: position vector for particle one
    :param Vector2: position vector for particle two
    :return: distance between two particles
    """
    Distance = sum([x * x for x in (Vector1 - Vector2)]) ** 0.5
    return Distance

def make_galaxy(Number_disk, Mass_disk, shift_x, shift_y, shift_z, shift_vx, shift_vy, shift_vz):
    Disk_potential_parameter = 1
    Add_Zero_For_Plot = True
    Minimum_Distance = 0.02
    Particleset = []
    for i in range(Number_disk):
        while True:
            m = np.random.uniform(0, Mass_disk)
            r = ((m / Mass_disk) ** (1 / 3.)) * Disk_potential_parameter * ((1 - (m / Mass_disk) ** (2 / 3.)) ** -0.5)
            if r < Max_R:
                break
        Too_Close = True
        while Too_Close:
            Too_Close = False
            theta = np.random.uniform(0, 2 * np.pi)
            for j in range(len(Particle)):

                if DistanceParticles(np.array(Particle[j][0][0:3]),
                                     np.array([r * np.cos(theta), r * np.sin(theta), 0])) < Minimum_Distance:
                    Too_Close = True

        Position_Mass = []
        Velocity = []
        Position_Mass.append(r * np.cos(theta) + shift_x)
        Position_Mass.append(r * np.sin(theta) + shift_y)
        Position_Mass.append(0 + shift_z)
        Position_Mass.append(Mass_disk / Number_disk)

        Vcirc = (G * m / r) ** 0.5
        Velocity.append(Vcirc * -np.sin(theta) + shift_vx)
        Velocity.append(Vcirc * np.cos(theta) + shift_vy)
        Velocity.append(0 + shift_vz)
        if Add_Zero_For_Plot:
            Velocity.append(0)
        Particleset.append([Position_Mass, Velocity])
    return Particleset

N1 = 1000
N2 = 100
M1 = 1
M2 = 1
disk1 = make_galaxy(N1, M1, 0, 0, 0, 0, 0, 0)
Softening_Parameter = 0.98 * N1 ** (-0.26)
# disk2 = make_galaxy(N2, M2, 10, 0, 0, 0, 0, 0)
print(disk1)
Particle = disk1
# Particle += disk2
PlotParticle1 = np.array(disk1)
# PlotParticle2 = np.array(disk2)
fig1 = plt.figure()
plt.scatter(PlotParticle1[:, 0, 0], PlotParticle1[:, 0, 1], marker='.')
# plt.scatter(PlotParticle2[:, 0, 0], PlotParticle2[:, 0, 1], marker='.')
my_x_ticks = np.arange(-10, 10, 1)
my_y_ticks = np.arange(-10, 10, 1)
plt.xticks(my_x_ticks)
plt.yticks(my_y_ticks)


plt.show()

import numpy as np
from tqdm import tqdm

G = 1  # km*(km/s)/M_Solar
Length = 100
Unit_Converter_Dist = 1
Unit_Converter_Mass = 1
Unit_Converter_Velocity = 1
Steps = 4000
Boxstructure = [[], [], [], [], [], [], [], [], [0, 0, 0, 0]]  # [[],[]] -> [[],[[pos1],[]]]
Center_Of_Mass = [0, 0, 0]
Fillername = Boxstructure
Looping = True
Angle_Criterion = 0.5
Timestep_Size = 0.5


def CoMUpdater(FirstCoM, SecondCoM):
    UpdatedCoM = [0, 0, 0, 0]
    Totalmass = FirstCoM[3] + SecondCoM[3]
    UpdatedCoM[0] = (FirstCoM[3] * FirstCoM[0] + SecondCoM[3] * SecondCoM[0]) / Totalmass
    UpdatedCoM[1] = (FirstCoM[3] * FirstCoM[1] + SecondCoM[3] * SecondCoM[1]) / Totalmass
    UpdatedCoM[2] = (FirstCoM[3] * FirstCoM[2] + SecondCoM[3] * SecondCoM[2]) / Totalmass
    UpdatedCoM[3] = Totalmass
    return UpdatedCoM


def DivisionIndex(Particle, Center_Of_Mass):
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


def TemporaryCoM(DivIndex, Center_Of_Mass, Depth):
    if DivIndex == 0:
        return [Center_Of_Mass[0] - Length / 2 ** (Depth), Center_Of_Mass[1] - Length / 2 ** (Depth),
                Center_Of_Mass[2] - Length / 2 ** (Depth)]
    elif DivIndex == 1:
        return [Center_Of_Mass[0] + Length / 2 ** (Depth), Center_Of_Mass[1] - Length / 2 ** (Depth),
                Center_Of_Mass[2] - Length / 2 ** (Depth)]
    elif DivIndex == 2:
        return [Center_Of_Mass[0] - Length / 2 ** (Depth), Center_Of_Mass[1] + Length / 2 ** (Depth),
                Center_Of_Mass[2] - Length / 2 ** (Depth)]
    elif DivIndex == 3:
        return [Center_Of_Mass[0] + Length / 2 ** (Depth), Center_Of_Mass[1] + Length / 2 ** (Depth),
                Center_Of_Mass[2] - Length / 2 ** (Depth)]
    elif DivIndex == 4:
        return [Center_Of_Mass[0] - Length / 2 ** (Depth), Center_Of_Mass[1] - Length / 2 ** (Depth),
                Center_Of_Mass[2] + Length / 2 ** (Depth)]
    elif DivIndex == 5:
        return [Center_Of_Mass[0] + Length / 2 ** (Depth), Center_Of_Mass[1] - Length / 2 ** (Depth),
                Center_Of_Mass[2] + Length / 2 ** (Depth)]
    elif DivIndex == 6:
        return [Center_Of_Mass[0] - Length / 2 ** (Depth), Center_Of_Mass[1] + Length / 2 ** (Depth),
                Center_Of_Mass[2] + Length / 2 ** (Depth)]
    elif DivIndex == 7:
        return [Center_Of_Mass[0] + Length / 2 ** (Depth), Center_Of_Mass[1] + Length / 2 ** (Depth),
                Center_Of_Mass[2] + Length / 2 ** (Depth)]


def DeeperArray(Div_Index, Temp_CoM, Fillername):
    Index = DivisionIndex(Fillername[Div_Index], Temp_CoM)
    Array = [[] for i in range(9)]
    Array[8] = Fillername[Div_Index]
    Array[Index] = Fillername[Div_Index]
    return Array


def Array_updater(Daughter_Array):
    Updated_CoM = [0, 0, 0, 0]
    for i in range(8):
        if len(Daughter_Array[i]) > 0:
            if isinstance(Daughter_Array[i][0], float):
                Updated_CoM = CoMUpdater(Updated_CoM, Daughter_Array[i])
            else:
                Daughter_Array[i] = Array_updater(Daughter_Array[i])
                Updated_CoM = CoMUpdater(Updated_CoM, Daughter_Array[i][8])
    Daughter_Array[8] = Updated_CoM
    return Daughter_Array


def TreeGenerator():
    Boxstructure = [[], [], [], [], [], [], [], [], [0, 0, 0, 0]]
    for i in range(len(Particle)):
        Looping = True
        Depth = 0
        Center_Of_Mass = [0, 0, 0]
        Fillername = Boxstructure
        for k in range(3):
            if abs(Particle[i][0][k]) > Length:
                Particle[i] = [[random.uniform(-Max_R / 2, Max_R / 2), random.uniform(-Max_R / 2, Max_R / 2),
                                random.uniform(-Max_R / 2, Max_R / 2), 0.0000001], [0.0, 0.0, 0.0]]
        while Looping:

            Div_Index = DivisionIndex(Particle[i][0], Center_Of_Mass)

            if len(Fillername[
                       Div_Index]) > 1:  # This means it is not an empty array [], here we need to update our array correspondingly
                if isinstance(Fillername[Div_Index][0],
                              float):  # If True, these are positions and we need to make new daughter array
                    Temp_CoM = TemporaryCoM(Div_Index, Center_Of_Mass, Depth + 1)
                    Array = DeeperArray(Div_Index, Temp_CoM, Fillername)
                    Fillername[Div_Index] = Array
                else:  # This is a daughter array so we update our CoM and Array and repeat to see where in this array our particle should go
                    Fillername[8] = CoMUpdater(Fillername[8], Particle[i][0])
                    Fillername = Fillername[Div_Index]
                    Center_Of_Mass = TemporaryCoM(Div_Index, Center_Of_Mass, Depth + 1)

                    Depth += 1
            else:  # It is empty so we just insert it
                Fillername[8] = CoMUpdater(Fillername[8], Particle[i][0])
                Fillername[Div_Index] = Particle[i][0]
                Looping = False
    return Boxstructure


def ForceParticles(Vector1, Vector2):
    Distance = DistanceParticles(np.array(Vector1), np.array(Vector2[0:3]))

    return Vector2[3] * G * (np.array(Vector2[0:3]) - np.array(Vector1)) / (Distance**2+Softening_Parameter**2)**1.5


def TotalForce(Fillername, Force, Length, Part_Index, Calculations):
    if (Length / DistanceParticles(np.array(Particle[Part_Index][0][0:3]),
                                   np.array(Fillername[8][0:3]))) < Angle_Criterion:
        Force += ForceParticles(Particle[Part_Index][0][0:3], Fillername[8])
        Calculations += 1
        # print(t)
        # print(Part_Index)
        # print(Length)
        # print(DistanceParticles(np.array(Particle[Part_Index][0][0:3]), np.array(Fillername[8][0:3])))
    else:
        for j in range(8):
            if len(Fillername[j]) > 0:  # Is empty? we skip
                if isinstance(Fillername[j][0], float):  # Is number? We cannot go deeper and do force with this
                    if not Fillername[j] == Particle[Part_Index][0]:  # Except if this was our original particle
                        Force += ForceParticles(Particle[Part_Index][0][0:3], Fillername[j])
                        Calculations += 1
                else:  # We can and must go deeper in the TreeStructure
                    AddedForce, Calculations = TotalForce(Fillername[j], 0, Length / 2, Part_Index, Calculations)
                    Force += AddedForce
    return Force, Calculations


Boxstructure = TreeGenerator()


def VelocityUpdate(Timestep):
    ForceCalculations = [() for i in range(len(Particle))]
    for i in range(len(Particle)):
        Force, ForceCalculations[i] = TotalForce(Boxstructure, 0, Length, i, 0)
        for k in range(3):
            Particle[i][1][k] += Timestep * Force[
                k] * Unit_Converter_Dist ** 2 * Unit_Converter_Mass / Unit_Converter_Velocity ** 2
    return Particle, ForceCalculations


def PositionUpdate(Timestep):
    for i in range(len(Particle)):
        for k in range(3):
            Particle[i][0][k] += Timestep * Particle[i][1][k] * Unit_Converter_Dist
    return Particle


ParticlePosHistory = [[] for k in range(Steps + 2)]
ParticlePosHistory[0] = [Particle[k][0][0:3] for k in range(len(Particle))]

Particle = VelocityUpdate(Timestep=Timestep_Size / 2)[0]
print("Vel update" + str(Particle))
Particle = PositionUpdate(Timestep=Timestep_Size)
print("pos update" + str(Particle))

ParticlePosHistory[1] = [Particle[k][0][0:3] for k in range(len(Particle))]
bar = tqdm(range(Steps))
ForceCalculations = [() for i in range(Steps)]
for t in bar:
    Particle, ForceCalculations[t] = VelocityUpdate(Timestep=Timestep_Size)
    Particle = PositionUpdate(Timestep=Timestep_Size)
    ParticlePosHistory[t + 2] = [Particle[k][0][0:3] for k in range(len(Particle))]
    Boxstructure = Array_updater(Boxstructure)

    if t % 1 == 0:
        Boxstructure = TreeGenerator()
        # print(Boxstructure)
ParticlePosHistory = np.array(ParticlePosHistory)
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation

fig, ax = plt.subplots()
ln1, = plt.plot([], [], 'ro')


def init():
    ax.set_xlim(-Length / 4, Length / 4)
    ax.set_ylim(-Length / 4, Length / 4)


print(ForceCalculations)


def update(q):
    ln1.set_data([ParticlePosHistory[int(10 * q)][:, 0]], [ParticlePosHistory[int(10 * q)][:, 1]])


ani = FuncAnimation(fig, update, frames=int(int(Steps / 10)), interval=20, init_func=init)
ani.save(folder + '\\' + 'animation.gif', writer='pillow')
plt.figure(2)
plt.plot(ParticlePosHistory[:, :, 0], ParticlePosHistory[:, :, 1])
ForceCalculations = np.array(ForceCalculations)
plt.figure(3)
plt.plot(range(Steps), ForceCalculations[:, :])
plt.show()
