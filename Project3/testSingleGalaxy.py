import os
import random
import datetime
import argparse
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# create folder
nowTime = datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S')
folder = os.getcwd() + '\\' + nowTime + '\\' + 'plot'

if not os.path.exists(folder):
    os.makedirs(folder)

folder1 = os.getcwd() + '\\' + nowTime + '\\' + 'data'

if not os.path.exists(folder1):
    os.makedirs(folder1)

# parser part
parser = argparse.ArgumentParser()
parser.add_argument('--ParticleNumber',
                    help='Value of particle number of bulge, default: 100',
                    nargs=1)

parser.add_argument('--Step',
                    help='Value of total steps, default: 100',
                    nargs=1)

parser.add_argument('--Length',
                    help='Value of length of boxsize, default: 100',
                    nargs=1)

parser.add_argument('--MultipleMass',
                    help='Nfold of the mass for galaxy(original mass per galaxy is 7.4*10^10 solar mass), default: 1',
                    nargs=1)

parser.add_argument('--EnergyTracking',
                    help='Analyse the systems energy',
                    action='store_true')

args = parser.parse_args()

Energy_Tracking = args.EnergyTracking

if args.ParticleNumber is None:
    N = 100
else:
    N = int(args.ParticleNumber[0])

if args.Step is None:
    Steps = 100
else:
    Steps = int(args.Step[0])

if args.Length is None:
    Length = 100
else:
    Length = int(args.Length[0])

if args.MultipleMass is None:
    NFold = 1
else:
    NFold = int(args.MultipleMass[0])

# Write a txt output profile
f = open(folder1 + '\\' + 'output.txt', 'w')
f.write('Initial condition: ParticleNumber = ' + str(N*4) + ' Steps = ' + str(Steps))
f.write('\nThe mass of the galaxy =' + str(NFold) + '*7.4*10^10 solar mass')
f.write('\nThe length of the box =' + str(Length))
f.close()


def make_galaxy(Number_disk, Mass_disk, Number_bulge, Mass_bulge):
    """
    :param Number_disk: galaxy disk particle number
    :param Mass_disk: galaxy disk mass
    :param Number_bulge: galaxy bulge particle number
    :param Mass_bulge: galaxy bulge mass
    :return: galaxy particle set
    """
    Disk_potential_parameter = 1
    bulge_scale_length = 0.1
    Add_Zero_For_Plot = True
    Particleset = []
    for i in range(Number_disk):
        while True:
            m = np.random.uniform(0, Mass_disk)
            r = ((m / Mass_disk) ** (1 / 3.)) * Disk_potential_parameter * ((1 - (m / Mass_disk) ** (2 / 3.)) ** -0.5)
            if r < Max_R:
                break

        theta = np.random.uniform(0, 2 * np.pi)
        Position_Mass = []
        Velocity = []
        Position_Mass.append(r * np.cos(theta))
        Position_Mass.append(r * np.sin(theta))
        Position_Mass.append(0)
        Position_Mass.append(Mass_disk / Number_disk)

        Vcirc = (G * m / r) ** 0.5
        Velocity.append(Vcirc * -np.sin(theta))
        Velocity.append(Vcirc * np.cos(theta))
        Velocity.append(0)
        if Add_Zero_For_Plot:
            Velocity.append(0)
        Particleset.append([Position_Mass, Velocity])

    for i in range(Number_bulge):
        while True:
            m = np.random.uniform(0, Mass_bulge)
            r = bulge_scale_length / ((Mass_bulge / m) ** 0.5 - 1)
            if r < Max_R:
                break

        theta = np.random.uniform(0, 2 * np.pi)
        phi = np.random.uniform(0, np.pi)
        Position_Mass = []
        Velocity = []
        Position_Mass.append(r * np.sin(phi) * np.cos(theta))
        Position_Mass.append(r * np.sin(phi) * np.sin(theta))
        Position_Mass.append(r * np.cos(phi))
        Position_Mass.append(Mass_bulge / Number_bulge)
        Vcirc = (G * m / r) ** 0.5
        tan_x = random.uniform(-1, 1)
        tan_y = random.uniform(-1, 1)

        tan_z = - np.tan(phi) * (tan_x * np.cos(theta) + tan_y * np.sin(theta)) #Particles still go clockwise and counter-clockwise, use cross of position and velocity
        #and take inner product with (0,0,1) I think and check the sign to determine clockwise or counterclockwise and then make sure all particles have same rotation. (r x v)·(0,0,1) > or < 0

        Norm = np.sqrt(tan_x ** 2 + tan_y ** 2 + tan_z ** 2)
        Velocity.append(Vcirc * tan_x / Norm)
        Velocity.append(Vcirc * tan_y / Norm)
        Velocity.append(Vcirc * tan_z / Norm)
        if Add_Zero_For_Plot:
            Velocity.append(0)
        Particleset.append([Position_Mass, Velocity])
    return Particleset



def DistanceParticles(Vector1, Vector2):
    """
    Takes 2 position vectors and calculates distance between them
    :param Vector1: position vector for particle one
    :param Vector2: position vector for particle two
    :return: distance between two particles
    """
    Distance = sum([x * x for x in (Vector1 - Vector2)]) ** 0.5
    return Distance


def CoMUpdater(FirstCoM, SecondCoM):
    """
    Takes 2 position + mass vectors and calculates the position and mass of their center of mass.
    :param FirstCoM: Position and mass vector for particle one
    :param SecondCoM: Position and mass vector for particle one
    :return: UpdatedCoM: Vector with the position of the center of mass and the combined mass
    """
    UpdatedCoM = [0, 0, 0, 0]
    Totalmass = FirstCoM[3] + SecondCoM[3]
    UpdatedCoM[0] = (FirstCoM[3] * FirstCoM[0] + SecondCoM[3] * SecondCoM[0]) / Totalmass
    UpdatedCoM[1] = (FirstCoM[3] * FirstCoM[1] + SecondCoM[3] * SecondCoM[1]) / Totalmass
    UpdatedCoM[2] = (FirstCoM[3] * FirstCoM[2] + SecondCoM[3] * SecondCoM[2]) / Totalmass
    UpdatedCoM[3] = Totalmass
    return UpdatedCoM


def DivisionIndex(Particle, Center_Of_Box):
    """
    Takes the position vector of a particle and a position vector of a center. Calculates the index of where in the eight quadrants our particle is in octree division
    :param Particle: Position of the particle whose subdivision we want to find
    :param Center_Of_Box: Position around which we perform octree division
    :return: An integer from 0 to 7 indicating the division,
    Indices numbering is in the structure: 0=(-,-,-) -> 1=(+,-,-) -> 2=(-,+,-) -> 3=(+,+,-) -> 4=(-,-,+) ... where - denotes below the CoM position and + above for (x,y,z)
    """
    if Particle[0] < Center_Of_Box[0]:
        if Particle[1] < Center_Of_Box[1]:
            if Particle[2] < Center_Of_Box[2]:
                return 0
            else:
                return 4
        else:
            if Particle[2] < Center_Of_Box[2]:
                return 2
            else:
                return 6
    else:
        if Particle[1] < Center_Of_Box[1]:
            if Particle[2] < Center_Of_Box[2]:
                return 1
            else:
                return 5
        else:
            if Particle[2] < Center_Of_Box[2]:
                return 3
            else:
                return 7


def Box_Center_Displacer(DivIndex, Box_Center, Depth):
    """
    Takes a Division index, position and Depth and returns the new position of the subbox made when making the eight children bins and moving to the corresponding index
    :param Div_Index: Index which indicates in which direction the CoM must be updated
    :param Box_Center: Position that will be updated to the center of a subbox
    :param Depth: Indicates how deep this child division is which determines how much Center_Of_Box moves
    :return: The position of the center of the subbox with Index = Div_Index
    """
    if DivIndex == 0:
        return [Box_Center[0] - Length / 2 ** (Depth), Box_Center[1] - Length / 2 ** (Depth),
                Box_Center[2] - Length / 2 ** (Depth)]
    elif DivIndex == 1:
        return [Box_Center[0] + Length / 2 ** (Depth), Box_Center[1] - Length / 2 ** (Depth),
                Box_Center[2] - Length / 2 ** (Depth)]
    elif DivIndex == 2:
        return [Box_Center[0] - Length / 2 ** (Depth), Box_Center[1] + Length / 2 ** (Depth),
                Box_Center[2] - Length / 2 ** (Depth)]
    elif DivIndex == 3:
        return [Box_Center[0] + Length / 2 ** (Depth), Box_Center[1] + Length / 2 ** (Depth),
                Box_Center[2] - Length / 2 ** (Depth)]
    elif DivIndex == 4:
        return [Box_Center[0] - Length / 2 ** (Depth), Box_Center[1] - Length / 2 ** (Depth),
                Box_Center[2] + Length / 2 ** (Depth)]
    elif DivIndex == 5:
        return [Box_Center[0] + Length / 2 ** (Depth), Box_Center[1] - Length / 2 ** (Depth),
                Box_Center[2] + Length / 2 ** (Depth)]
    elif DivIndex == 6:
        return [Box_Center[0] - Length / 2 ** (Depth), Box_Center[1] + Length / 2 ** (Depth),
                Box_Center[2] + Length / 2 ** (Depth)]
    elif DivIndex == 7:
        return [Box_Center[0] + Length / 2 ** (Depth), Box_Center[1] + Length / 2 ** (Depth),
                Box_Center[2] + Length / 2 ** (Depth)]


def DeeperArray(Div_Index, Displaced_Center, Particle_Box):
    """
    This function is called when 2 particles need to go in a box that has not been split yet, this divides it into 8 subboxes and puts the original particle
    into it's new corresponding subbox and returns the array of subboxes
    :param Div_Index: Original index of the unsplit array
    :param Displaced_Center: Center of Box that needs to be used to find the new division index of the original particle
    :param Particle_Box: This is the unsplit array with the original particle
    :return: An array with 8 daughter arrays, one of which has the original particle in it, and a center of mass on it's 9th position
    """
    Index = DivisionIndex(Particle_Box[Div_Index], Displaced_Center)
    Array = [[] for i in range(9)]
    Array[8] = Particle_Box[Div_Index]
    Array[Index] = Particle_Box[Div_Index]
    return Array


def Array_updater(Daughter_Array):
    """
    This function allows the 9th element (the center of masses) of all arrays and their subarrays to be updated without rebuilding the tree structure
    :param Daughter_Array: Array whose Center of masses of itself and it's daughter arrays we want to update
    :return: Array with Center of masses updated to the current known positions
    """
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
    """
    Creates an array of nested arrays using the current particle positions and masses, the particle values get updated since they are saved as a reference,
    however the center of masses don't update (see Array_updater()) and when a particle leaves it's box this doesn't get updated either.
    This groups together nearby particles allowing for the usage of a tree code algorithm but it needs to be called multiple times during running to rebuild.
    :return: Array of nested arrays with particles in it and the center of mass and mass in position 9 (index 8) of all daughter arrays on the same depth
    """
    Boxstructure = [[], [], [], [], [], [], [], [], [0, 0, 0, 0]]
    for i in range(len(Particle)):
        Looping = True
        Depth = 0
        Center_Of_Box = [0, 0, 0]
        Particle_Box = Boxstructure
        for k in range(3):
            if abs(Particle[i][0][k]) > Length:
                Particle[i][0][k] = ((Particle[i][0][k] + Length) % (2 * Length) - Length)

        while Looping:

            Div_Index = DivisionIndex(Particle[i][0], Center_Of_Box)

            if len(Particle_Box[Div_Index]) > 1:  # This means it is not an empty array [], here we need to update our array correspondingly
                if isinstance(Particle_Box[Div_Index][0],
                              float):  # If True, these are positions and we need to make new daughter array
                    Displaced_Center = Box_Center_Displacer(Div_Index, Center_Of_Box, Depth + 1)
                    Array = DeeperArray(Div_Index, Displaced_Center, Particle_Box)
                    Particle_Box[Div_Index] = Array
                else:  # This is a daughter array so we update our CoM and Array and repeat to see where in this array our particle should go
                    Particle_Box[8] = CoMUpdater(Particle_Box[8], Particle[i][0])
                    Particle_Box = Particle_Box[Div_Index]
                    Center_Of_Box = Box_Center_Displacer(Div_Index, Center_Of_Box, Depth + 1)

                    Depth += 1
            else:  # It is empty so we just insert it
                Particle_Box[8] = CoMUpdater(Particle_Box[8], Particle[i][0])
                Particle_Box[Div_Index] = Particle[i][0]
                Looping = False
    return Boxstructure


def ForceParticles(Vector1, Vector2):
    """
    Calculates the gravitational force with a damped effect by using distance R = sqrt(r²+eps²) instead of r between 2 particles
    :param Vector1: Position of particle 1
    :param Vector2: Position of particle 2
    :return: Gravitational force vector on particle 1 due to particle 2
    """
    Distance = DistanceParticles(np.array(Vector1), np.array(Vector2[0:3]))

    return Vector2[3] * G * (np.array(Vector2[0:3]) - np.array(Vector1)) / (
            Distance ** 2 + Softening_Parameter ** 2) ** 1.5


def TotalForce(Boxstructure, Force, Length, Part_Index, Calculations):
    """
    Calculates the gravitational force with a damped effect by using distance R = sqrt(r²+eps²) instead of r on particle with index Part_Index and keeps track
    of the amount of force calculations performed
    :param Boxstructure: Tree structure with the particles and their center of masses and mass
    :param Force: Input force, this is added on top, used in this function when called recursively
    :param Length: Length of the box, when recursively calling we have to divide length by 2 every time we do
    :param Part_Index: Index of the particle in the Particles array
    :param Calculations: Amount of force calculations performed, used when called recursively
    :return: Gravitational force vector on particle 'Part_Index' and the amount of force calculations performed for this particle
    """
    if (Length / DistanceParticles(np.array(Particle[Part_Index][0][0:3]),
                                   np.array(Boxstructure[8][0:3]))) < Angle_Criterion:
        Force += ForceParticles(Particle[Part_Index][0][0:3], Boxstructure[8])
        Calculations += 1

    else:
        for j in range(8):
            if len(Boxstructure[j]) > 0:  # Is empty? we skip
                if isinstance(Boxstructure[j][0], float):  # Is number? We cannot go deeper and do force with this
                    if not Boxstructure[j] == Particle[Part_Index][0]:  # Except if this was our original particle
                        Force += ForceParticles(Particle[Part_Index][0][0:3], Boxstructure[j])
                        Calculations += 1
                else:  # We can and must go deeper in the TreeStructure
                    AddedForce, Calculations = TotalForce(Boxstructure[j], 0, Length / 2, Part_Index, Calculations)
                    Force += AddedForce
    return Force, Calculations





def VelocityUpdate(Timestep):
    """
    Calculates the new velocity of all particles by calculating the acceleration on all particles, also gives amount of force calculations for all particles to do this
    :param Timestep: Size of the timestep
    :return: The updated Particle array with new velocities and an array with the amount of force calculations for all particles
    """
    ForceCalculations = [() for i in range(len(Particle))]
    for i in range(len(Particle)):
        Force, ForceCalculations[i] = TotalForce(Boxstructure, 0, Length, i, 0)
        for k in range(3):
            Particle[i][1][k] += Timestep * Force[
                k]
    return Particle, ForceCalculations


def PositionUpdate(Timestep):
    """
    Calculates the new position of all particles by using the known velocities of the particles
    :param Timestep: Size of the timestep
    :return: The updated Particle array with new positions
    """
    for i in range(len(Particle)):
        for k in range(3):
            Particle[i][0][k] += Timestep * Particle[i][1][k]
    return Particle


Energy = []


def Kin_Energy():
    """
    Calculates the combined kinetic energy of all particles
    :return: Scalar value of the kinetic energy
    """
    E_Kin = 0
    for i in range(len(Particle)):
        E_Kin += 0.5 * Particle[i][0][3] * (Particle[i][1][0] ** 2 + Particle[i][1][1] ** 2 + Particle[i][1][2] ** 2)
    return E_Kin


def Pot_Energy():
    """
    Calculates the combined potential energy of all particles, note this scales as O(N²), don't use this for large simulations!
    :return: Scalar value of the potential energy
    """
    E_Pot = 0
    for i in range(len(Particle)):
        for j in range(len(Particle) - i - 1):
            E_Pot += -G * Particle[i][0][3] * Particle[i + j + 1][0][3] / DistanceParticles(
                np.array(Particle[i][0][0:3]), np.array(Particle[i + j + 1][0][0:3]))
    return E_Pot



Particle = []
Energy_Rate = 10
Energy_Tracking = False
G = 1
Max_R = 15
N1d = 3 * N
N1b = N
N_Galaxy1 = N1d + N1b
M1d = 1 * NFold
M1b = 0.333 * NFold
Softening_Parameter = 0.98 * (N1d + N1b) ** (-0.26)
galaxy1 = make_galaxy(N1d, M1d, N1b, M1b)
Particle = galaxy1
PlotParticle1 = np.array(galaxy1)




Boxstructure = [[], [], [], [], [], [], [], [], [0, 0, 0, 0]]  #This is not needed but to show how the treestructure is built and evolves when adding particles in the next line
# [...,[Part1,M1],...[COM1,M1]]  -> [...,[[],...,[Part1,M1],...[Part2,M2],...[COM2],[M2]],...[COM1,M1]]
Center_Of_Box = [0, 0, 0]
Particle_Box = Boxstructure
Looping = True
Angle_Criterion = 0.5 #In this code our value has to be multiplied by 2 to correspond with values in literature
Timestep_Size = 0.76










Boxstructure = TreeGenerator()
Energy = []
Energy.append([Kin_Energy(), Pot_Energy()])

ParticlePosHistory = [[] for k in range(Steps + 2)]
ParticlePosHistory[0] = [Particle[k][0][0:3] for k in range(len(Particle))]

Particle = VelocityUpdate(Timestep=Timestep_Size / 2)[0]
Particle = PositionUpdate(Timestep=Timestep_Size)
print('Start running')
ParticlePosHistory[1] = [Particle[k][0][0:3] for k in range(len(Particle))]
bar = tqdm(range(Steps))
ForceCalculations = [() for i in range(Steps)]
for t in bar:
    Particle, ForceCalculations[t] = VelocityUpdate(Timestep=Timestep_Size)
    Particle = PositionUpdate(Timestep=Timestep_Size)
    ParticlePosHistory[t + 2] = [Particle[k][0][0:3] for k in range(len(Particle))]


    if Energy_Tracking:
        if t % Energy_Rate == Energy_Rate - 1:
            Energy.append([Kin_Energy(), Pot_Energy()])

    #Boxstructure = Array_updater(Boxstructure)
    if t % 1 == 0: #This can be used to not remake tree structure every iteration by changing the 1.
        # The above line with Array_updater needs to be taken out of comment then however
        Boxstructure = TreeGenerator()

ParticlePosHistory = np.array(ParticlePosHistory)

fig, ax = plt.subplots()
ln1, = plt.plot([], [], 'go')


def init():
    ax.set_xlim(-Length, Length)
    ax.set_ylim(-Length, Length)
    ax.set_xlabel('x (unit=3.5kpc)')
    ax.set_ylabel('y (unit=3.5kpc)')
    ax.set_title('Time:' + str(Steps * 10) + 'Myr')
    ax.legend(["galaxy"], loc="upper right")


def update(q):
    ln1.set_data([ParticlePosHistory[int(q)][0:N1d + N1b, 0]], [ParticlePosHistory[int(q)][0:N1d + N1b, 1]])


ani = FuncAnimation(fig, update, frames=int(Steps), interval=20, init_func=init)
ani.save(folder + '\\' + 'animation.gif', writer='pillow')
plt.figure(2)
plt.scatter(PlotParticle1[:, 0, 0], PlotParticle1[:, 0, 1], marker='.')
plt.savefig(folder + '\\' + 'galaxy.png')
plt.figure(3)
plt.xlabel('x position (Unit = 3.5kpc)')
plt.ylabel('y position (Unit = 3.5kpc)')
plt.plot(ParticlePosHistory[:, :, 0], ParticlePosHistory[:, :, 1])
plt.savefig(folder + '\\' + 'galaxyTrajectory.png')
ForceCalculations = np.array(ForceCalculations)
plt.figure(4)
plt.xlabel('Time (Unit=10Myr)')
plt.ylabel('Number of force calculations')
plt.plot(range(Steps), ForceCalculations[:, :])
plt.savefig(folder + '\\' + 'ForceCalculation.png')

if Energy_Tracking:
    plt.figure(5)
    Energy = np.array(Energy)
    plt.plot(Energy[:, 0])
    plt.plot(Energy[:, 1])
    plt.plot(Energy[:, 0] + Energy[:, 1])
    plt.legend(["E_Kin", "E_Pot", "Total Energy"])
    plt.xlabel('Time (Unit=10 Myr)')
    plt.ylabel('Energy')
