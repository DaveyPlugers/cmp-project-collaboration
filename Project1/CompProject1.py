"""
# Framework of Allposition array: Allposition[i][j][k] =>
# i gives timestep [Timestep0,Timestep1,Timestep2,...] Will be appended like that where =>
# j gives particle number Timestep0 = [Particle1,Particle2] where =>
# k=0 gives position vector, k=1 gives velocity vector Particlej = [posjvector,veljvector]
"""

# import part
import os
import math
import random
import datetime
import argparse
import numpy as np
import matplotlib.pyplot as plt
from numpy import ndarray
from tqdm import tqdm
from matplotlib.animation import FuncAnimation
from scipy.interpolate import make_interp_spline

# creat folder
nowTime = datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S')
folder = os.getcwd() + '\\' + nowTime + '\\' + 'plot'

if not os.path.exists(folder):
    os.makedirs(folder)

# parser part
parser = argparse.ArgumentParser()
parser.add_argument('--Density',
                    help='Value of box density, default: 1',
                    nargs=1)
parser.add_argument('--Temperature',
                    help='Value of system temperature, default: 1',
                    nargs=1)
args = parser.parse_args()

if args.Density is None:
    Density = 1.
else:
    Density = float(args.Density[0])

if args.Temperature is None:
    Temperature = 1.
else:
    Temperature = float(args.Temperature[0])

# set initial condition
# Temperature = 3.0  # Make this an input variable later
# Density = 0.3  # In units, atoms/sigma**3, make input variable later
Mass = 1  # Mass in atomic mass units
BoxSize = 3 * (4 * Mass / Density) ** (1 / 3)  # Times 4 since 4 particles per cube
print("Boxsize = " + str(BoxSize))
TimeSteps = 500
TimeStepLength = 0.001

HistBins = 20
HistStart = 100
HistTimes = np.linspace(HistStart, TimeSteps, num=int(((TimeSteps - HistStart) / 50)) + 1)
RescaleTimes = np.linspace(100, TimeSteps, num=int(TimeSteps / 100))
# RescaleTimes = [15, 50, 80, 150]
# HistTimes = [HistStart + (1 + x) * 50 for x in range(int((TimeSteps - HistStart) / 50))]
PressureTimes = [200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950]

RandomInitialisation = False
CreatePlots = True

Particles = []
if RandomInitialisation:
    ParticleAmount = 8
    Dimension = 2
    ParticleVelocity = 0.001
    for l in range(ParticleAmount):
        Partic = [[BoxSize * random.random() for d in range(Dimension)],
                  [2 * ParticleVelocity * (random.random() - 0.5) for d in range(Dimension)]]
        Particles.append(Partic)
else:
    Dimension = 3
    ParticleAmount = 108
    ParticleNumber = 0
    VelVariance = Temperature ** (1 / 2)
    Velocities = np.random.normal(0, VelVariance, 324)
    for i in range(6):
        for j in range(6):
            for k in range(3):
                Partic = [[k * BoxSize / 3 + (j + i) % 2 * BoxSize / 6, j * BoxSize / 6, i * BoxSize / 6],
                          # Modulo operators create FCC lattice structure
                          [Velocities[3 * ParticleNumber], Velocities[3 * ParticleNumber + 1],
                           Velocities[3 * ParticleNumber + 2]]]
                ParticleNumber += 1
                Particles.append(Partic)
Particles = np.array(Particles).astype('float64')
AllPositions = [[] for x in range(TimeSteps + 1)]  # Premaking is faster than appending in Matlab, not sure about Python
AllPositions[0] = Particles.copy()


# function part
def DistancePoints(Vector1, Vector2):
    """
    Takes 2 position vectors and calculates distance between them
    :param Vector1: position vector for particle one
    :param Vector2: position vector for particle two
    :return: distance between two particles
    """
    Distance = sum([x * x for x in ((Vector1 - Vector2 + BoxSize / 2) % BoxSize - BoxSize / 2)]) ** 0.5
    return Distance


def VecPairLJForce(Vector1, Vector2):
    """
    Takes 2 position vectors and returns the Lennard-Jones force on Vector1 due to Vector2
    :param Vector1: position vector for particle one
    :param Vector2: position vector for particle two
    :return: Lennard-Jones force vector between two particles
    """
    Dist = DistancePoints(Vector1, Vector2)
    Constant = 24 * (2 / Dist ** 14 - 1 / Dist ** 8)  # Extra 1/r for normalisation Vector1 - Vector2 in next line
    Force = [Constant * x for x in ((Vector1 - Vector2 + BoxSize / 2) % BoxSize - BoxSize / 2)]
    return Force


def TotalForce(j, tstep):
    """
    This calculates the total force vector experienced by particle j
    :param j: index of a particle
    :param tstep: time step with a special case for tstep==-1 for the updated values that aren't saved yet
    :return: Total force on particle j
    """

    TotalForce = np.array([0 for l in range(Dimension)], dtype=np.float64)
    if tstep == -1:  # This is a workaround to get the updated position and thus updated force
        for i in range(ParticleAmount):
            if i != j:
                TotalForce += VecPairLJForce(Particles[j, 0], Particles[i, 0])
    else:
        for i in range(ParticleAmount):
            if i != j:
                TotalForce += VecPairLJForce(AllPositions[tstep][j, 0], AllPositions[tstep][i][
                    0])  # We use the Allpositions instead of Particles since we need the old positions

    return TotalForce


Epot = []
Ekin = []
Etot = []


def Energy():
    """
    Can be called to save the potential,kinetic and total energy of the system at the current timestep
    """
    pot = 0
    kin = 0
    for j in range(ParticleAmount):
        pot1 = 0
        kin += (0.5 * Mass) * (sum([x * x for x in Particles[j, 1]]))
        for k in range(ParticleAmount):
            if k != j:
                r = DistancePoints(Particles[j, 0], Particles[k, 0])
                pot1 += 4 * ((1 / r) ** 12 - (1 / r) ** 6)
        pot += pot1 / 2  # *0.5 to prevent double counting
    Epot.append(pot)
    Ekin.append(kin)
    Etot.append(pot + kin)


def Rescale():
    """
    Can be called to rescale the kinetic energy of the system
    :return: Scale factor for the velocities
    """
    VelSquaredSummed = 0
    for j in range(ParticleAmount):
        VelSquaredSummed += sum([x * x for x in Particles[j, 1]])
    Lambda = (((ParticleAmount - 1) * 3 * Temperature) / VelSquaredSummed) ** (1 / 2)
    return Lambda


def Histogram(Bins):
    """
    Can be called to calculate a histogram giving number of particle pairs at certain distances
    :param Bins: Amount of intervals of the histogram
    :return: An array of Histogram values
    """
    Histo = np.histogram([BoxSize], bins=[BoxSize / (2 * Bins) * x for x in range(Bins + 1)])[0]
    for i in range(ParticleAmount):
        Distances = [DistancePoints(Particles[1 + i + j, 0], Particles[i, 0]) for j in range(ParticleAmount - i - 1)]
        Histo += np.histogram(Distances, bins=[BoxSize / (2 * Bins) * x for x in range(Bins + 1)])[0]
    # Histo = 2*Histo
    return Histo


DistTimesForce = []


def SumDistTimesForce():
    """
    Can be called to save the sum of the Distance times the Force of all particle pairs, this will later be averaged to calculate the pressure
    """
    DistTimesForce = 0
    for i in range(ParticleAmount):
        for j in range(ParticleAmount - i - 1):
            Dist = DistancePoints(Particles[1 + i + j, 0], Particles[i, 0])
            DistTimesForce += (12 * (-2 / Dist ** 12 + 1 / Dist ** 6))
    return DistTimesForce


def Pressure():
    print(DistTimesForce)
    P = Temperature * Density - (Density * np.average(np.array(DistTimesForce))) / (3 * ParticleAmount)
    return P


def PairCorrelation():
    PairCorrel = [(2 * BoxSize ** 3 * Histo[i]) / (
            ParticleAmount * (ParticleAmount - 1) * 4 * math.pi * ((i + 1) * BinSize) ** 2 * BinSize) for i in
                  range(len(Histo))]
    print("This is PairCorrel" + str(PairCorrel))
    return PairCorrel

ForcePrevious = [np.array(TotalForce(k,0)) for k in range(ParticleAmount)]
# running simulation
bar = tqdm(range(TimeSteps))
for tstep in bar:
    bar.set_description(f"running simulation")
    # print("timestep = " + str(tstep))
    for j in range(ParticleAmount):
        Particles[j, 0] = Particles[j, 0] + Particles[j, 1] * TimeStepLength + TimeStepLength ** 2 / (
                2 * Mass) * ForcePrevious[j]
        Particles[j, 0] = Particles[j, 0] % BoxSize
    for j in range(ParticleAmount):
        NewForce = TotalForce(j, -1)
        Particles[j, 1] = Particles[j, 1] + TimeStepLength / (2 * Mass) * ((NewForce + ForcePrevious[j])
        ForcePrevious[j] = NewForce
                                                                           
    if tstep % 10 == 0:
        Energy()

    if tstep in RescaleTimes:
        ScaleFactor = Rescale()
        print("Scalefactor= " + str(ScaleFactor))
        for k in range(ParticleAmount):
            Particles[k, 1] = ScaleFactor * Particles[k, 1]

    if tstep == HistStart:
        Histo = Histogram(HistBins)
    elif tstep in HistTimes:
        Histo += Histogram(HistBins)
        print(Histo)

    if tstep in PressureTimes:
        DistTimesForce.append(SumDistTimesForce())

    AllPositions[tstep + 1] = (Particles.copy())

if TimeSteps > HistStart:
    Histo = Histo * (1 / (len(HistTimes) - 1))
    BinSize = BoxSize / (2 * HistBins)
    g = PairCorrelation()
    print("(r+Binsize)^2 instead " + str(PairCorrelation()))

    # PairCorrelation = [(2*BoxSize**3*Histo[i+1]) / (ParticleAmount * (ParticleAmount - 1) * 4 * math.pi * ((i + 1) * BinSize) ** 2 * BinSize) for i in range(len(Histo) - 1)]
    # print("Skipped first one" + str(PairCorrelation))

print("Pressure of the system = " + str(Pressure()))

# print(Allpositions[0:5]) #Used for checking mistakes in simulation, remove later

# make plots
if CreatePlots:
    plot1 = plt.figure(1)
    plt.plot(range(0, TimeSteps, 10), Epot)
    plt.plot(range(0, TimeSteps, 10), Ekin)
    plt.plot(range(0, TimeSteps, 10), Etot)
    plt.xlabel('timestep')
    plt.ylabel('dimensionless energy')
    plt.legend(["Epot", "Ekin", "Etot"])
    plt.title('energy')
    plt.savefig(folder + '\\' + 'energy.png')

    plot2 = plt.figure(2)
    xnew = np.linspace(0, HistBins - 1, 300)
    g_smooth = make_interp_spline(range(HistBins), g)(xnew)
    plt.plot(xnew, g_smooth)
    plt.xlabel('r')
    plt.ylabel('g(r)')
    plt.title('pair correlation function')
    #plt.legend()
    plt.savefig(folder + '\\' + 'pair_correlation.png')
    # plt.show()

    fig, ax = plt.subplots()
    ln1, = plt.plot([], [], 'ro')


    def init():
        ax.set_xlim(0, BoxSize)
        ax.set_ylim(0, BoxSize)


    def update(q):
        ln1.set_data([AllPositions[int(10 * q)][:, 0, 0]], [AllPositions[int(10 * q)][:, 0, 1]])


    ani = FuncAnimation(fig, update, frames=int(TimeSteps / 10 + 1), interval=10, init_func=init)
    ani.save(folder + '\\' + 'animation.gif', writer='pillow')
    #plt.show()
