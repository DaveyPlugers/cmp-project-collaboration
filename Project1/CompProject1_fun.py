"""


"""

#Framework of Allposition array: Allposition[i][j][k] =>
#i gives timestep [Timestep0,Timestep1,Timestep2,...] Will be appended like that where =>
#j gives particle number Timestep0 = [Particle1,Particle2] where =>
#k=0 gives position vector, k=1 gives velocity vector Particlej = [posjvector,veljvector]

import math
import numpy as np
import random
import argparse
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from Force import DistancePoints
from Force import VecPairLJForce
from Force import TotalForce

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


#Temperature = 0.5  #Make this an input variable later
#Density= 1.2 #In units, atoms/sigma**3, make input variable later
Mass=1 #Mass in atomic mass units
BoxSize = 3 * (4 * Mass / Density) ** (1 / 3)  #Times 4 since 4 particles per cube
print("Boxsize = " + str(BoxSize))
TimeSteps = 400
TimeStepLength = 0.001

HistBins = 20
HistStart = 100
HistTimes = [HistStart + (1+x) * 50 for x in range(int((TimeSteps - HistStart) / 50))]
RescaleTimes = [15,50,80,150]
PressureTimes = [20, 80, 200, 250, 300,350]


RandomInitialisation = False
CreatePlots = True

Particles = []
if RandomInitialisation:
    ParticleAmount = 8
    Dimension = 2
    Particlevelocity = 0.001
    for l in range(ParticleAmount):

        Partic = [[BoxSize * random.random() for d in range(Dimension)], [2 * Particlevelocity * (random.random() - 0.5) for d in range(Dimension)]]
        Particles.append(Partic)
else:
    Dimension=3
    ParticleAmount = 108
    ParticleNumber = 0
    VelVariance = (Temperature) ** (1 / 2)
    Velocities = np.random.normal(0, VelVariance, 324)
    for i in range(6):
        for j in range(6):
            for k in range(3):
                Partic = [[k * BoxSize / 3 + (j + i) % 2 * BoxSize / 6, j * BoxSize / 6, i * BoxSize / 6], #Modulo operators create FCC lattice structure
                          [Velocities[3 * ParticleNumber], Velocities[3 * ParticleNumber + 1], Velocities[3 * ParticleNumber + 2]]]
                ParticleNumber += 1
                Particles.append(Partic)
Particles = np.array(Particles).astype('float64')
AllPositions = [[] for x in range(TimeSteps + 1)]  #Premaking is faster than appending in Matlab, not sure about Python
AllPositions[0] = Particles.copy()




Epot = []
Ekin = []
Etot = []
def Energy():
    '''
    Can be called to save the potential,kinetic and total energy of the system at the current timestep
    '''
    pot=0
    kin=0
    for j in range(ParticleAmount):
        pot1 = 0
        kin += (0.5 * Mass) * (sum([x * x for x in Particles[j, 1]]))
        for k in range(ParticleAmount):
            if k !=j:
                r = DistancePoints(Particles[j,0],Particles[k,0], BoxSize)
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
    for j in range(ParticleAmount):
        VelSquaredSummed += sum([x*x for x in Particles[j,1]])
    Lambda = (((ParticleAmount - 1) * 3 * Temperature) /VelSquaredSummed) ** (1 / 2)
    return Lambda

def Histogram(Bins):
    '''
    Can be called to calculate a histogram giving number of particle pairs at certain distances
    :param Bins: Amount of intervals of the histogram
    :return: An array of Histogram values
    '''
    Histo = np.histogram([BoxSize], bins=[BoxSize / (2 * Bins) * x for x in range(Bins + 1)])[0]
    for i in range(ParticleAmount):
        Distances = [DistancePoints(Particles[1+i+j,0],Particles[i,0]) for j in range(ParticleAmount - i - 1)]
        Histo += np.histogram(Distances, bins=[BoxSize/(2*Bins)*x for x in range(Bins+1)])[0]
    #Histo = 2*Histo
    return Histo


def Pressure():
    '''
    Can be called to save the pressure of the system at the current timestep
    '''
    DistTimesForce = []
    for i in range(ParticleAmount):
        for j in range(ParticleAmount - i - 1):
            Dist = DistancePoints(Particles[1+i+j,0],Particles[i,0])
            DistTimesForce.append(12*(-2/Dist**12 + 1/Dist**6))
    Pressure = 1 - (np.average(np.array(DistTimesForce)) / (3 * ParticleAmount * Temperature))
    print(np.average(np.array(DistTimesForce)))
    print(Dist)
    return Pressure






for tstep in range(TimeSteps):
    print("timestep = " + str(tstep))
    for j in range(ParticleAmount):
        Particles[j,0] = Particles[j,0] + Particles[j,1] * TimeStepLength + TimeStepLength ** 2 / (2 * Mass) * TotalForce(j, tstep, Dimension, ParticleAmount, AllPositions, BoxSize, Particles)
        Particles[j,0] = Particles[j,0] %BoxSize
    for j in range(ParticleAmount):
        Particles[j,1] = Particles[j,1] + TimeStepLength / (2 * Mass) * (TotalForce(j, -1, Dimension, ParticleAmount, AllPositions, BoxSize, Particles) + TotalForce(j, tstep, Dimension, ParticleAmount, AllPositions, BoxSize, Particles))


    if tstep%10==0:
        Energy()

    if tstep in RescaleTimes:
        ScaleFactor = Rescale()
        print("Scalefactor= " + str(ScaleFactor))
        for k in range(ParticleAmount):
            Particles[k,1] = ScaleFactor*Particles[k,1]

    if tstep ==HistStart:
        Histo = Histogram(HistBins)
    elif tstep in HistTimes:
        Histo += Histogram(HistBins)
        print(Histo)

    if tstep in PressureTimes:
        print("Pressure = " + str(Pressure()))

    AllPositions[tstep + 1] = (Particles.copy())


Histo = Histo*(1/(len(HistTimes)))
print(Histo)

BinSize = BoxSize/(2 * HistBins)
PairCorrelationDifferent = [(2*BoxSize**3*Histo[i]) / (ParticleAmount * (ParticleAmount - 1) * 4 * math.pi * ((i + 1) * BinSize) ** 2 * BinSize) for i in range(len(Histo))]
print("(r+Binsize)^2 instead " + str(PairCorrelationDifferent))

PairCorrelation = [(2*BoxSize**3*Histo[i+1]) / (ParticleAmount * (ParticleAmount - 1) * 4 * math.pi * ((i + 1) * BinSize) ** 2 * BinSize) for i in range(len(Histo) - 1)]
print("Skipped first one" + str(PairCorrelation))



#print(Allpositions[0:5]) #Used for checking mistakes in simulation, remove later

if CreatePlots:
    plot2 = figure(2)
    plot(range(0, TimeSteps, 10), Epot)
    plot(range(0, TimeSteps, 10), Ekin)
    plot(range(0, TimeSteps, 10), Etot)
    legend(["Epot", "Ekin", "Etot"])


    fig, ax = subplots()
    ln1, = plot([], [], 'ro')
    def init():
        ax.set_xlim(0, BoxSize)
        ax.set_ylim(0, BoxSize)
    def update(q):
        ln1.set_data([AllPositions[int(10 * q)][:, 0, 0]], [AllPositions[int(10 * q)][:, 0, 1]])
    ani = FuncAnimation(fig, update, frames=int((TimeSteps) / 10 + 1), interval=10, init_func=init)
    show()

