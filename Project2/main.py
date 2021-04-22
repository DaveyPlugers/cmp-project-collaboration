"""
2D Ising Model code written by Zhen Xiang & Davey Plugers

Can be called in Correlation mode with a specified temperature, lattice size, equilibration time and Total run time.
Or it can be called in Thermodynamic mode with a specified temperature, lattice size, equilibration time, correlation time
and the amount of requested data points.

Afterwards in Correlation mode it will return a plot of the final spin alignment, 4 plots of the mean magnetisation of
different simulations and 4 estimates of the correlation time.
In Thermodynamic mode it will return a plot of the final spin alignment,
a plot of the mean magnetisation until equilibration time and an estimate of the mean and the standard deviation
magnetisation per spin, energy per spin, magnetic susceptibility per spin and specific heat per spin
_________________________________________________________

Functions available:

Random_Initialisation()

Total_Energy()

Energy_Difference(x_index, y_index)

Flip_Site(Delta_Energy)

Auto_Correlation(t_equilibrium, Spin_Magn)

Correlation_Time(Auto_Corr)

Magnetic_Suscept(Total_Magnetisations)

Specific_Heat(Total_Energy)

Standard_deviation(dataset, t_Max)
"""

import os
import math
import datetime
import argparse
import random
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

nowTime = datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S')
folder = os.getcwd() + '\\' + nowTime + '\\' + 'plot'

if not os.path.exists(folder):
    os.makedirs(folder)

folder1 = os.getcwd() + '\\' + nowTime + '\\' + 'data'

if not os.path.exists(folder1):
    os.makedirs(folder1)

parser = argparse.ArgumentParser()
parser.add_argument('--Temperature',
                    help='Value of system temperature, default: 1',
                    nargs=1)

parser.add_argument('--Lattice',
                    help='Width of the grid in spins sites, default: 50',nargs=1)
parser.add_argument('--Equilibrium_time',
                    help='Time duration until system is considered in equilibrium: 1000',nargs=1)
parser.add_argument('--Correlation_time',
                    help='Value of system correlation time, default: 1',
                    nargs=1)
parser.add_argument('--Total_Duration',
                    help='Time duration the total simulation will run in correlation mode, default: 3000',nargs=1)
parser.add_argument('--Data_Points',
                    help='Amount of independent blocks to calculate physical quantities: 20',nargs=1)

parser.add_argument('--Correlation_Mode',
                    help='Analysis system to get the correlation time',
                    action='store_true')

parser.add_argument('--Thermodynamic_Mode',
                    help='Analysis system to get the Thermodynamic properties',
                    action='store_true')

parser.add_argument('--No_Plot',
                    help='Allows plots to be turned off',
                    action='store_true')


args = parser.parse_args()

Correlation_Mode = args.Correlation_Mode
Thermodynamic_Mode = args.Thermodynamic_Mode
No_Plot = args.No_Plot

if args.Temperature is None:
    Temperature = 1.
else:
    Temperature = float(args.Temperature[0])

if args.Lattice is None:
    Lattice = 50
else:
    Lattice = int(args.Lattice[0])

if args.Equilibrium_time is None:
    t_equilibrium = 1000
else:
    t_equilibrium = int(args.Equilibrium_time[0])

if args.Correlation_time is None:
    t_correlation = 1
else:
    t_correlation = int(args.Correlation_time[0])

if args.Total_Duration is None:
    Repetitions = Lattice**2 * 3000
else:
    Repetitions = Lattice**2 *int(args.Total_Duration[0])

if args.Data_Points is None:
    Amount_Blocks = 20
else:
    Amount_Blocks = int((args.Data_Points[0]))



J = 1
KB = 1
LatticeSquared = Lattice ** 2


f = open(folder1 + '\\' + 'output.txt', 'w')
f.write('Initial condition: Temperature = ' + str(Temperature))
f.close()


def Random_Initialisation():
    """
    Can be called to return an initialisation with spin up or down randomly distributed
    :return: np.array consisting of -1 or 1 representing the spins
    """
    return np.array([[int(2 * (random.randint(0, 1) - 0.5)) for x in range(Lattice)] for y in range(Lattice)])


Spins_System = Random_Initialisation()


def Total_Energy():
    """
    Can be called to calculate the total energy in the system
    :return: Scalar value of the total energy of the system
    """
    Energy = 0
    for x in range(Lattice):
        for y in range(Lattice):
            Energy += -0.5 * J * Spins_System[y, x] * (
                    Spins_System[y, (x + 1) % Lattice] + Spins_System[y, (x - 1) % Lattice] + Spins_System[
                (y + 1) % Lattice, x] + Spins_System[(y - 1) % Lattice, x])
    return Energy


Energy = Total_Energy()


def Energy_Difference(x_index, y_index):
    """
    Takes index of a spin, returns the energy difference if we would flip this spin
    :param x_index: x position of the spin
    :param y_index: y position of the spin
    :return: Scalar energy difference if the spin is flipped
    """
    Energy_Diff = 2 * J * Spins_System[y_index, x_index] * (
            Spins_System[y_index, (x_index + 1) % Lattice] + Spins_System[y_index, (x_index - 1) % Lattice] +
            Spins_System[(y_index + 1) % Lattice, x_index] + Spins_System[(y_index - 1) % Lattice, x_index])
    return Energy_Diff


Probabilities = [math.exp(-4 / (KB * Temperature)), math.exp(-8 / (KB * Temperature))]


def Flip_Site(Delta_Energy):
    """
    Takes energy difference, calculates whether the spin flips with probability e^(-Delta_Energy/(kb*T))
    :param Delta_Energy: Energy difference for if it is flipped, this value can only be 4 or 8 since 2d grid and
    values =< 0 always flip
    :return: Boolean for if it should flip
    """
    if Delta_Energy == 4:
        Boolean = random.random() < Probabilities[0]
    else:
        Boolean = random.random() < Probabilities[1]
    return Boolean


def Auto_Correlation(t_equilibrium, Spin_Magn):
    """
    Takes a time value after which the system is in equilibrium and an array of mean magn. per spin
    Returns an array of values of the autocorrelation function of the magnetisation
    :param t_equilibrium: Time value that indicates when the system is in equilibrium
    :param Spin_Magn: Array of mean magn. per spin
    :return: Array that corresponds with the discrete autocorrelation function X(t)
    """
    t_max = len(Spin_Magn) - t_equilibrium
    Auto_Corr = np.zeros(t_max)
    for t in range(t_max):
        Auto_Corr[t] = 1 / (t_max - t) * (np.sum(Spin_Magn[t_equilibrium:t_equilibrium + t_max - t] *
                                                 Spin_Magn[t_equilibrium + t: t_equilibrium + t_max])) - (
                               1 / (t_max - t) ** 2 *
                               np.sum(Spin_Magn[t_equilibrium: t_equilibrium + t_max - t]) *
                               np.sum(Spin_Magn[t_equilibrium + t: t_equilibrium + t_max]))

    return Auto_Corr


def Correlation_Time(Auto_Corr):
    """
    Takes an array of autocorrelation values, returns the correlation time
    :param Auto_Corr: Array that corresponds with the discrete autocorrelation function X(t)
    :return: Scalar value of the correlation time
    """
    Tau = 0
    for t in range(len(Auto_Corr)):
        if Auto_Corr[t] > 0:
            Tau += Auto_Corr[t] / Auto_Corr[0]
        else:
            return Tau
    return Tau


def Magnetic_Suscept(Total_Magnetisations):
    """
    Takes an array of total magnetisation, returns the Magnetic Susceptibility calculated from these values
    :param Total_Magnetisations: Array that has the total magnetisation of the system at different times (in our case 8
    values)
    :return: Scalar value of the approximated Magnetic Susceptibility
    """
    return 1 / (KB * Temperature * LatticeSquared) * (
            np.mean(Total_Magnetisations ** 2) - np.mean(Total_Magnetisations) ** 2)


def Specific_Heat(Total_Energy):
    """
    Takes an array of total energy, returns the Specific Heat calculated from these values
    :param Total_Energy: Array that has the total energy of the system at different times (in our case 8 values)
    :return: Scalar value of the approximated Specific Heat
    """
    return 1 / (KB * LatticeSquared * Temperature ** 2) * (np.mean(Total_Energy ** 2) - np.mean(Total_Energy) ** 2)


def Standard_deviation(dataset, t_Max):
    """
    Calculate the standard deviation for Thermodynamic properties
    :param dataset: array of the Thermodynamic properties
    :param t_Max: total number of sweeps after equilibrium
    :return: Standard deviation value for the Thermodynamic properties
    """
    return ((2 * t_correlation / t_Max) * (np.mean(dataset ** 2) - np.mean(dataset) ** 2)) ** 0.5





if Correlation_Mode:
    bar = tqdm(range(Repetitions))  # Running full simulation
    Amount_Initialisations = 4
    Multiple_Spin_Magn = [0, 0, 0, 0]
else:
    bar = tqdm(range(
        t_equilibrium * LatticeSquared))  # Running untill equilibrium then start a continuation for physical values later
    Amount_Initialisations = 1

for z in range(Amount_Initialisations):
    Spin_Magn = np.array([])
    for a in bar:
        x_index = random.randint(0, Lattice - 1)
        y_index = random.randint(0, Lattice - 1)

        Delta_Energy = Energy_Difference(x_index, y_index)
        if Delta_Energy <= 0:
            Spins_System[y_index, x_index] = -Spins_System[y_index, x_index]
            Energy += Delta_Energy
        elif Flip_Site(Delta_Energy):
            Spins_System[y_index, x_index] = -Spins_System[y_index, x_index]
            Energy += Delta_Energy
        if a % LatticeSquared == 0:
            Spin_Magn = np.append(Spin_Magn, [np.sum(Spins_System) / Lattice ** 2])

    if Correlation_Mode:
        Multiple_Spin_Magn[z] = Spin_Magn
        if not z == 3:
            Spins_System = Random_Initialisation()
            bar = tqdm(range(Repetitions))

if Thermodynamic_Mode:
    Magn_PSpin = np.array([])
    Energy_PSpin = np.array([])
    Estimated_Magnet_Suscept = np.array([])
    Estimated_Spec_Heat = np.array([])

    Total_Energy = np.array([])  # Can use Energy_PSpin too but this is easier
    Total_Magnetisations = np.array([])


    t_max = Amount_Blocks * 16 * t_correlation
    bar = tqdm(range(8 * Amount_Blocks))
    Measurement_Rate = math.ceil(2 * t_correlation * LatticeSquared)
    for i in bar:  # Create blocks of length 16*t_corr each consisting of 8 individual measurements every 2*t_corr
        for a in range(Measurement_Rate):
            x_index = random.randint(0, Lattice - 1)
            y_index = random.randint(0, Lattice - 1)

            Delta_Energy = Energy_Difference(x_index, y_index)
            if Delta_Energy <= 0:
                Spins_System[y_index, x_index] = -Spins_System[y_index, x_index]
                Energy += Delta_Energy
            elif Flip_Site(Delta_Energy):
                Spins_System[y_index, x_index] = -Spins_System[y_index, x_index]
                Energy += Delta_Energy
        Sum_Of_Spins = np.sum(Spins_System)
        Magn_PSpin = np.append(Magn_PSpin, abs(Sum_Of_Spins) / Lattice ** 2)
        Energy_PSpin = np.append(Energy_PSpin, Energy / Lattice ** 2)
        Total_Energy = np.append(Total_Energy, Energy)
        Total_Magnetisations = np.append(Total_Magnetisations, Sum_Of_Spins)
        if (i + 1) % 8 == 0:
            Estimated_Magnet_Suscept = np.append(Estimated_Magnet_Suscept, Magnetic_Suscept(Total_Magnetisations))
            Estimated_Spec_Heat = np.append(Estimated_Spec_Heat, Specific_Heat(Total_Energy))
            Total_Energy = np.array([])
            Total_Magnetisations = np.array([])
    f = open(folder1 + '\\' + 'output.txt', 'a')
    f.write('\nFor magnetization per spin, the mean is ' + str(np.mean(Magn_PSpin))
            + '\nthe standard deviation is ' + str(Standard_deviation(Magn_PSpin, t_max)))
    f.write('\nFor energy per spin, the mean is ' + str(np.mean(Energy_PSpin))
            + '\nthe standard deviation is ' + str(Standard_deviation(Energy_PSpin, t_max)))
    f.write('\nFor magnetic susceptibility per spin, the mean is ' + str(np.mean(Estimated_Magnet_Suscept))
            + '\nthe standard deviation is ' + str(Standard_deviation(Estimated_Magnet_Suscept, t_max/8)))
    f.write('\nFor specific heat per spin, the mean is ' + str(np.mean(Estimated_Spec_Heat))
            + '\nthe standard deviation is ' + str(Standard_deviation(Estimated_Spec_Heat, t_max/8)))
    f.close()

if Correlation_Mode:
    Calculated_Correlation_Tau = [0,0,0,0]
    if t_equilibrium < Repetitions / LatticeSquared:
        for i in range(4):
            Auto_Corr = Auto_Correlation(t_equilibrium, Multiple_Spin_Magn[i])
            Calculated_Correlation_Tau[i] = ", " + str(Correlation_Time(Auto_Corr))
        f = open(folder1 + '\\' + 'output.txt', 'a')
        for j in range(4):
            Values = ''.join(Calculated_Correlation_Tau)
        f.write('\nThe correlation times of the system are: ' + Values)
        f.close()
    if not No_Plot:
        plot2 = plt.figure(2)
        Colours = ['r', 'b', 'g', 'k']
        for z in range(4):
            plt.plot(Multiple_Spin_Magn[z], Colours[z])
        plt.savefig(folder + '\\' + 'magnetization.png')


if not No_Plot:
    if Thermodynamic_Mode:
        plt.figure(2)
        plt.plot(Spin_Magn)
        plt.savefig(folder + '\\' + 'magnetization.png')

    plt.figure(1)
    imgplot = plt.imshow(Spins_System)
    plt.colorbar()
    plt.savefig(folder + '\\' + 'spin_alignment.png')

