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

parser.add_argument('--Thermodynamics_Mode',
                    help='Analysis system to get the thermodynamics properties',
                    action='store_true')

args = parser.parse_args()

Correlation_Mode = args.Correlation_Mode
Thermodynamics_Mode = args.Thermodynamics_Mode

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
# Temperature = 1.5
Lattice = 50
LatticeSquared = Lattice ** 2
Repetitions = LatticeSquared * 5000
# Time in units of steps per site
t_equilibrium = 3000
# t_correlation = 5.11
Amount_Blocks = 50
# Correlation_Mode = True
# Thermodynamics_Mode = True

f = open(folder1 + '\\' + 'output.txt', 'w')
f.write('Initial condition: Temperature = ' + str(Temperature))
f.close()


def Random_Initialisation():
    """
    Can be called to return an initialisation with spin up or down randomly distributed
    :return: np.array consisting of -1 or 1 representing the spins
    """
    return np.array([[int(2 * (random.randint(0, 1) - 0.5)) for x in range(Lattice)] for y in range(Lattice)])


Spin_Array = Random_Initialisation()


def Energy_Calculation():
    """
    Can be called to calculate the total energy in the system
    :return: Scalar value of the total energy of the system
    """
    Energy = 0
    for x in range(Lattice):
        for y in range(Lattice):
            Energy += -0.5 * J * Spin_Array[y, x] * (
                    Spin_Array[y, (x + 1) % Lattice] + Spin_Array[y, (x - 1) % Lattice] + Spin_Array[
                (y + 1) % Lattice, x] + Spin_Array[(y - 1) % Lattice, x])
    return Energy


Energy = Energy_Calculation()


def Energy_Difference(x_index, y_index):
    """
    Takes index of a spin, returns the energy difference if we would flip this spin
    :param x_index: x position of the spin
    :param y_index: y position of the spin
    :return: Scalar energy difference if the spin is flipped
    """
    Energy_Diff = 2 * J * Spin_Array[y_index, x_index] * (
            Spin_Array[y_index, (x_index + 1) % Lattice] + Spin_Array[y_index, (x_index - 1) % Lattice] +
            Spin_Array[(y_index + 1) % Lattice, x_index] + Spin_Array[(y_index - 1) % Lattice, x_index])
    return Energy_Diff


Probabilities = [math.exp(-4 / (KB * Temperature)), math.exp(-8 / (KB * Temperature))]


def Flip_Calculator(Delta_Energy):
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


def Spin_Magnetization(Spin_Magn):
    """
    Takes the array of Spin Magnetisations, returns a version appended with the new value of mean magn. per spin
    :param Spin_Magn: Array of mean magn. per spin
    :return: Same array with additionally the current mean magn. per spin
    """
    return np.append(Spin_Magn, [np.sum(Spin_Array) / Lattice ** 2])


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


def Magnetic_Suscept(Magnetisation_Array):
    """
    Takes an array of total magnetisation, returns the Magnetic Susceptibility calculated from these values
    :param Magnetisation_Array: Array that has the total magnetisation of the system at different times (in our case 8
    values)
    :return: Scalar value of the approximated Magnetic Susceptibility
    """
    return 1 / (KB * Temperature * LatticeSquared) * (
            np.mean(Magnetisation_Array ** 2) - np.mean(Magnetisation_Array) ** 2)


def Specific_Heat(Energy_Array):
    """
    Takes an array of total energy, returns the Specific Heat calculated from these values
    :param Energy_Array: Array that has the total energy of the system at different times (in our case 8 values)
    :return: Scalar value of the approximated Specific Heat
    """
    return 1 / (KB * LatticeSquared * Temperature ** 2) * (np.mean(Energy_Array ** 2) - np.mean(Energy_Array) ** 2)


def Standard_deviation(dataset, t_Max):
    """
    Calculate the standard deviation for thermodynamics properties
    :param dataset: array of the thermodynamics properties
    :param t_Max: total number of sweeps after equilibrium
    :return: Standard deviation value for the thermodynamics properties
    """
    return ((2 * t_correlation / t_Max) * (np.mean(dataset ** 2) - np.mean(dataset) ** 2)) ** 0.5


# New structureplan: Have 2 different simulation possibilities:
# 1) Run multiple simulations on the same temperature, different initials, use these to calculate equilibrium and correlation time
# 2) Simulation where you give t_equi and t_corr, then make it run until t_equi, then start dividing in k blocks of 16t_corr
# Every 2t_corr we calculate quantities: {m, e,X, C}, we average X and C over the next 8 measurements
# And then use multiple blocks of these 16t_corr to do error analysis, for m and e we just do regular average with STD given


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
            Spin_Array[y_index, x_index] = -Spin_Array[y_index, x_index]
            Energy += Delta_Energy
        elif Flip_Calculator(Delta_Energy):
            Spin_Array[y_index, x_index] = -Spin_Array[y_index, x_index]
            Energy += Delta_Energy
        if a % LatticeSquared == 0:
            Spin_Magn = Spin_Magnetization(Spin_Magn)

    if Correlation_Mode:
        Multiple_Spin_Magn[z] = Spin_Magn
        if not z == 3:
            Spin_Array = Random_Initialisation()
            bar = tqdm(range(Repetitions))

if Thermodynamics_Mode:
    Magn_PSpin = np.array([])
    Energy_PSpin = np.array([])
    Magnet_Suscept_Array = np.array([])
    Spec_Heat_Array = np.array([])

    Energy_Array = np.array([])  # Can use Energy_PSpin too but this is easier
    Magnetisation_Array = np.array([])

    Measurement_Rate = math.ceil(2 * t_correlation * LatticeSquared)
    t_max = Amount_Blocks * 16 * t_correlation
    bar = tqdm(range(8 * Amount_Blocks))
    for i in bar:  # Create blocks of length 16*t_corr each consisting of 8 individual measurements every 2*t_corr
        for a in range(Measurement_Rate):
            x_index = random.randint(0, Lattice - 1)
            y_index = random.randint(0, Lattice - 1)

            Delta_Energy = Energy_Difference(x_index, y_index)
            if Delta_Energy <= 0:
                Spin_Array[y_index, x_index] = -Spin_Array[y_index, x_index]
                Energy += Delta_Energy
            elif Flip_Calculator(Delta_Energy):
                Spin_Array[y_index, x_index] = -Spin_Array[y_index, x_index]
                Energy += Delta_Energy
        Sum_Of_Spins = np.sum(Spin_Array)
        Magn_PSpin = np.append(Magn_PSpin, abs(Sum_Of_Spins) / Lattice ** 2)
        Energy_PSpin = np.append(Energy_PSpin, Energy / Lattice ** 2)
        Energy_Array = np.append(Energy_Array, Energy)
        Magnetisation_Array = np.append(Magnetisation_Array, Sum_Of_Spins)
        if (i + 1) % 8 == 0:
            Magnet_Suscept_Array = np.append(Magnet_Suscept_Array, Magnetic_Suscept(Magnetisation_Array))
            Spec_Heat_Array = np.append(Spec_Heat_Array, Specific_Heat(Energy_Array))
            Energy_Array = np.array([])
            Magnetisation_Array = np.array([])
    # print(Magn_PSpin)
    # print(Energy_PSpin)
    # print("Magn suscept" + str(Magnet_Suscept_Array))
    # print("Specific heat" + str(Spec_Heat_Array))
    # print(np.mean(Magn_PSpin))
    # print(np.mean(Energy_PSpin))
    # print(np.mean(Magnet_Suscept_Array))
    # print(np.mean(Spec_Heat_Array))
    f = open(folder1 + '\\' + 'output.txt', 'a')
    f.write('\nFor magnetization per spin, the mean is ' + str(np.mean(Magn_PSpin))
            + '\nthe standard deviation is ' + str(Standard_deviation(Magn_PSpin, t_max)))
    f.write('\nFor energy per spin, the mean is ' + str(np.mean(Energy_PSpin))
            + '\nthe standard deviation is ' + str(Standard_deviation(Energy_PSpin, t_max)))
    f.write('\nFor magnetic susceptibility per spin, the mean is ' + str(np.mean(Magnet_Suscept_Array))
            + '\nthe standard deviation is ' + str(Standard_deviation(Magnet_Suscept_Array, t_max)))
    f.write('\nFor specific heat per spin, the mean is ' + str(np.mean(Spec_Heat_Array))
            + '\nthe standard deviation is ' + str(Standard_deviation(Spec_Heat_Array, t_max)))
    f.close()

if Correlation_Mode:
    Calculated_Correlation_Tau = [0,0,0,0]
    if t_equilibrium < Repetitions / LatticeSquared:
        for i in range(4):
            Auto_Corr = Auto_Correlation(t_equilibrium, Multiple_Spin_Magn[i])
            Calculated_Correlation_Tau[i] = ", " + str(Correlation_Time(Auto_Corr))
        # print(Calculated_Correlation_Tau)
        f = open(folder1 + '\\' + 'output.txt', 'a')
        for j in range(4):
            Values = ''.join(Calculated_Correlation_Tau)
        f.write('\nThe correlation times of the system are: ' + Values)
        f.close()
        # plt.figure(3)
        # plt.plot(Auto_Corr)
    plot2 = plt.figure(2)
    Colours = ['r', 'b', 'g', 'k']
    for z in range(4):
        plt.plot(Multiple_Spin_Magn[z], Colours[z])
    plt.savefig(folder + '\\' + 'magnetization.png')
if Thermodynamics_Mode:
    plt.figure(2)
    plt.plot(Spin_Magn)

plt.figure(1)
imgplot = plt.imshow(Spin_Array)
plt.colorbar()
plt.savefig(folder + '\\' + 'spin_alignment.png')
# plt.show()
