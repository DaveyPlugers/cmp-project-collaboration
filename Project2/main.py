import math
import random
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

J=1
KB = 1
T=3
Number = 50
NumberSquared = Number**2
Repetitions = NumberSquared*800
#Time in units of steps per lattice
t_equilibrium = 500
t_correlation = 5
Amount_Blocks = 4
Correlation_Mode = False



Spin_Array = np.array([[int(2*(random.randint(0,1)-0.5)) for x in range(Number)]for y in range(Number)])

def Energy_Calculation():
    """
    Can be called to calculate the total energy in the system
    :return: Scalar value of the total energy of the system
    """
    Energy = 0
    for x in range(Number):
        for y in range(Number):
            Energy += -0.5*J*Spin_Array[y,x]*(Spin_Array[y,(x+1)%Number] + Spin_Array[y,(x-1)%Number] + Spin_Array[(y+1)%Number,x] + Spin_Array[(y-1)%Number,x])
    return Energy
Energy = Energy_Calculation()

def Energy_Difference(x_index,y_index):
    """
    Takes index of a spin, returns the energy difference if we would flip this spin
    :param x_index: x position of the spin
    :param y_index: y position of the spin
    :return: Scalar energy difference if the spin is flipped
    """
    Energy_Diff = 2*J*Spin_Array[y_index,x_index]*(Spin_Array[y_index,(x_index+1)%Number] + Spin_Array[y_index,(x_index-1)%Number] + Spin_Array[(y_index+1)%Number,x_index] + Spin_Array[(y_index-1)%Number,x_index])
    return Energy_Diff

Probabilities = [math.exp(-4/(KB*T)),math.exp(-8/(KB*T))]
def Flip_Calculator(Delta_Energy):
    """
    Takes energy difference, calculates whether the spin flips with probability e^(-Delta_Energy/(kb*T))
    :param Delta_Energy: Energy difference for if it is flipped, this value can only be 4 or 8 since 2d grid and values =< 0 always flip
    :return: Boolean for if it should flip
    """
    if Delta_Energy ==4:
        Boolean = random.random() < Probabilities[0]
    else:
        Boolean = random.random() < Probabilities[1]
    return Boolean

Spin_Magn = np.array([])
def Spin_Magnetization(Spin_Magn):
    """
    Takes the array of Spin Magnetisations, returns a version appended with the new value of mean magn. per spin
    :param Spin_Magn: Array of mean magn. per spin
    :return: Same array with additionally the current mean magn. per spin
    """
    return np.append(Spin_Magn,[np.sum(Spin_Array)/Number**2])


def Auto_Correlation(t_equilibrium,Spin_Magn):
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
        Auto_Corr[t] = 1 / (t_max - t) * (np.sum(Spin_Magn[t_equilibrium:t_equilibrium + t_max - t] * \
        Spin_Magn[t_equilibrium + t: t_equilibrium + t_max ])) - (1 / (t_max - t) ** 2 * \
        np.sum(Spin_Magn[t_equilibrium: t_equilibrium + t_max - t ]) * \
        np.sum(Spin_Magn[t_equilibrium + t: t_equilibrium + t_max ]))

    return Auto_Corr

def Correlation_Time(Auto_Corr):
    """
    Takes an array of autocorrelation values, returns the correlation time
    :param Auto_Corr: Array that corresponds with the discrete autocorrelation function X(t)
    :return: Scalar value of the correlation time
    """
    Tau = 0
    for t in range(len(Auto_Corr)):
        if Auto_Corr[t]>0:
            Tau += Auto_Corr[t]/Auto_Corr[0]
        else:
            return Tau
    return Tau

def Magnetic_Suscept(Magnetisation_Array):
    return 1/(KB*T*NumberSquared)*(np.mean(Magnetisation_Array**2)-np.mean(Magnetisation_Array)**2)


def Specific_Heat(Energy_Array):
    return 1 / (KB * NumberSquared * T**2) * (np.mean(Energy_Array ** 2) - np.mean(Energy_Array) ** 2)
#New structureplan: Have 2 different simulation possibilities:
# 1) Run multiple simulations on the same temperature, different initials, use these to calculate equilibrium and correlation time
# 2) Simulation where you give t_equi and t_corr, then make it run until t_equi, then start dividing in k blocks of 16t_corr
# Every 2t_corr we calculate quantities: {m, e,X, C}, we average X and C over the next 8 measurements
# And then use multiple blocks of these 16t_corr to do error analysis, for m and e we just do regular average with STD given




if Correlation_Mode:
    bar = tqdm(range(Repetitions)) #Running full simulation
else:
    bar = tqdm(range(t_equilibrium * NumberSquared)) #Running untill equilibrium then start a continuation for physical values later

for a in bar:
    x_index = random.randint(0,Number-1)
    y_index = random.randint(0,Number-1)

    Delta_Energy = Energy_Difference(x_index,y_index)
    if Delta_Energy <= 0:
        Spin_Array[y_index,x_index] = -Spin_Array[y_index,x_index]
        Energy += Delta_Energy
    elif Flip_Calculator(Delta_Energy):
        Spin_Array[y_index,x_index] = -Spin_Array[y_index,x_index]
        Energy += Delta_Energy
    if a%NumberSquared == 0:
        Spin_Magn = Spin_Magnetization(Spin_Magn)


if not Correlation_Mode:
    Magn_PSpin = []
    Energy_PSpin = []
    Magnet_Suscept_Array = []
    Spec_Heat_Array = []

    Energy_Array = np.array([]) #Can use Energy_PSpin too but this is easier
    Magnetisation_Array = np.array([])

    Measurement_Rate = math.ceil(2*t_correlation*NumberSquared)
    bar = tqdm(range(8*Amount_Blocks))
    for i in bar: #Create blocks of length 16*t_corr each consisting of 8 individual measurements every 2*t_corr
        for a in range(Measurement_Rate):
            x_index = random.randint(0, Number - 1)
            y_index = random.randint(0, Number - 1)

            Delta_Energy = Energy_Difference(x_index, y_index)
            if Delta_Energy <= 0:
                Spin_Array[y_index, x_index] = -Spin_Array[y_index, x_index]
                Energy += Delta_Energy
            elif Flip_Calculator(Delta_Energy):
                Spin_Array[y_index, x_index] = -Spin_Array[y_index, x_index]
                Energy += Delta_Energy
        Sum_Of_Spins = np.sum(Spin_Array)
        Magn_PSpin.append(abs(Sum_Of_Spins)/Number**2)
        Energy_PSpin.append(Energy/Number**2)
        Energy_Array = np.append(Energy_Array,Energy)
        Magnetisation_Array = np.append(Magnetisation_Array,Sum_Of_Spins)
        if (i+1)%8 ==0:
            print("We are in the i+1%8 thing" + str(i))
            Magnet_Suscept_Array.append(Magnetic_Suscept(Magnetisation_Array))
            Spec_Heat_Array.append(Specific_Heat(Energy_Array))
            Energy_Array = np.array([])
            Magnetisation_Array = np.array([])
    print(Magn_PSpin)
    print(Energy_PSpin)
    print("Magn suscept" + str(Magnet_Suscept_Array))
    print("Specific heat" + str(Spec_Heat_Array))



print(len(Spin_Magn))
if Correlation_Mode:
    if t_equilibrium < Repetitions/NumberSquared:
        Auto_Corr = Auto_Correlation(t_equilibrium,Spin_Magn)
        Calculated_Correlation_Tau = Correlation_Time(Auto_Corr)
        print(Calculated_Correlation_Tau)
        plt.figure(3)
        plt.plot(Auto_Corr)


plt.figure(1)
imgplot = plt.imshow(Spin_Array)
plt.colorbar()
plt.figure(2)
plt.plot(Spin_Magn)

plt.show()



