import math
import random
import matplotlib.pyplot as plt
import numpy as np
J=1
KB = 1
T=1.4
Number = 50
NumberSquared = Number**2
Repetitions = 2000000
#  [[10 spins],repeat this for a total of 10 times] All values 1 or -1


Spin_Array = np.array([[int(2*(random.randint(0,1)-0.5)) for x in range(Number)]for y in range(Number)])

def Energy_Calculation():
    Energy = 0
    for x in range(Number):
        for y in range(Number):
            Energy += -0.5*J*Spin_Array[y,x]*(Spin_Array[y,(x+1)%Number] + Spin_Array[y,(x-1)%Number] + Spin_Array[(y+1)%Number,x] + Spin_Array[(y-1)%Number,x])
    return Energy
Energy = Energy_Calculation()

def Energy_Difference(x_index,y_index):
    Energy_Diff = 2*J*Spin_Array[y_index,x_index]*(Spin_Array[y_index,(x_index+1)%Number] + Spin_Array[y_index,(x_index-1)%Number] + Spin_Array[(y_index+1)%Number,x_index] + Spin_Array[(y_index-1)%Number,x_index])
    return Energy_Diff

Probabilities = [math.exp(-4/(KB*T)),math.exp(-8/(KB*T))]
def Flip_Calculator(Delta_Energy):
    if Delta_Energy ==4:
        Boolean = random.random() < Probabilities[0]
    else:
        Boolean = random.random() < Probabilities[1]
    return Boolean

Spin_Magn = np.array([])
def Spin_Magnetization(Spin_Magn):
    return np.append(Spin_Magn,[np.sum(Spin_Array)/Number**2])


def Auto_Correlation(t_equilibrium,Spin_Magn):
    t_max = len(Spin_Magn) - t_equilibrium
    Auto_Corr = np.zeros(t_max)
    for t in range(t_max):
        Auto_Corr[t] = 1 / (t_max - t) * (np.sum(Spin_Magn[t_equilibrium:t_equilibrium + t_max - t] * \
        Spin_Magn[t_equilibrium + t: t_equilibrium + t_max ])) - (1 / (t_max - t) ** 2 * \
        np.sum(Spin_Magn[t_equilibrium: t_equilibrium + t_max - t ]) * \
        np.sum(Spin_Magn[t_equilibrium + t: t_equilibrium + t_max ]))

    return Auto_Corr

def Correlation_Time(Auto_Corr):
    Tau = 0
    for t in range(len(Auto_Corr)):
        if Auto_Corr[t]>0:
            Tau += Auto_Corr[t]/Auto_Corr[0]
        else:
            return Tau
    return Tau



for a in range(Repetitions):
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



print(len(Spin_Magn))
Auto_Corr = Auto_Correlation(500,Spin_Magn)
Correlation_Tau = Correlation_Time(Auto_Corr)
print(Correlation_Tau)


plt.figure(1)
imgplot = plt.imshow(Spin_Array)
plt.colorbar()
plt.figure(2)
plt.plot(Spin_Magn)
plt.figure(3)
plt.plot(Auto_Corr)
plt.show()



