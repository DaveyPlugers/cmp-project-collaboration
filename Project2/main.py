import math
import random
import matplotlib.pyplot as plt
import numpy as np
J=1
KB = 1
T=2.2
Number = 30
Repetitions = 200000
#  [[10 spins],repeat this for a total of 10 times] All values 1 or -1


Spin_Array = np.array([[int(2*(random.randint(0,1)-0.5)) for x in range(Number)]for y in range(Number)])


def Energy_Difference(x_index,y_index):
    Energy_Diff = 2*J*Spin_Array[y_index,x_index]*(Spin_Array[y_index,(x_index+1)%Number] + Spin_Array[y_index,(x_index-1)%Number] + Spin_Array[(y_index+1)%Number,x_index] + Spin_Array[(y_index-1)%Number,x_index])
    return Energy_Diff

def Flip_Calculator(Delta_Energy):
    probability = math.exp(-Delta_Energy/(KB*T))
    Boolean = random.random() < probability
    return Boolean

Spin_Magn = []
def Spin_Magnetization():
    Spin_Magn.append(np.sum(Spin_Array)/Number**2)


for a in range(Repetitions):
    x_index = random.randint(0,Number-1)
    y_index = random.randint(0,Number-1)

    Delta_Energy = Energy_Difference(x_index,y_index)
    if Delta_Energy <= 0:
        Spin_Array[y_index,x_index] = -Spin_Array[y_index,x_index]
    elif Flip_Calculator(Delta_Energy):
        Spin_Array[y_index,x_index] = -Spin_Array[y_index,x_index]
    if a%500 == 0:
        Spin_Magnetization()


plt.figure(1)
imgplot = plt.imshow(Spin_Array)
plt.colorbar()
plt.show()

plt.figure(2)
plt.plot(Spin_Magn)
plt.show()



