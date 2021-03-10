Allpositions = [[] for x in range(2000)]

from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
BoxSize = 60


Particles = []
for i in range(6):
    for j in range(6):
        for k in range(3):
            Partic = [k*BoxSize/3 + (j+i)%2*BoxSize/6,j*BoxSize/6,i*BoxSize/6]
            Particles.append(Partic)


ax = plt.axes(projection='3d')

zdata = [x[2] for x in Particles]
xdata = [x[0] for x in Particles]
ydata = [x[1] for x in Particles]


ax.scatter3D(xdata, ydata, zdata)
plt.show()