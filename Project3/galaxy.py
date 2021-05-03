import numpy as np
import random
import matplotlib.pyplot as plt

N = 100
r = 10
x = []
y = []
z = []
vx = []
vy = []
vz = []

for i in range(N):
    theta = random.random()*2*np.pi
    phi = random.random()*np.pi
    x.append(r * np.sin(phi) * np.cos(theta))
    y.append(r * np.sin(phi) * np.sin(theta))
    z.append(r * np.cos(theta))

fig1 = plt.figure(1)
plt.scatter(x, y)
my_x_ticks = np.arange(-50, 50, 10)
my_y_ticks = np.arange(-50, 50, 10)
plt.xticks(my_x_ticks)
plt.yticks(my_y_ticks)
plt.show()
