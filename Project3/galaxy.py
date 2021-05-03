import numpy as np
import random
import matplotlib.pyplot as plt

N = 1000
R = 10
x = []
y = []
z = []
vx = []
vy = []
vz = []

for i in range(N):
    theta = random.random()*2*np.pi
    phi = random.random()*np.pi
    r = random.random()*R
    x.append(r * np.sin(phi) * np.cos(theta))
    y.append(r * np.sin(phi) * np.sin(theta))
    z.append(r * np.cos(phi))

fig1 = plt.figure()
plt.scatter(x, y)
my_x_ticks = np.arange(-50, 50, 10)
my_y_ticks = np.arange(-50, 50, 10)
plt.xticks(my_x_ticks)
plt.yticks(my_y_ticks)

fig2 = plt.figure()
ax = fig2.add_subplot(1, 1, 1, projection='3d')
ax.scatter(x, y, z)
plt.show()
