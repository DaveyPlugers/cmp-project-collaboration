import numpy as np
import random
import matplotlib.pyplot as plt

N = 100
M = 0.333
R = 10
a = 0.1
G = 1
x = []
y = []
z = []
vx = []
vy = []
vz = []
rs = []
V = []

for i in range(N):
    while True:
        m = np.random.uniform(0, M)
        r = a / ((M / m) ** 0.5 - 1)
        # print(r)
        if r < R:
            break
    rs.append(r)
    theta = random.random()*2*np.pi
    phi = random.random()*np.pi
    x.append(r * np.sin(phi) * np.cos(theta))
    y.append(r * np.sin(phi) * np.sin(theta))
    z.append(r * np.cos(phi))
    Vcirc = (G*m/r) ** 0.5
    V.append(Vcirc)
    vx.append(Vcirc * np.sin(phi) * np.cos(theta))
    vy.append(Vcirc * np.sin(phi) * np.sin(theta))
    vz.append(Vcirc * np.cos(phi))

fig1 = plt.figure()
plt.scatter(x, y, marker='.')
my_x_ticks = np.arange(-5, 5, 1)
my_y_ticks = np.arange(-5, 5, 1)
plt.xticks(my_x_ticks)
plt.yticks(my_y_ticks)

fig2 = plt.figure()
ax = fig2.add_subplot(1, 1, 1, projection='3d')
ax.scatter(x, y, z)

fig3 = plt.figure()
plt.scatter(rs, V)
plt.show()
