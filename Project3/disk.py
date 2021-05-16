import numpy as np
import random
import matplotlib.pyplot as plt

Number_disk = 100
Mass_disk = 1
Max_R = 15
Disk_potential_parameter = 1
G = 1
x = []
y = []
z = []
vx = []
vy = []
vz = []
rs = []
V = []

for i in range(Number_disk):
    while True:
        m = np.random.uniform(0, Mass_disk)
        r = ((m/Mass_disk)**(1/3.))*Disk_potential_parameter*((1-(m/Mass_disk)**(2/3.))**-0.5)
        print(r)
        if r < Max_R:
            break
    rs.append(r)
    theta = np.random.uniform(0, 2*np.pi)
    # phi = random.random()*np.pi
    x.append(r * np.cos(theta))
    y.append(r * np.sin(theta))
    z.append(0)
    Vcirc = (G*m/r) ** 0.5
    V.append(Vcirc)
    vx.append(Vcirc * -np.sin(theta))
    vy.append(Vcirc * np.cos(theta))
    vz.append(0)


fig1 = plt.figure()
plt.scatter(x, y, marker='.')
my_x_ticks = np.arange(-10, 10, 1)
my_y_ticks = np.arange(-10, 10, 1)
plt.xticks(my_x_ticks)
plt.yticks(my_y_ticks)

fig2 = plt.figure()
ax = fig2.add_subplot(1, 1, 1, projection='3d')
ax.scatter(x, y, z)

fig3 = plt.figure()
plt.scatter(rs, V)
plt.show()