import os
import numpy as np
import matplotlib.pyplot as plt

folder = os.getcwd() + '/data_analyse'

if not os.path.exists(folder):
    os.makedirs(folder)

dataset = np.loadtxt('correlation-time.txt')
T = dataset[0]
Tau = dataset[1]

plt.scatter(T, Tau, label='data points', c='b')
plt.plot(T, Tau)
plt.xlabel('Temperature [J/$k_B$]')
plt.ylabel('correlation time ' + r'$\tau$' + '[Monte Carlo Steps per Site]')
my_x_ticks = np.arange(1, 4, 0.2)
plt.xticks(my_x_ticks)
plt.legend()
plt.savefig('./data_analyse/correlation_time.png')
plt.show()


