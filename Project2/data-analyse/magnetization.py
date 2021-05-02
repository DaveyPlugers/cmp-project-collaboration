import os
import numpy as np
import matplotlib.pyplot as plt

folder = os.getcwd() + '/data_analyse'

if not os.path.exists(folder):
    os.makedirs(folder)

dataset1 = np.loadtxt('magnetization_per_spin.txt')
T = dataset1[0]
T1 = np.linspace(1, 2.25, 1000)
M = (1 - (1 / 2 * (np.exp(2 / T1) - np.exp(-2 / T1))) ** (-4)) ** (1 / 8)
m = dataset1[1]
m_error = np.array([dataset1[2], dataset1[2]])

dataset2 = np.loadtxt('energy_per_spin.txt')
e = dataset2[0]
e_error = np.array([dataset2[1], dataset2[1]])

dataset3 = np.loadtxt('magnetic_susceptibility_per_spin.txt')
X = dataset3[0]
dataset3[1] = (8**0.5)*dataset3[1]
X_error = np.array([dataset3[1], dataset3[1]])

dataset4 = np.loadtxt('specific_heat_per_spin.txt')
C = dataset4[0]
dataset4[1] = (8**0.5)*dataset4[1]
C_error = np.array([dataset4[1], dataset4[1]])


plot1 = plt.figure(1)
plt.scatter(T, m, label='data points', c='b')
plt.plot(T, m, linestyle='--')
plt.plot(T1, M, label='theoretical curve')
plt.xlabel('Temperature [J/$k_B$]')
plt.ylabel(r'$\langle|m|\rangle$'+'per site')
my_x_ticks = np.arange(1, 4, 0.2)
plt.xticks(my_x_ticks)
plt.errorbar(T, m, yerr=m_error, label='error bar', ecolor='r', elinewidth=15, linestyle='none')
plt.legend()
plt.savefig('./data_analyse/Magnetization_per_spin.png')

plot2 = plt.figure(2)
plt.scatter(T, e, label='data points', c='b')
plt.plot(T, e, linestyle='--')
plt.xlabel('Temperature [J/$k_B$]')
plt.ylabel(r'$\langle e\rangle$'+'per site[J]')
my_x_ticks = np.arange(1, 4, 0.2)
plt.xticks(my_x_ticks)
plt.errorbar(T, e, yerr=e_error, label='error bar', ecolor='r', elinewidth=15, linestyle='none')
plt.legend()
plt.savefig('./data_analyse/Energy_per_spin.png')

plot3 = plt.figure(3)
plt.scatter(T, X, label='data points', c='b')
plt.plot(T, X, linestyle='--')
plt.xlabel('Temperature [J/$k_B$]')
plt.ylabel(r'$\langle\chi_M\rangle$'+'['+r'$\mu/k_B$]')
my_x_ticks = np.arange(1, 4, 0.2)
plt.xticks(my_x_ticks)
plt.errorbar(T, X, yerr=X_error, label='error bar', ecolor='r', elinewidth=3, linestyle='none')
plt.legend()
plt.savefig('./data_analyse/magnetic_susceptibility_per_spin.png')

plot4 = plt.figure(4)
plt.scatter(T, C, label='data points', c='b')
plt.plot(T, C, linestyle='--')
plt.xlabel('Temperature [J/$k_B$]')
plt.ylabel(r'$\langle C\rangle$'+'['+r'$J/k_B^2$]')
my_x_ticks = np.arange(1, 4, 0.2)
plt.xticks(my_x_ticks)
plt.errorbar(T, C, yerr=C_error, label='error bar', ecolor='r', elinewidth=2, linestyle='none')
plt.legend()
plt.savefig('./data_analyse/specific_heat_per_spin.png')
plt.show()
