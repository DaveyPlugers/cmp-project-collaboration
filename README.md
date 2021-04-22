# CMP-project2
## MOSIM
This is the public repository for the code MOSIM (MOnte carlo Simulation of Ising Model), which can be use to analyse the following quantities of 2 dimensional Ising model:
- Mean magnetization plot of system
- Spin alignment plot of system
- Correlation time of system

<br/>And it can also get the following mean and standard deviation values of the thermodynamics properties:
- Magnetization per spin
- Energy per spin
- Magnetic susceptibility per spin
- Specific heat per spin
### Installation
Now MOSIM can be worked on **Windows** system. Before installation, make sure your system meets the prerequisites listed in **Dependencies**, list below.
<br/>To install, there are two ways:
1. **Download this repository**
<br/>You can go to our [github page](https://github.com/boson112358/cmp-project) to download it and save it on your machine
2. **Clone this repository**
<br/>Clone it on your machine:
```
git clone https://github.com/boson112358/cmp-project
```
#### Dependencies
The code requires the following(some of the functionality might work with older versions):
- Python version 3.9+
- argparse >= 1.4.0
- matplotlib >= 3.3.4
- numpy >= 1.20.1
- tqdm >= 4.59.0
### Running
For our code, we contain 2 mode to analyse different quantities, **correlation mode** and **thermodynamic mode**.

For the **correlation mode**, it can be used to generate a final spin alignment plot, mean magnetization plot(can be used to estimate equilibrium time), and get 4 estimates of correlation time. For the value of correlation time, it is output in the folder `./data`, named `output.txt`. And find the plot in the folder `./plot`.
In our repository, we provide you a bat profile `test_correlation.bat` with some arguments. You can use it to test the correlation mode. And you have to wait to get 4 simulations of the system. And for this mode, you can use following arguments: Temperature, Lattice, Total_Duration, Correlation_Mode, No_plot.
```
python main.py --Correlation_Mode --Temperature 1 --Total_Duration 2000
```

For the **thermodynamic mode**, It can output mean and standard deviation values of the thermodynamics properties in folder `./data`, named `output.txt`. And with a final spin alignment plot, in the folder `./plot`. In our repository, we provide you a bat profile `test_thermo.bat` with some arguments. You can use it to test the thermodynamic mode. And you first need to wait the system to reach equilibrium and then do the thermo anlysis. And for this mode, you need find the system's correlation time and equilibrium time first in order to calculate thermodynamic properties. We provide you a list for correlation time for some temperature, in `correlation.txt`. For this mode, you can use following arguments: Temperature, Lattice, Equilibrium_time, Correlation_time, Total_Duration, Data_Points, Thermodynamic_Mode, No_Plot.
```
python main.py --Temperature 1 --Data_Points 40 --Correlation_time 1.54 --Thermodynamic_Mode
```

In `run.bat` profile, you can manually change the argument we provided to get quantities and plot of system you want. The following is the explaination of argument.
1. --Temperature
<br/>Set the initial temperature of system. The default value is 1.
2. --Lattice
<br/>Set the width of the grid in spin sites. The default value is 50.
3. --Equilibrium_time
<br/>Set the time duration until system considered in equilibrium. The default value is 1000.
4. --Correlation_time
<br/>Set the system correlation time. The default value is 1.
5. --Total_Duration
<br/>Set the time duration the total simulation will run in correlation mode. The default value is 3000, but in system the value is `Lattice ** 2 * 3000` .
6. --Data_Points
<br/>The amount of independent blocks to calculate thermal properties. The default value is 20.
7. --Correlation_Mode
<br/>Turn on the correlation mode.
8. --Thermodynamic_Mode
<br/>Turn on the thermodynamic mode. 
9. --No_Plot
<br/>Allows plots to be turned off. Otherwise, when using these two modes, you can find the corresponding plots in folder `./plot`.
### Relevant report
For more detail about MOSIM, you can check out our report in this repository.
### Authors
- [Davey Plugers](https://github.com/DaveyPlugers)
- [Zhen Xiang](https://github.com/boson112358)
