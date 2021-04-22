# CMP-project2
## MOSIM
This is the public repository for the code MOSIM (MOnte carlo Simulation of Ising Model), which can be use to analyse the following quantities of 2 dimensional Ising model:
- Magnetization plot of system
- Spin alignment plot of system
- Correlation time of system

<br/>And it can also get the following mean and standard deviation values of the thermodynamics properties:
- Magnetization per spin
- Energy per spin
- Magnetic susceptibility per spin
- Specific heat per spin
### Installation
Now MOSA can be worked on **Windows** system. Before installation, make sure your system meets the prerequisites listed in **Dependencies**, list below.
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
- scipy >= 0.16.0
- tqdm >= 4.59.0
### Running
In our repository, we provide you a bat profile `test.bat` with some arguments. You can directly run this profile on your machine to test, and it will take you around 2 minutes.
```
python simulation.py --Density 0.3 --Temperature 3 --Pressure --Timestep 500 --Energy --Correlation
```
In `run.bat` profile, you can manually change the argument we provided to get quantities and plot of system you want. The following is the explaination of argument.
1. --Density
<br/>Set the initial density of system. The default value is 1.
2. --Temperature
<br/>Set the initial temperature of system. The default value is 1.
3. --Timestep
<br/>Set the timestep of system. The default value is 10000. Runing 10000 timesteps will cost you 35 minutes. Make sure this is a multiple of 100 if you want to calculate the Correlation.
4. --Energy
<br/>Can generate system energy plot in folder `./plot`, named `energy.png`.
5. --Correlation
<br/>Can generate system pair correlation function plot in folder `./plot`, named `pair_correlation.png`.
6. --Pressure
<br/>Can output the system pressure in folder `./data`, named `output.txt`.
7. --Animation2d
<br/>Can generate system 2d animation in folder `./plot`, named `animation.gif`
8. --RandomInitialisation
<br/>This is a bonus argument we made to start with a 2D system with 8 argon particles of random velocity.
### Relevant report
For more detail about MOSA, you can check out our report in this repository.
### Authors
- [Davey Plugers](https://github.com/DaveyPlugers)
- [Zhen Xiang](https://github.com/boson112358)
