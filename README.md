# mtbCfuCeqDynamics
These are codes used in our study of dynamics of CFUs and CEQs of Mycobacterium tuberculosis in Mice and Monkeys.

## stochasticSimulations ##
For improved flexibility and speed, we implemented the adaptive tau leaping algorithm with efficient time step selection (Cau, Gillespie, and Petzold, J Chem Phys 2006) in Fortran 08. We developed classes in Python 3 for designing kinetic models, setting up simulation parmeters, and automatically generating input files, running the Fortran code, then parsing the Fortran output into a Python dictionary.

### Files ###
* stochastic.f08 - This is the Fortran code that actually executes the stochastic simulations.\
* StochasticSimulation.py - This Python 3 class builds input files, runs the fortran program, and returns the trajectories in a dictionary.\
* KineticModel.py - This Python 3 class helps build the kinetic model of the reaction network to be simulated.\
* ddModelStochastic.py - This Python 3 script sets up and runs simulations of the dependent dynamics model of CFU and CEQ dynamics described in our paper, and generates some basic plots of the output.

### Notes ###

We used Pythone 3.11 for our simulations; we are not aware of anything in the code that uniquely requires this version. On the other hand, the Fortran code requires Fortran 08 or later, as it uses the log_gamma function, which became standard in Fortran 08. We compiled our code with gfortran. The code has been tested on Ubuntu 22.04 and macOS 12.5 Monterey.

Instructions to run the simulations with these codes:
1. Compile stochastic.f08 with a Fortran compiler. We used gfortran. Call the execuatable stochastic.exe.\
      code: gfortran stochastic.f08 -o stochastic.exe
2. Decide on a default location for the input file. As you will see in our code, we created a special directory for the file in ~/local/stochastic.
3. Correct the path to the input file and the executable on lines 6 and 10 respectively, in stochasticSimulation.py (lines setting the variables infile and cmd). 
4. Be sure your Python 3 installation has the necessary packages/modules:
  i. subprocess (native to Python)
  ii. numpy (3rd party but included with many bundled distributions, including Anaconda)
  iii. matplotlib (3rd party but included with many bundled distributions, including Anaconda)
5. With all the python files in the same directory, run ddModelStochastic.py.

## powerAnalysis ##

### File ###
* powerAnalysis.ipynb - This is a simplified version of the code we used for our power analysis.

### Notes ###

We performed our power analysis in Python 3.11, with a Jupyter notebook, which can be run through the web browser using Jupyter, or in an IDE with a Jupyter extension (such as VS code).

This Python 3 code requires the following modules/packages:
* numpy
* scipy
* math (native)

We ran our analysis using Python 3.11; it is likely that it will function on other recent versions of Python 3, though we have not tested this.

## fitting ##

### Files ###

* fittingModelsToCFUCEQ.ipynb - Jupyter notebook with Python 3 code to fit DD and ID models to mean or median CFU and CEQ data
* fig2a.csv  - digitized data from Fig. 2A from Munoz Elias et al Infection and Immunity 2005.
* munozEliasFig2aMaxZ.csv - digitized data plus assumed values of Z at 1 and 14 days post-infection.

### Notes ###

Although the models can in principle be fit to the data by minimizing residuals (for example, using SciPy.minimize), we instead produced exact fits to the data over each time interval between pairs of experimental time points, since we are interested in changes in rates over time.

This Python 3 code requires the following modules/packages:
* numpy
* matplotlib
* pandas
* math (native)
