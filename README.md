# Project overview
This project evaluates the sources of arsenic exposure for villagers in rural Bangladesh who consume arsenic-contaminated well water.  A mass balance approach paired with linear regressions is used to estimate the average fraction of a Bangladesh villager’s drinking water that comes from the individual's primary well, other wells within the individual's household compound, and other wells throughout the study area.

# Required packages
Install from Pipfile.

# Code structure
* The functions in 'run_all.py' set up the input parameters, load the data, and run the analysis.
* The functions in 'regressions.py' run linear regressions on the input data for two different mass-balance models of water and arsenic consumption and excretion.
* The functions in 'solve_mass_balance.py' use the input parameters and the parameters from the linear regressions to solve for the estimated fractions of water consumed from different sources and the uncertainties on these fractions, saving the results to csv files.
* The functions in 'plots.py' make scatter plots of the original data compared with the model predictions.

# Running the code
The code can be run with 'python3 run_all.py'. However, because the data used in this project contain sensitive information, they have not been included in this repo.

