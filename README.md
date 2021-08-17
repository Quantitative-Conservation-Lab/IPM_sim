
<!-- README.md is generated from README.Rmd. Please edit that file -->

Something here about background?

# Scripts, current version

## 0 - preparing scenarios

### compute\_time\_calc.R

Code for approximating the total run time, given number of simulated
datasets, number of parameter scenarios, computer cores, and estimated
run time for each model.

### generate\_scenarios.R

Function ‘getNviable’ for finding the combinations of parameters that
give a population growth rate (*λ*) within certain bounds, using the
eigenvalue from the Leslie matrix. The output from this function
produces the parameter scenarios used to simulate population
trajectories, observation data from the trajectories, and fit the the
IPM models.

## 1 - simulating data

### IPM\_sim\_2.0function.R

Script contains two functions for simulating data: ‘simPopTrajectory’
and ‘simData’. The ‘simPopTrajectory’ function simulates the true
population trajectory by individual and age class, starting from the
stable age distribution. Individual fates are retained throughout total
number of years.

The ‘simData’ function uses the output from the ‘simPopTrajectory’
function to simulate the observation process for each data type. The
true individual data is broken down into independent sets, each is
subject to observation error. Mark-resight, count, and reproductive
success data are output from this function.

### simulateTrajectories.R

Script that pulls from the scenarios that were created in ‘0 - preparing
scenarios’ and the functions in IPM\_sim\_2.0function.R to simulate
observation data sets for model fitting.

## 2 - models

### IPMinitvalues.R

Functions for putting data in the form required by the survival model.

### IPM\_marray.R

Script contains 4 IPM NIMBLE models.

1.  full IPM model (‘IPMmod’) where all three datasets are included
2.  model without productivity dataset (‘nonests’)
3.  model without survival model (‘nomr’)
4.  abundance only model without productivity or survival (‘abundonly’)

## 3 - run models

### run\_scenarios\_helperFns.R

Script with functions and code to run NIMBLE models. Function ‘marray’
transforms the survival data capture histories to m-array format. MCMC
settings in this script, with functions for each of the 4 models. Each
model has a run model function that defines constants, initial values,
data, parameters to monitor, and code to build/run the models in NIMBLE.

### run\_scenarios.R

Script that pulls from the scripts above to run models in parallel and
output results.

## 4 - process results

### 01\_load\_and\_process\_data.R

Script that checks the Gelman-Rubin convergence statistic for each
output MCMC sample, then thins the chains to reduce file size and places
the output into categories based on the population growth rate, *λ*.

### 02\_compute\_geom\_means.R

Computes the geometric means for *λ* in each MCMC output.

# Data

Contains true parameter scenarios used to simulate population
trajectories and observation data, all viable parameter combinations in
files ‘-.lam.combos.RDS’ and the selected parameter combinations used to
simulate data are in files ‘-.lam.params.RDS’ *(where “-” is low, med,
or high)*. Folders labeled ‘-Trajectories’ contain the simulated
population trajectories for each of the scenarios found in
‘-.lam.params.RDS’. Each simulated trajectory is labeled with the
scenario and the simulation number for that scenario. For example, in
the ‘highTrajectory’ folder: ‘highpopTraj-1-2.RDS’ corresponds to the
first parameter scenario in ‘high.lam.param.RDS’ and is the second
population trajectory simulated using those parameter values.

File ‘true.vals.csv’ is the true parameter values used in each scenario.
‘scenario\_ID.csv’ tracks the combinations of detection probabilities
and IPM models.
