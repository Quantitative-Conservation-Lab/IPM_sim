
<!-- README.md is generated from README.Rmd. Please edit that file -->

Understanding the processes that control population dynamics in wild
populations is a critical component of conservation and management but
is frequently limited by missing demographic information, uncertainty,
and temporally or spatially misaligned datasets (Schaub & Abadi 2011,
Zipkin & Saunders 2017). Integrated population models (IPMs; Besbeas et
al. 2002, Brooks et al. 2004) are an increasingly popular tool in
ecology that can help overcome these limitations by combining disparate
datasets in a single, unified analysis (Schaub & Abadi 2011). This
approach can provide several advantages, such as improved precision
(Schaub et al. 2007, Tavecchia 2009, Abadi et al. 2010) and the
estimation of parameters that would otherwise be unidentifiable (Besbeas
2005, Schaub et al. 2007, Veran & Lebreton 2008), which is often the
case for key demographic rates such as recruitment or immigration. These
advances have practical benefits for resource managers that must make
decisions based on limited available information. As IPMs become more
widely used, however, it is critical that they are robustly evaluated.
In particular, the limitations of the IPM framework are still
underexplored, particularly in situations where available data are
sparse, of poor quality, or when data to inform specific parameters are
entirely absent.

Through a simulation analysis, we will identify monitoring programs and
demographic scenarios where integrated population models may not lead to
improved precision or return unbiased estimates of parameters that would
otherwise be unidentifiable when certain datasets are excluded. We apply
this simulation study to to a passerine bird life history. We assume all
model assumptions are met and compare parameter estimation for a
traditional IPM with data informing abundance (count data), survival
(mark-recapture data), and reproductive output (nest monitoring data) to
situations when either of the latter two datasets are omitted. In doing
so, we will refine our understanding of the mechanisms underlying IPMs
and can better identify when they may serve as a useful tool or might
not be necessary or warrant the collection of additional data. Improving
our understanding of the strengths and weaknesses of the IPM framework
can provide helpful insights for the design and implementation of field
survey programs so as to most efficiently make use of limited resources
to inform management efforts aimed at monitoring or recovering species
that are deemed ecosystem indicators or depleted and in need of
conservation intervention.

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
and ‘simData’.

The ‘simPopTrajectory’ function simulates the true population trajectory
by individual and age class, starting from the stable age distribution.
Individual fates are retained throughout total number of years.

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

Main workhorse script that pulls from the scripts above to run models in
parallel and output results.

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
