## Table 2

The folder contains the code to generate the results for Scenario (a) presented in Table 2. The simulation aims to show the validity and efficiency of the proposed design and inference procedure for the scenario where we are interested in sampling based on both outcomes.

To run the simulation, load the scripts `ACML.R`, `DatGenFunBiv.R`and `FitImputeMI.R`reported [here](https://github.com/ChiaraDG/MultivariateODS_LMM/tree/main/Simulation%20Studies), load the `setup.R` script in this folder, and run the code in `runBiv.R` script in this folder. For all the covariates included in the model, the resulting table will contain the 1) estimated coefficients, 2) estimated variance, 2) indicator variable looking at whether the 95% confidence interval contains the true parameter value.
