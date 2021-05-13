### Simulation Studies

The folder contains .R files to replicate the results from the simulation studies reported in the paper:

* [ACML.R](https://github.com/ChiaraDG/MultivariateODS_LMM/blob/main/Simulation%20Studies/ACML.R). Contains the code to run ODS and BDS using ascertainment corrected likelihood

* [DatGenFunBiv.R](https://github.com/ChiaraDG/MultivariateODS_LMM/blob/main/Simulation%20Studies/DatGenFunBiv.R). Contains the code to generate bivariate LMM data, and identify "informative subjects" for the two-phase design under ODS and BDS

* [FitImputeFunsIIA.R](https://github.com/ChiaraDG/MultivariateODS_LMM/blob/main/Simulation%20Studies/FitImputeFunsIIA.R). Contains the code necessary to run IIA

* [Table 2](https://github.com/ChiaraDG/MultivariateODS_LMM/tree/main/Simulation%20Studies/Table%202). Contains the simulation setup used for Table 2 (`setup.R`) as well as the R file (`run.R`) needed to run all the code under different designs and inference procedures (ACML or IIA)

* [Table 3](https://github.com/ChiaraDG/MultivariateODS_LMM/tree/main/Simulation%20Studies/Table%203). Contains the simulation setup used for Table 3 (`setup.R`) as well as the R file (`run.R`) needed to run all the code under different designs and inference procedures (ACML or IIA)

* [Web Tables](https://github.com/ChiaraDG/MultivariateODS_LMM/tree/main/Simulation%20Studies/Web%20Tables). Contains the code to replicate the results from Web Tables 1-4. Each sub-folder contains the files `setup.R` and `run.R` used for the simulation.