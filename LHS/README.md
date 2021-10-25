## Lung Health Study: Data Analysis

The files `ACML.R` and `FitImputeMI.R` contain the code to perform ACML and MI. The file `runBivLHSMI.R` contains the code to fit the six designs + procedures reported in the paper (Figure 2 and Web Table 8). 

To be able to reproduce the results one would need the original data from the Lung Health Study, load the files `ACML.R` and `FitImputeMI.R`, and run the code in `runBivLHSMI.R`. The resulting table will contain the estimated coefficients and standard errors for all the covariates included in the model.

The code to reproduce Figure 2 can be found in `plotLHS.R`.  
