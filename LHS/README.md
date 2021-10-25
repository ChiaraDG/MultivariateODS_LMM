## Lung Health Study: Data Analysis

The files `ACML.R` and `FitImputeMILHS.R` contain the code to perform ACML and MI. The file `runBivLHS.R` contains the code to fit the six designs + procedures reported in the paper (Figure 2 and Web Table 8). 

To be able to reproduce the results one would need the original data from the Lung Health Study, load the files `ACML.R` and `FitImputeFunsMILHS.R`, and run the code in `runBivLHS.R`.

After running `runBivLHS.R`, the resulting table will contain:

1. Estimated coefficients and standard errors for all the covariates includes in the model.
