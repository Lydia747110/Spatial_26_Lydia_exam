# Spatial_26_Lydia_exam
Repository for spatial analytics f26

## Structure of repository
The folder "Scripts" contains 1 R file and 2 STAN files. The STAN files contain the code for the statistical models, and are loaded by R. All analysis was run from the R file.

The folder "Raw data" contains the raw data files necessary to generate all results. (Also BBR?)

The folder "Posterior draws" contains the draws from the model posteriors used for further analysis. It also contains the full model fits

Note that the R script generates results but doesn't save them, as they are quick to run. If you decide to fit the models yourself, I advise saving the fits as rds files using writeRDS(). as they are time consuming to fit. 
