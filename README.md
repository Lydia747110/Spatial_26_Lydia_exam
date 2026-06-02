# Spatial_26_Lydia_exam
Repository for spatial analytics f26.

This code should render my analysis fully reproducible.

## Structure of repository
The folder "Scripts" contains 1 R file and 2 STAN files. The STAN files contain the code for the statistical models, and are loaded by R. All analysis was run from the R file. All paths in the R file are relative to how files are structured in this repository

The folder "Raw_data" contains the raw data files necessary to generate all results.

The folder "small_posterior draws" contains 300 draws from the model posteriors used for further analysis. Note that the results in the paper are based on 2000 draws for better accuracy, but these files were much to big to be uploaded

Note that the R script generates results but doesn't save them, as they are quick to run. If you decide to fit the models yourself, I advise saving the fits as rds files using writeRDS(). as they are time consuming to fit. 

Indeholder data fra Styrelsen for Dataforsyning og Infrastruktur, Danmarks Administrative Geografiske Inddelinger , Maj 2026”
