# BayesianQRE
Data and replication files for "Bayesian Inference for Quantal Response Equilibrium in Normal-Form Games"

This repository contains data and replication files for "Bayesian Inference for Quantal Response Equilibrium in Normal-Form Games" by James Bland. 

Please direct any questions to james dot bland at utoledo dot edu

These files are orgainized as follows:

* `code` Mostly contains scripts for replicating the figures and tables in the paper. You should set your working directory to this folder
		* `ExampleStan` contains code for running the expamples and simulation in the online appendix
		* `SimulationOutputs` is a blank folder that will be populated when you run `EstimateAll.R`
		* `StanFiles` contains the *Stan* programs used to estimate the models in Table 1
		* `2x2_all_data.csv`: Data from Selten & Chmura (2008) at the individual decision level (i.e. no aggregation)
		* `Appendix.Rmd`: Markdown file used to make the online appendix
		* `CompileModels.R` Compiles all of the *Stan* programs (run this first)
		* `EstimateAll.R` Estimates all of the models (run this second)
		* `FiguresAndTables.Rmd` produces some of the figures and tables in the manuscript
		* `KFOLD.R` runs the k-fold cross-validation (run this third)
		* `LoadSeltenCommon.R` is a common script to most other scripts which loads the data and some packages
		* `LuceSelectedModel.Rmd` produces the figures in Section 6.3 of the manuscript
		* `PriorsCheckSelten2008NoStan.R` produces some plots used in the prior calibration.
		* `SelectedModel.Rmd` produces some plots of the selected model in Table 1

* `data`
		* `Selten2008.xlsx` contains aggregated data from Selten & Chmura (2008)

* `outputs` is a blank folder that is populated with plots and tables when you run some of the above scripts


