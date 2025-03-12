## README.md

This folder includes the code necessary to produce supplementary figure 5.

The following files are included: 

- [Supp_run_simulations.R](Supp_run_simulations.R). This R script runs the symmetric and assymetric simulations as defined in functions in the file *SuppFig-functions.R*.
- [SuppFig-AsymmetricCase_R1.R](SuppFig-AsymmetricCase_R1.R). This R script produces the figure using the results in the folder simResults.
- [SuppFig_Assym_Arc_n50.pdf](SuppFig_Assym_Arc_n50.pdf). This is the figure.
- [suppFig-functions.R](suppFig-functions.R) This Rscript contains all the function needed by *SuppFig-AsymmetricCase_R1.R* to generate the figure. It also includes the functions that generate the symmetric and asymmetric results in *simResults/*. 
- [simResults](simResults). Symmetric and assymetric simulations produced by the code commented out lines 606-616 of *Supp_run_simulations.R*
