# genetic_trade_offs_simulations
This repository contains the code necessary to replicate the simulations in Connallon et al. : Predicting the prevalence of genetic trade-offs among adaptive substitutions.

The goal of this repository is to enhance reproducibility and to enable future research. For any questions, please post [here](https://github.com/ldutoit/genetic_trade_offs_simulations/issues).

The following files are included:
  - [core_simulations.R](core_simulations.R) . This file runs the simulations and produces the visualisation included in the manuscript.
  - [fixSA.netbenefit.R](fixSA.netbenefit.R) .This file defines a function that runs the Wright-Fisher distribution of fitness effects with sexes
  - [num_pred&visualise.R](num_pred&visualise.R) . This file contains the numerical prediction function num_pred() and the visualise() function producing the simulations figure in the manuscript.  

A set of example outputs is included in the folder [example_results](example_results).
