# Impact of vaccinations, boosters and lockdowns on COVID-19 incidence 

This repository contains code and data for the analyses in the preprint 'Impact of vaccinations, boosters and lockdowns on COVID-19 waves in French Polynesia'. The aim of this work is to estimate the counterfactual impact of vaccinations, boosters and lockdowns on incidence of COVID-19 cases, hospitalisations and deaths in French Polynesia. We do this by fitting an age-structured multi-strain model of SARS-CoV-2 transmission to data on daily numbers of COVID-19 confirmed cases, hospitalisations and hospital deaths in French Polynesia from July 2020 to May 2022, and data from two sero-surveys conducted in February 2021 and November-December 2021 [1]. The model is written in the [odin](https://mrc-ide.github.io/odin/) syntax and compiled to work with the [dust](https://mrc-ide.github.io/dust/) discrete-time C++ simulation engine using [odin.dust](https://mrc-ide.github.io/odin.dust/) [2]. The code also uses the particle filter tools of the [mcstate](https://mrc-ide.github.io/mcstate/) R package, but implements a bespoke adaptive MCMC algorithm for fitting the model.

## Requirements
The [eigen1](https://github.com/mrc-ide/eigen1/), odin, dust, odin.dust and mcstate R packages are required to run the code. Since these packages are in constant development it is probably best to use the same versions of these packages as were used for this analysis:
```
eigen1 0.1.1
odin 1.3.2
dust 0.11.24
odin.dust 0.2.16 
mcstate 0.9.0
```
These can be installed in R with:
```r
remotes::install_github(c(
  "mrc-ide/eigen1@v0.1.1",
  "mrc-ide/odin@v1.3.2",
  "mrc-ide/dust@v0.11.24",
  "mrc-ide/odin.dust@v0.2.16",
  "mrc-ide/mcstate@v0.9.0"))
```

## Installation

Clone/download this project onto your machine.

## Data
All data required to run the code are contained in the [data](data) folder and are also available on Zenodo: [doi:10.5281/zenodo.7790839](https://doi.org/10.5281/zenodo.7790839)

## Running the code

The odin model code is contained in [inst/odin/covid_multi_strain.R](inst/odin/covid_multi_strain.R). If required, the model can be modified by editing this code, and then recompiled by calling 
```R
odin_dust("inst/odin/covid_multi_strain.R")
```
as in [fit_covid_multi_strain.R](R/fit_covid_multi_strain.R).

The workflow for fitting the model and running counterfactual simulations is as follows. The model can be fit to the data for French Polynesia using the adaptive MCMC algorithm in [pmcmc.R](R/pmcmc.R) by running:
```R
source("R/fit_FP.R")
```
This produces output called `MCMCoutput<run>.RData` in the `output` folder, where `<run>` is the MCMC `run` number set in [fit_FP.R](R/fit_FP.R), along with trace and posterior plots of the parameters produced by [plot_fit.R](R/plot_fit.R).

The breakdown of population immunity over time (Fig. 3 in the preprint) can be plotted with:
```R
source("R/plot_immunity.R")
```
setting the `run` number in [plot_immunity.R](R/plot_immunity.R) to match that in [fit_FP.R](R/fit_FP.R).

The counterfactual simulations without vaccinations, boosters, and lockdowns, and with different lockdown timings can then be run with:
```R
source("R/run_simulations.R")
```
setting `run` in [run_simulations.R](R/run_simulations.R) to match `run` in [fit_FP.R](R/fit_FP.R).

Note that the fitting code takes ~1.5hrs to run (in R 4.1.0 on an 8-core MacBook Pro with an M1 chip and 16GB RAM).

## Built With

* [R version 4.1.0 (2021-05-18)](https://www.r-project.org/)

## Authors

* Lloyd Chapman: <l.chapman4@lancaster.ac.uk> 

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE.txt](LICENSE.txt) file for details

## Citation
[Chapman LAC, Aubry M, Maset N, Russell TW, Knock ES, Lees JA, Mallet H-P, Cao-Lormeau V-M, Kucharski AJ. Impact of vaccinations, boosters and lockdowns on COVID-19 waves in French Polynesia. medRxiv. 2023.](https://doi.org/10.1101/2023.03.29.23287906)

## References
1. Aubry M, Maset N, Chapman L, Simon A, Olivier S, Bos R, Chung K, Teiti I, Kucharski A, Mallet H-P, Cao-Lormeau V-M. Seroprevalence of sars-cov-2 antibodies in french polynesia and perspective for vaccine strategies. Preprints, 2022. [doi:10.20944/preprints202212.0386.v1](https://doi.org/10.20944/preprints202212.0386.v1)

2. FitzJohn RG, Knock ES, Whittles LK, Perez-Guzman PN, Bhatia S, Guntoro F, Watson OJ, Whittaker C, Ferguson NM, Cori A, Baguelin M, and Lees JA. Reproducible parallel inference and simulation of stochastic state space models using odin, dust, and mcstate. Wellcome Open Research, 5:288, 12 2021. [doi:10.12688/wellcomeopenres.16466.2](https://doi.org/10.12688/wellcomeopenres.16466.2)

