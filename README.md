# Predicting fine-scale distributions and emergent spatiotemporal patterns from temporally dynamic step selection simulations

The code, data and R objects in this repository accompany a paper titled 'Predicting fine-scale distributions and emergent spatiotemporal patterns from temporally dynamic step selection simulations', which is currently available as a preprint at: [https://www.biorxiv.org/content/10.1101/2024.03.19.585696v4](https://www.biorxiv.org/content/10.1101/2024.03.19.585696v4).

In this paper we used harmonic terms to estimate temporally dynamic coefficients from step selection models, from which we simulated animal movement trajectories. The simulations with temporal dynamics gave informative hourly predictions of expected buffalo distribution (animations below), and also gave more accurate long-term predictions. 

The R Markdown and accompanying knitted html files are numbered by the order of analysis, with an additional 'walkthrough' to build some intuition around fitting models with harmonic terms. In the bioRxiv_v3 folder there are scripts for estimating previous space use density using kernel density estimation with temporally-decaying weights, which are accompanied by a preprint of that version: [https://www.biorxiv.org/content/10.1101/2024.03.19.585696v3?versioned=true](https://www.biorxiv.org/content/10.1101/2024.03.19.585696v3?versioned=true). The memory process was not included in the models in the current version of the preprint, as the focus was to generate landscape-scale predictions, and the memory process is not required. For generating realistic geographic space use (i.e. home ranges) then the memory process helps to constrain the extent of the simulated animal's movement. To see the results of the simulations with the memory process view Version 3 of the bioRxiv preprint.

The larger files such as spatial covariates can be found in a Zenodo folder at: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10838069.svg)](https://doi.org/10.5281/zenodo.10838069). In this folder, the `data` folder contains the data used in the first script to generate random steps and sample from the covariates (which are included in the `mapping` folder). All the code should run with just these two sets of inputs (buffalo.csv for GPS data and rasters of spatial covariates in the mapping folder).

The details of the inputs and outputs of each script are below. Creating a project or setting a working directory in the `DynSSF` folder should ensure that all data and files are read in without having to change path names in the scripts.

## Animations of hourly distributions ##

The observed buffalo locations for a given hour are shown as the white locations. There are locations from several individual buffalo in this landscape extent.

![](https://github.com/swforrest/dynamic_SSF_sims/blob/main/sim_preds_0p_hourly.gif)
![](https://github.com/swforrest/dynamic_SSF_sims/blob/main/sim_preds_1p_hourly.gif)
![](https://github.com/swforrest/dynamic_SSF_sims/blob/main/sim_preds_2p_hourly.gif)
![](https://github.com/swforrest/dynamic_SSF_sims/blob/main/sim_preds_3p_hourly.gif)


## Walkthrough script - DynSSF_Walkthrough_Harmonics_and_selection_surfaces.Rmd ##

This script is a walkthrough to build intuition around fitting models with harmonic interaction terms (harmonic regression)

**Inputs**

* Data file that contains random steps, sampled covariate values and KDE density at the end of each used and random step (when using memory parameters that were estimated for all individuals simultaneously)
  * buffalo_popn_GvM_covs_ST_KDEmem1000_allOPTIM_10rs_2024-02-05.csv

**Outputs**

* Fitted model objects (2 pairs of harmonics as an example)
  * model_twostep_2p_harms_dry.rds


## Script 1 - DynamicSSF_1_step_generation.qmd ##

Generate a track object to sample random step lengths and turning angles, and sample covariate values at the end of each step.

**Inputs**

* Buffalo GPS data
  * buffalo.csv
* spatial covariate rasters
  * canopy_cover.tif
  * ndvi_GEE_projected_watermask20230207.tif
  * slope_raster.tif
  * veg_herby.tif

**Outputs**

* Data file that contains random steps and sampled covariate values at the end of each used and random step
  * buffalo_parametric_popn_covs_GvM_10rs_2024-09-04.csv



## Script 2a/b - DynamicSSF_2a_Model_fit_dry_season.qmd/DynamicSSF_2b_Model_fit_wet_season.qmd ##

This script fits step selection models with varying numbers of harmonic terms to the buffalo's dry (a) or wet (b) season data. 

**Inputs**

* Data file that contains random steps, sampled covariate values and KDE density at the end of each used and random step (when using memory parameters that were estimated for all individuals simultaneously)
  * buffalo_parametric_popn_covs_GvM_10rs_2024-09-04.csv

**Outputs**

* Fitted model objects (0, 1, 2 or 3 pairs of harmonics)
* Tables of hourly coefficient values from fitted models 



## Script 3 - DynamicSSF_3_Simulating_trajectories.qmd ##

This script takes the hourly coefficients of the fitted models and simulates (dynamic) animal movement trajectories.

**Inputs**

* Tables of hourly coefficient values from fitted models (change the input file to change which model's coefficients are used to simulate trajectories) 
* Spatial covariate rasters
  * ndvi_GEE_projected_watermask20230207.tif
  * canopy_cover.tif
  * veg_herby.tif
  * slope_raster.tif

**Outputs**

* Simulated data (animal movement trajectories)
* 

## Script 4a - DynamicSSF_4a_Aggregating_simulations.R ##

This script takes animal movement trajectories (observed and simulated), and calculates the landscape-scale prediction maps by aggregating the simulated locations into raster. We ran this script on the QUT HPC.


## Script 4b - DynamicSSF_4b_Simulation_convergence.R ##

This script assesses the convergence of the simulated trajectories, to ensure that the landscape-scale prediction maps are stable. We ran this script on the QUT HPC.


## Script 5 - DynamicSSF_5_Hourly_summaries.qmd ##

We calculate hourly summary statistics which we used to compare between the observed data and simulated data.

**Inputs**

* Simulated data (animal movement trajectories)
* Spatial covariate rasters
  * ndvi_GEE_projected_watermask20230207.tif
  * canopy_cover.tif
  * veg_herby.tif
  * slope_raster.tif
* Data file of buffalo GPS data that contains random steps and sampled covariate values

**Outputs**

* Hourly summary statistic values for observed and simulated data (with number of harmonics depending on which trajectories are read in)


## Script 6 - DynamicSSF_6_Assessing_predictions.qmd ##

This script assesses the long-term and hourly prediction maps, and validates them against the observed buffalo data using the continuous Boyce index, which is appropriate for presence-only GPS data.

**Inputs**

* Prediction maps from the `DynamicSSF_4a_Aggregating_simulations.R` script.

**Outputs**

* Plots of the landscape-scale predictions for the late dry season (~ 4 months) and for each hour of the day.
* Validation results from the continuous Boyce index

For more details, don't hesitate to get in contact at scott.forrest@qut.edu.au.
