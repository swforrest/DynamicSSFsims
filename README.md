# Simulating animal movement trajectories from temporally dynamic step selection functions

The code, data and R objects in this repository accompany a paper titled 'Predicting fine-scale distributions and emergent spatiotemporal patterns from temporally dynamic step selection simulations', which is currently available as a preprint at: [https://www.biorxiv.org/content/10.1101/2024.03.19.585696v4](https://www.biorxiv.org/content/10.1101/2024.03.19.585696v4).

In this paper we used harmonic terms to estimate temporally dynamic coefficients from step selection models, from which we simulated animal movement trajectories. The simulations with temporal dynamics gave informative hourly predictions of expected buffalo distribution (animations below), and also gave more accurate long-term predictions. 

The R Markdown and accompanying knitted html files are numbered by the order of analysis, with an additional 'walkthrough' to build some intuition around fitting models with harmonic terms. In the bioRxiv_v3 folder there are scripts for estimating previous space use density using kernel density estimation with temporally-decaying weights, which are accompanied by a preprint of that version: [https://www.biorxiv.org/content/10.1101/2024.03.19.585696v3?versioned=true](https://www.biorxiv.org/content/10.1101/2024.03.19.585696v3?versioned=true). The memory process was not included in the models in the current version of the preprint, as the focus was to generate landscape-scale predictions, and the memory process is not required. For generating realistic geographic space use (i.e. home ranges) then the memory process helps to constrain the extent of the simulated animal's movement. To see the results of the simulations with the memory process view Version 3 of the bioRxiv preprint.

The `data` folder contains the data used in the first script to generate random steps and sample from the covariates (which are included in the `mapping` folder). All the code should run with just these two sets of inputs (buffalo.csv for GPS data and rasters of spatial covariates in the mapping folder), although I have also included the intermediate outputs in the `outputs` folder that are produced by scripts and are used as inputs in the succeeding scripts. 

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


## Script 1 - DynSSF_1_step_generation.Rmd ##

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
  * buffalo_parametric_popn_covs_GvM_10rs_2024-02-18.csv



## Script 2a - DynSSF_2a_memory_parameter_estimation_hour_subset.Rmd ##

This script estimates the most likely parameters for spatial (using Kernel Density Estimation - KDE) memory with a temporal decay (negative exponential) component. It subsets the previous locations by time (rather than a fixed number of locations as in the following script), which may be more appropriate for irregularly sampled data.

**Inputs**

* Data file that contains random steps and sampled covariate values at the end of each used and random step
  * buffalo_parametric_popn_covs_GvM_10rs_2024-02-18.csv

**Outputs**

* Estimated temporal decay parameters when optimising over each individual separately (as a list of optim objects)
  * temporal_decay_param_hours_list.rds
* Estimated KDE bandwidth and temporal decay parameters when optimising over each individual separately (as a csv)
  * memory_params_KDE_exp_decay_hour_subset_2024-02-20.csv
* Estimated temporal decay parameter when using the KDE bandwiddth mean and optimising over all individuals simultanously (as an optim object)
  * optim_space_time_Gamma_param_ALLoptim_hour_subset.rds
* Estimated mean KDE bandwidth and temporal decay parameter when using the KDE bandwiddth mean and optimising over all individuals simultanously (as a csv)
  * memory_params_ALLoptim_hour_subset_2024-02-22.csv
* Data file that contains random steps, sampled covariate values and KDE density at the end of each used and random step (when using memory parameters that were estimated for each individual separately)
  * buffalo_popn_GvM_KDEmem_ID_hour_subset_10rs_2024-02-20.csv
* Data file that contains random steps, sampled covariate values and KDE density at the end of each used and random step (when using memory parameters that were estimated for all individuals simultaneously)
  * buffalo_popn_GvM_KDEmem_allOPTIM_hour_subset_10rs_2024-02-20.csv



## Script 2b - DynSSF_2b_memory_parameter_estimation_loc_subset.Rmd ##

This script estimates the most likely parameters for spatial (using Kernel Density Estimation - KDE) memory with a temporal decay (negative exponential) component. It subsets the previous locations by a fixed number of locations.

**Inputs**

* Data file that contains random steps and sampled covariate values at the end of each used and random step
  * buffalo_parametric_popn_covs_GvM_10rs_2024-02-18.csv

**Outputs**

* Estimated temporal decay parameters when optimising over each individual separately (as a list of optim objects)
  * temporal_decay_param_loc_list.rds
* Estimated KDE bandwidth and temporal decay parameters when optimising over each individual separately (as a csv)
  * memory_params_KDE_exp_decay_2024-02-20.csv
* Estimated temporal decay parameter when using the KDE bandwiddth mean and optimising over all individuals simultanously (as an optim object)
  * optim_space_time_Gamma_param_ALLoptim.rds
* Estimated mean KDE bandwidth and temporal decay parameter when using the KDE bandwiddth mean and optimising over all individuals simultanously (as a csv)
  * memory_params_ALLoptim_2024-02-22.csv
* Data file that contains random steps, sampled covariate values and KDE density at the end of each used and random step (when using memory parameters that were estimated for each individual separately)
  * buffalo_popn_GvM_KDEmem_ID_10rs_", Sys.Date(), ".csv (not used)
* Data file that contains random steps, sampled covariate values and KDE density at the end of each used and random step (when using memory parameters that were estimated for all individuals simultaneously)
  * buffalo_popn_GvM_covs_ST_KDEmem1000_allOPTIM_10rs_2024-02-05.csv



## Script 3 - DynSSF_3_daily_model_fit.Rmd ##

This script fits step selection models with varying numbers of harmonic terms. We used the memory parameters estimated with the location subset when optimisaing over all individuals simultaneously.

**Inputs**

* Data file that contains random steps, sampled covariate values and KDE density at the end of each used and random step (when using memory parameters that were estimated for all individuals simultaneously)
  * buffalo_popn_GvM_covs_ST_KDEmem1000_allOPTIM_10rs_2024-02-05.csv

**Outputs**

* Fitted model objects (0, 1, 2 or 3 pairs of harmonics)
  * model_twostep_0p_harms_dry.rds
  * model_twostep_1p_harms_dry.rds
  * model_twostep_2p_harms_dry.rds
  * model_twostep_3p_harms_dry.rds
* Tables of hourly coefficient values from fitted models 
  * TwoStep_0pDaily_coefs_dry_2024-02-20.csv
  * TwoStep_1pDaily_coefs_dry_2024-02-20.csv
  * TwoStep_2pDaily_coefs_dry_2024-02-20.csv
  * TwoStep_3pDaily_coefs_dry_2024-02-20.csv



## Script 4 - DynSSF_4_simulations.Rmd ##

This script takes the hourly coefficients of the fitted models and simulates (dynamic) animal movement trajectories, including memory with a warm-up period.

**Inputs**

* Tables of hourly coefficient values from fitted models (change the input file to change which model's coefficients are used to simulate trajectories) 
  * TwoStep_0pDaily_coefs_dry_2024-02-20.csv
  * TwoStep_1pDaily_coefs_dry_2024-02-20.csv
  * TwoStep_2pDaily_coefs_dry_2024-02-20.csv
  * TwoStep_3pDaily_coefs_dry_2024-02-20.csv
* Memory parameters from population level estimation
  * memory_params_ALLoptim_2024-02-05.csv
* Buffalo starting locations to start the simulations at
  * buffalo_starting_locations_13inds_2024-02-06.csv
* Spatial covariate rasters
  * ndvi_GEE_projected_watermask20230207.tif
  * canopy_cover.tif
  * veg_herby.tif
  * slope_raster.tif

**Outputs**

* Simulated data (animal movement trajectories)
In the 'simulated trajectories' folder (an example for a single starting location (corresponding to ID 2005, with a memory period of 500 locations, memory delay of 24 hours, 10 trajectories, 50 proposed steps at each time-point and 3000 steps in the entire trajectory)
  * 0p_id_2005_GvM_memALLoptim500_24_10ind_50ch_3000stp_2024-02-05.csv
  * ... other simulated trajectories



## Script 5 - DynSSF_5_trajectory_validation.Rmd ##

This script takes animal movement trajectories (observed and simulated), and calculates summary statistics which can be used to compare between the observed data and simulated data

**Inputs**

* Simulated data (animal movement trajectories)
In the 'simulated trajectories' folder (an example for a single starting location (corresponding to ID 2005, with a memory period of 500 locations, memory delay of 24 hours, 10 trajectories, 50 proposed steps at each time-point and 3000 steps in the entire trajectory)
  * 0p_id_2005_GvM_memALLoptim500_24_10ind_50ch_3000stp_2024-02-05.csv
  * ... other simulated trajectories
* Spatial covariate rasters
  * ndvi_GEE_projected_watermask20230207.tif
  * canopy_cover.tif
  * veg_herby.tif
  * slope_raster.tif
* Data file of buffalo GPS data that contains random steps, sampled covariate values and KDE density at the end of each used and random step (when using memory parameters that were estimated for all individuals simultaneously)
  * buffalo_popn_GvM_covs_ST_KDEmem1000_allOPTIM_10rs_2024-02-05.csv

**Outputs**

* Hourly summary statistic values for observed and simulated data (with number of harmonics depending on which trajectories are read in)
  * buffalo_summaries_hourly_habitat_2024-02-07.csv
  * sim_0p_memALL_summaries_hourly_habitat_2024-02-07.csv
  * sim_1p_memALL_summaries_hourly_habitat_2024-02-07.csv
  * sim_2p_memALL_summaries_hourly_habitat_2024-02-07.csv
  * sim_3p_memALL_summaries_hourly_habitat_2024-02-07.csv
* Entire trajectory summary statistic values for observed and simulated data (with number of harmonics depending on which trajectories are read in)
  * buffalo_summary_statistics_df_2024-02-07.csv
  * sim_0p_memALL_daily_summary_statistics_df_2024-02-07.csv
  * sim_1p_memALL_daily_summary_statistics_df_2024-02-07.csv
  * sim_2p_memALL_daily_summary_statistics_df_2024-02-07.csv
  * sim_3p_memALL_daily_summary_statistics_df_2024-02-07.csv



## Script 6 - DynSSF_6_comparing_summaries.Rmd ##

This script compares the summary statistics (hourly and entire-trajectory) between the observed data and simulated data

**Inputs**

* Hourly summary statistic values for observed and simulated data (with number of harmonics depending on which trajectories are read in)
  * buffalo_summaries_hourly_habitat_2024-02-07.csv
  * sim_0p_memALL_summaries_hourly_habitat_2024-02-07.csv
  * sim_1p_memALL_summaries_hourly_habitat_2024-02-07.csv
  * sim_2p_memALL_summaries_hourly_habitat_2024-02-07.csv
  * sim_3p_memALL_summaries_hourly_habitat_2024-02-07.csv
* Entire trajectory summary statistic values for observed and simulated data (with number of harmonics depending on which trajectories are read in)
  * buffalo_summary_statistics_df_2024-02-07.csv
  * sim_0p_memALL_daily_summary_statistics_df_2024-02-07.csv
  * sim_1p_memALL_daily_summary_statistics_df_2024-02-07.csv
  * sim_2p_memALL_daily_summary_statistics_df_2024-02-07.csv
  * sim_3p_memALL_daily_summary_statistics_df_2024-02-07.csv

**Outputs**

* Plots only

For more details, don't hesitate to get in contact.
