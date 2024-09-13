#### Script to read in the simulated trajectories on the HPC 
# The files in total are massive (~50 GB per model for 100,000 simulated trajectories),
# so it's easier to keep everything on the HPC

# Here we assess how many simulated trajectories are required to get a stable map of 
# predicted distribution, which is the number of simulated locations in each cell

library(tidyverse)
packages <- c("amt", "lubridate", "terra", "tictoc", "beepr", "ks", 
              "adehabitatLT", "adehabitatHR", "ggpubr", "patchwork",
              "sf", "spatstat", "maptools")
walk(packages, require, character.only = T)

# setwd("Z:/ssa_simulations")

# subset rasters
ndvi_stack <- rast("subset_landscape_sims/ndvi_20kmx20km_cropped.tif")
canopy <- rast("subset_landscape_sims/canopy_cover_20kmx20km_cropped.tif")
herby <- rast("subset_landscape_sims/veg_herby_20kmx20km_cropped.tif")
slope <- rast("subset_landscape_sims/slope_20kmx20km_cropped.tif")

ndvi_2018_dry <- ndvi_stack[[8:10]]
ndvi_2019_dry <- ndvi_stack[[19:21]]
ndvi_late_dry <- terra::mean(c(ndvi_2018_dry, ndvi_2019_dry))
names(ndvi_late_dry) <- "ndvi_late_dry"

plot(ndvi_late_dry)
plot(canopy)
plot(herby)
plot(slope)


# Import the simulated trajectories ---------------------------------------

# read in multiple csv files with similar filenames and bind them together
sim_data_full_list <- 
  list.files("subset_NOmemFIT_landscape_sims/2p_sims/",
             pattern = "*.csv", full.names = T)

sim_data_all <- grep("2p", sim_data_full_list, value = T) %>% 
  map_dfr(read_csv, show_col_types = FALSE)

str(sim_data_all)
sim_data_ids <- unique(sim_data_all$id)



# Create steps with step information (amt could be used here as well) --------

sim_data_all <- sim_data_all %>% mutate(
  
  x1_ = x_,
  y1_ = y_,
  x2_ = lead(x1_, n = 1, default = NA),
  y2_ = lead(y1_, n = 1, default = NA),
  t2_ = lead(t_, n = 1, default = NA),
  hour_t1 = ifelse(lubridate::hour(t_) == 0, 24, lubridate::hour(t_)),
  yday_t1 = lubridate::yday(t_),
  hour_t2 = ifelse(lubridate::hour(t2_) == 0, 24, lubridate::hour(t2_)),
  yday_t2 = lubridate::yday(t2_),
  
  sl = c(sqrt(diff(y_)^2 + diff(x_)^2), NA),
  bearing = c(atan2(diff(y_), diff(x_)), NA),
  ta = c(NA, ifelse(
    diff(bearing) > pi, diff(bearing)-(2*pi), ifelse(
      diff(bearing) < -pi, diff(bearing)+(2*pi), diff(bearing))))
  
)


# Assess the stability of the simulated predictions -----------------------

# number of resamples for each sample size
n_iterations <- 10

# template raster
template_raster <- rast(resolution = 50, extent = ext(ndvi_late_dry))
names(template_raster) <- "template_raster" 

# number of random cells distributed in the landscape - we assess the standard 
# deviation of predicted values (number of simulated locations) in each of these 
# cells from subsamples of the simulated trajectories

n_cells <- 1000
# keep cells constant between different samples of points
set.seed(123)
random_cells <- sample(1:ncell(ndvi_late_dry), n_cells, replace = FALSE)
random_cells_coords <- terra::xyFromCell(ndvi_late_dry, random_cells)
plot(random_cells_coords)

# for local testing
# n_samples <- seq(100, 1000, 100)
# for HPC
# create a sequence of the number of simulated trajectories used to create 
# the prediction maps, with more samples at the lower end
n_samples <- c(seq(10, 100, 10), seq(100, 1000, 100), seq(1000, 100000, 1000))



# Sample simulated trajectories and assess convergence --------------------

cell_summaries_list <- vector(mode = "list", length(n_samples))

tic("total time for all sims")

for(i in 1:length(n_samples)) {
  
  tic(paste0("iteration: ", n_samples[i], " sims"))
  
  # Initialize an empty data frame with the specified number of rows
  cell_values_df <- data.frame(matrix(ncol = 0, nrow = n_cells))
  
  for(j in 1:n_iterations) {
    
    # set the template raster values to 0
    template_raster[] <- 0
    
    # take a subset of the trajectories corresponding to the n_samples vector
    random_ids <- sample(1:length(sim_data_ids), n_samples[i], replace = FALSE)
    sim_data_random_subset <- sim_data_all[which(sim_data_all$id %in% sim_data_ids[random_ids]),]
    
    # create a matrix of the x and y coordinates
    coords_random_subset <- matrix(c(sim_data_random_subset$x_, 
                                     sim_data_random_subset$y_), 
                                   ncol = 2)
    
    # Rasterizing with points - using a matrix
    freq_points_2p_subset <- terra::rasterize(coords_random_subset, template_raster, fun = "sum")
    png(paste0("stability_plots/2p_stability/2p_subset_freq_points_", n_samples[i], "_samples.png"))
    plot(freq_points_2p_subset)
    dev.off()
  
    #normalizing the raster from 0 to 1
    freq_points_2p_subset_stretch <- terra::stretch(freq_points_2p_subset, minv = 0, maxv = 1)
    
    # Add the new row to the data frame
    cell_values <- as.vector(terra::extract(freq_points_2p_subset_stretch, random_cells_coords))
    cell_values_df[,j] <- cell_values
    
  }
  
  # calculate summaries across all of the cells
  cell_summaries_list[[i]] <- data.frame(n_sims = n_samples[i],
                                         mean = apply(cell_values_df, 1, mean, na.rm = TRUE),
                                         sd = apply(cell_values_df, 1, sd, na.rm = TRUE),
                                         var = apply(cell_values_df, 1, var, na.rm = TRUE),
                                         q025 = apply(cell_values_df, 1, quantile, probs = 0.025, na.rm = TRUE),
                                         q975 = apply(cell_values_df, 1, quantile, probs = 0.975, na.rm = TRUE))
  
  toc()
  
  cell_summaries_all <- do.call(rbind, cell_summaries_list)
  saveRDS(cell_summaries_all, "2p_cell_summaries_all_10iter_hpc.rds")
    
}

toc()



ggplot() +
  geom_jitter(data = cell_summaries_all, 
              aes(x = n_sims, y = sd), 
              na.rm = TRUE, width = 10, height = 0, 
              alpha = 0.1, size = 0.5) +
  scale_x_continuous("Number of simulated trajectories") +
  scale_y_continuous("Standard deviation") +
  ggtitle("Convergence of the simulated predictions - 2p") +
  theme_bw()

ggsave(paste0("stability_plots/2p_cell_summaries_all_10iter_hpc.png"),
       width=150, height=90, units="mm", dpi = 1000)

