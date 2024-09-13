#### Script to read in the simulated trajectories on the HPC 
# The files in total are massive (~50 GB per model for 100,000 simulated trajectories),
# so it's easier to keep everything on the HPC

# Here we show an example for the model with 2 pairs of harmonics (2p)
# To run with other models read in the files with the appropriate model name

library(tidyverse)
packages <- c("lubridate", "terra")
walk(packages, require, character.only = T)


# Import the simulated trajectories ---------------------------------------

# read in multiple csv files with similar filenames and bind them together
sim_data_full_list <- 
  list.files("subset_landscape_sims/2p_sims/",
             pattern = "*.csv", full.names = T)

sim_data_all <- grep("2p", sim_data_full_list, value = T) %>% 
  map_dfr(read_csv, show_col_types = FALSE) 

str(sim_data_all)

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


# If using connectivity approaches, then it might be helpful to split the track
# up into bursts and exclude the locations that jump from one side of the extent
# to the other, due to the wrapped landscape. As we are only counting the number
# of simulated locations in each cell, we don't need to do that, but the function
# is included below for if it is required.

# Identify jumps by checking if the distance is greater than the threshold
# jumps <- which(sim_data_all$sl > 10000)
# jumps <- c(1, jumps)
# 
# bursts <- rep(NA, nrow(sim_data_all))
# 
# for(i in 1:(length(jumps)-1)) {
#   bursts[jumps[i]:jumps[i+1]] <- i
# }
# 
# sim_data_all <- sim_data_all %>% mutate(burst = bursts)
# 
# # remove any sections that have less than 3 points
# sim_data_all <- sim_data_all %>%
#   group_by(burst) %>%
#   filter(n() > 2) %>%
#   ungroup()

# sim_data_all_nested <- sim_data_all %>% arrange(id) %>% nest



# Count the number of simulated locations in each cell (using terra::rasterize) --------------------

raster_resolution <- 25

# template raster
template_raster <- rast(resolution = raster_resolution, extent = ext(ndvi_dry))
names(template_raster) <- "template_raster" 
template_raster[] <- 0

# create a mask layer to filter out NAs
mask_layer <- terra::resample(ndvi_dry, template_raster)


# Rasterizing with all simulated locations -------------------------------------------

coords_all <- matrix(c(sim_data_all$x_, 
                       sim_data_all$y_), 
                     ncol = 2)

# Rasterizing using a matrix
freq_points_2p_subset <- terra::rasterize(coords_all, template_raster, fun = "sum")

# to save raster images
# png(paste0("2p_subset_freq_points_res_", raster_resolutions[k], "_", Sys.Date(), ".png"))
# plot(freq_points_2p_subset)
# dev.off()

writeRaster(freq_points_2p_subset, 
            paste0("2p_subset_freq_points_res_", raster_resolution, "m_", Sys.Date(), ".tif", sep = ""),
            overwrite = TRUE)

coords_all <- NULL
  
  
  
# Rasterizing with hourly points ------------------------------------------

for(i in 1:24) {
  
  sim_data_hour <- sim_data_all %>% filter(hour == i)
  
  coords_hour <- matrix(c(sim_data_hour$x_, 
                          sim_data_hour$y_), 
                        ncol = 2)
  
  # Rasterizing with points only - using a matrix
  freq_points_2p_subset <- terra::rasterize(coords_hour, template_raster, fun = "sum")
  
  # to save raster images
  # png(paste0("2p_subset_freq_points_hour_", i, "_res_", raster_resolutions[k], "_", Sys.Date(), ".png"))
  # plot(freq_points_2p_subset)
  # dev.off()
  
  writeRaster(freq_points_2p_subset, 
              paste0("2p_subset_freq_points_hour_", i, "_res_", raster_resolution, "_", Sys.Date(), ".tif", sep = ""),
              overwrite = TRUE)
  
}
