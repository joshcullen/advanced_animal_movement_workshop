
### Calculate dynamic Brownian Bridge Movement Model (dBBMM) ###

library(tidyverse)
library(move)
library(sf)
library(rnaturalearth)
library(tictoc)

source("R/utils.R")  #load in custom functions


###################
#### Load data ####
###################

dat <- read_csv('processed_data/cleaned_tracks.csv')

glimpse(dat)
summary(dat)



####################################################
#### Wrangle and prep data for dBBMM estimation ####
####################################################

# Split data into list by ID
dat_list <- dat |>
  add_trans_coords(coords = c("lon","lat"), proj = 4326, new_proj = 32736) |>  #add UTM coords (labeled 'x' & 'y')
  split(~id)





#########################
#### Run dBBMM model ####
#########################


dbbmm_list <- vector("list", length(dat_list))  #to store dBBMM results
contours <- vector("list", length(dat_list))  #to store resulting 50% and 95% UD contours


# Estimate separately by ID
for (i in seq_along(dat_list)) {
  message(paste("ID:", names(dat_list)[i]))  #print current ID

  # Create 'move' object
  dat_mov <- move(x = dat_list[[i]]$x, y = dat_list[[i]]$y, time = dat_list[[i]]$date, data = dat_list[[i]],
                  proj = "EPSG:32736",
                  animal = dat_list[[i]]$id)

  # Conditionally define extent; necessary for turtles that don't migrate
  # x.ext <- diff(dat.mov@bbox[1,])
  # rast.ext <- ifelse(x.ext < 100, 3, 0.3)


  ## Run dBBMM (this will take a little while to run)

  tic()
  dbbmm_list[[i]] <- brownian.bridge.dyn(object = dat_mov, raster = 500, location.error = 30,
                                       margin = 9, window.size = 29)
  toc()


  ## Extract 50 and 95% contours of space-use
  res <- raster2contour(dbbmm_list[[i]], levels = c(0.5, 0.95))

  contours[[i]] <- st_as_sf(res)

  if (st_geometry_type(contours[[i]])[1] == 'LINESTRING') {
    contours[[i]] <- st_cast(contours[[i]], "POLYGON")
  } else {
    contours[[i]] <- st_cast(contours[[i]], "MULTIPOLYGON") |>
      st_make_valid()  #fixes issue w/ negative areas being calculated
  }

}
BRRR::skrrrahh('liljon')
# took 2.5 min


# Quick plot of UDs at defined kernel volume (in raster form)
ud95 <- getVolumeUD(dbbmm_list[[1]])
ud95[ud95 > 0.95] <- NA
plot(ud95, main = paste0("ID ", names(dat_list)[1], ": 95% UD"))

ud50 <- getVolumeUD(dbbmm_list[[1]])
ud50[ud50 > 0.50] <- NA
plot(ud50, main = paste0("ID ", names(dat_list)[1], ": 50% UD"))


# Merge all UD contours together
names(contours) <- names(dat_list)
map(contours, st_area)

contours2 <- bind_rows(contours, .id = "id")




###############################
#### Let's viz the results ####
###############################

africa <- ne_countries(scale = 10, continent = c("Africa"), returnclass = "sf") |>
  st_transform(crs = 32736)

ggplot() +
  geom_sf(data = africa |> 
            dplyr::select(-level)) +
  geom_sf(data = contours2, aes(color = level), fill = NA, linewidth = 0.5) +
  scale_color_brewer(palette = 'Set1') +
  theme_bw() +
  coord_sf(xlim = st_bbox(contours2)[c(1,3)],
           ylim = st_bbox(contours2)[c(2,4)]) +
  facet_grid(level ~ id)
  

ggplot() +
  geom_sf(data = africa |> 
            dplyr::select(-level)) +
  geom_sf(data = contours2, aes(color = level), fill = NA, linewidth = 0.5) +
  scale_color_brewer(palette = 'Set1') +
  theme_bw() +
  coord_sf(xlim = st_bbox(contours2)[c(1,3)],
           ylim = st_bbox(contours2)[c(2,4)]) +
  facet_wrap(~ id)




##################################
#### Highlight IDs separately ####
##################################

#5605
ggplot() +
  geom_path(data = dat_list$`5605`, aes(x, y), linewidth = 0.5, alpha = 0.5) +
  geom_sf(data = contours2 |>
            filter(id == 5605), aes(color = level), fill = NA, linewidth = 0.75) +
  scale_color_brewer(palette = 'Set1') +
  labs(title = 'ID 5605') +
  theme_bw() +
  theme(panel.grid = element_blank())

#6470
ggplot() +
  geom_path(data = dat_list$`6470`, aes(x, y), linewidth = 0.5, alpha = 0.5) +
  geom_sf(data = contours2 |>
            filter(id == 6470), aes(color = level), fill = NA, linewidth = 0.75) +
  scale_color_brewer(palette = 'Set1') +
  labs(title = 'ID 6470') +
  theme_bw() +
  theme(panel.grid = element_blank())

#6471
ggplot() +
  geom_path(data = dat_list$`6471`, aes(x, y), linewidth = 0.5, alpha = 0.5) +
  geom_sf(data = contours2 |>
            filter(id == 6471), aes(color = level), fill = NA, linewidth = 0.75) +
  scale_color_brewer(palette = 'Set1') +
  labs(title = 'ID 6471') +
  theme_bw() +
  theme(panel.grid = element_blank())




# Plot all tracks w/ UD contours
ggplot() +
  geom_sf(data = africa |> 
            dplyr::select(-level)) +
  geom_path(data = bind_rows(dat_list), aes(x, y, color = factor(id), group = id), linewidth = 0.5, alpha = 0.5) +
  geom_sf(data = contours2, fill = NA, linewidth = 0.5, color = 'black') +
  scale_color_brewer("ID", palette = "Dark2") +
  theme_bw() +
  coord_sf(xlim = st_bbox(contours2)[c(1,3)],
           ylim = st_bbox(contours2)[c(2,4)]) +
  facet_wrap(~ level, ncol = 2)





##########################################
#### Export datasets for easy loading ####
##########################################

save(contours2, file = "processed_data/dBBMM_fits.RData")
