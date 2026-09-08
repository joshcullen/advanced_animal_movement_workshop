
### Calculate dynamic Brownian Bridge Movement Model (dBBMM) ###

library(tidyverse)
library(move)
library(sf)
library(rnaturalearth)
library(tictoc)
library(terra)

source("R/utils.R")  #load in custom functions


###################
#### Load data ####
###################

dat <- read_csv('processed_data/Session_1/cleaned_tracks.csv')

glimpse(dat)
summary(dat)




#########################
#### Run dBBMM model ####
#########################

# Split data into list by ID and add coords in UTM
dat_list <- dat |>
  add_trans_coords(coords = c("lon","lat"), proj = 4326, new_proj = 32736) |>  #add UTM coords (labeled 'x' & 'y')
  split(~id)


dbbmm_list <- vector("list", length(dat_list))  #to store dBBMM results
contours <- vector("list", length(dat_list))  #to store resulting 50% and 95% UD contours


## Fit dBBMM separately by ID
for (i in seq_along(dat_list)) {
  message(paste("ID:", names(dat_list)[i]))  #print current ID

  # Create 'move' object
  dat_mov <- move(x = dat_list[[i]]$x, y = dat_list[[i]]$y, time = dat_list[[i]]$date, data = dat_list[[i]],
                  proj = "EPSG:32736",
                  animal = dat_list[[i]]$id)

  
  # Run dBBMM (this will take a little while to run)
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
lines(contours[[1]][2,], col = "black", bg = NA)

ud50 <- getVolumeUD(dbbmm_list[[1]])
ud50[ud50 > 0.50] <- NA
plot(ud50, main = paste0("ID ", names(dat_list)[1], ": 50% UD"))
lines(contours[[1]][1,], col = "black", bg = NA)


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
  geom_path(data = dat_list$`5605`, aes(x, y), linewidth = 0.5, alpha = 0.25) +
  geom_sf(data = contours2 |>
            filter(id == 5605), aes(color = level), fill = NA, linewidth = 0.75) +
  scale_color_brewer(palette = 'Set1') +
  labs(title = 'ID 5605') +
  theme_bw() +
  theme(panel.grid = element_blank())

#6469
ggplot() +
  geom_path(data = dat_list$`6469`, aes(x, y), linewidth = 0.5, alpha = 0.25) +
  geom_sf(data = contours2 |>
            filter(id == 6469), aes(color = level), fill = NA, linewidth = 0.75) +
  scale_color_brewer(palette = 'Set1') +
  labs(title = 'ID 6469') +
  theme_bw() +
  theme(panel.grid = element_blank())

#6471
ggplot() +
  geom_path(data = dat_list$`6471`, aes(x, y), linewidth = 0.5, alpha = 0.25) +
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







###############################################################
### Measure space use overlap among IDs and protected areas ###
###############################################################

# Convert dBBMMs to SpatRasters
rasters <- map(dbbmm_list, rast) |> 
  # Add ID to each layer name
  map2(.y = names(dat_list),
       .f = ~{
         names(.x) <- .y
         return(.x)
       })

# Modify raster to share common grid
rasters2 <- create_shared_grid(data = bind_rows(dat_list),
                               rasters = rasters,
                               crs = "epsg:32736",
                               res = 500,
                               ext = 0.3)


## Quanitfy overlap among UDs

# Calculate Volume of Intersection (VI) index; symmetric, so we only need one set of pairs
calc_ud_overlap(rasters2, index = "vi")

# Calculate Bhattacharyya's Affinity (BA) index; symmetric, so we only need one set of pairs
calc_ud_overlap(rasters2, index = "ba")



## Calculate overlap with National Parks
gl_pa <- st_read("raw_data/gltfca_protectedAreasDetailed.shp") |> 
  st_transform(crs = 32736)

ggplot() +
  geom_sf(data = africa) +
  geom_sf(data = gl_pa, color = "black", fill = NA, linewidth = 1) +
  # geom_path(data = dat, aes(lon, lat, group = id, color = factor(id)), alpha = 0.5, linewidth = 0.3) +
  # Plot 95% isopleths
  geom_sf(data = contours2 |> 
            filter(level == 0.95),
          aes(color = factor(id), fill = factor(id)), linewidth = 0.75, alpha = 0.5) +
  # Plot 50% isopleths
  # geom_sf(data = contours2 |> 
  #           filter(level == 0.5),
  #         aes(color = factor(id), fill = factor(id)), linewidth = 0.75, alpha = 0.7) +
  scale_color_brewer("ID", palette = "Dark2") +
  scale_fill_brewer("ID", palette = "Dark2") +
  theme_bw() +
  coord_sf(xlim = ext(rasters2)[1:2],
           ylim = ext(rasters2)[3:4])



# All UDs and NPs
NPs <- gl_pa |> 
  filter(Designatio == 'National Park')


# Calc overlap with NPs
for (i in 1:nrow(NPs)) {
  message(paste("Overlap with", NPs$Name[i],"\n"))
  
  print(calc_ud_overlap(rasters2, index = "feature", feature = NPs[i,]))
}







##########################################
#### Export datasets for easy loading ####
##########################################

save(contours2, file = "processed_data/Session_2/dBBMM_fits.RData")
