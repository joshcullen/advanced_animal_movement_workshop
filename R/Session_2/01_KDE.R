
### Calculate kernel density estimates (KDE) ###

library(tidyverse)
library(amt)
library(sf)
library(terra)
library(tidyterra)
library(rnaturalearth)
library(plotly)
library(tictoc)
library(MetBrewer)
library(units)


###################
#### Load data ####
###################

dat <- read_csv('processed_data/cleaned_tracks.csv')

glimpse(dat)
summary(dat)




#######################################
#### Wrangle and prep data for KDE ####
#######################################

# Creat 'amt' track object
dat_track <- make_track(dat, lon, lat, date, crs = 4326, all_cols = TRUE)

# Convert to projected CRS
dat_track2 <- transform_coords(dat_track,
                               crs_to = 32736,  #UTM Zone 36 S; units in meters
                               crs_from = 4326)

# Create template raster for KDE
trast <- make_trast(dat_track2, res = 1000)  #resolution is 5 x 5 km




######################################
#### Calculate KDE across all IDs ####
######################################

#-- This may be of interest when generating pseudo-absences for an RSF --#
#-- Also could be helpful for measuring space use over short time intervals --#


## Href (reference bandwidth method); generally estimates larger areas
h_ref <- hr_kde_ref(dat_track2)
dat_kde_ref <- hr_kde(dat_track2, trast = trast, h = h_ref, levels = c(0.5, 0.95))
dat_kde_ref
plot(dat_kde_ref, col = c("red", "blue"))

# Check the estimated smoothing bandwidth
dat_kde_ref$h  #both values equal; denotes isotropic model

# Extract contours and raster for estimated UD
kde_href_contours <- hr_isopleths(dat_kde_ref)  #pulls the levels supplied to `hr_kde`
kde_href_rast <- dat_kde_ref$ud
kde_href_df <- as.data.frame(kde_href_rast, xy = TRUE)  #convert SpatRaster to data.frame

africa <- ne_countries(scale = 10, continent = c("Africa"), returnclass = "sf")


# Plot returned raster layer and contours
ggplot() +
  geom_sf(data = st_transform(africa, crs = 32736)) +
  geom_raster(data = kde_href_df, aes(x, y, fill = lyr.1)) +
  scale_fill_viridis_c(option = "rocket") +
  # geom_point(data = dat_track2, aes(x_, y_), alpha = 0.1, size = 0.01, color = "chartreuse") +
  geom_sf(data = kde_href_contours, aes(color = factor(level)), fill = NA, linewidth = 0.5) +
  scale_color_brewer("Level", palette = "Set1") +
  theme_bw() +
  coord_sf(xlim = ext(kde_href_rast)[1:2],
           ylim = ext(kde_href_rast)[3:4])


# Plot raster layer directly (via tidyterra::geom_spatraster)
ggplot() +
  geom_sf(data = st_transform(africa, crs = 32736)) +
  geom_spatraster(data = kde_href_rast) +  #replacing geom_raster()
  scale_fill_viridis_c(option = "rocket") +
  # geom_point(data = dat_track2, aes(x_, y_), alpha = 0.1, size = 0.01, color = "chartreuse") +
  geom_sf(data = kde_href_contours, aes(color = factor(level)), fill = NA, linewidth = 0.5) +
  scale_color_brewer("Level", palette = "Set1") +
  theme_bw() +
  coord_sf(xlim = ext(kde_href_rast)[1:2],
           ylim = ext(kde_href_rast)[3:4])


# Plot contours by themselves
ggplot() +
  geom_sf(data = st_transform(africa, crs = 32736)) +
  # geom_point(data = dat_track2, aes(x_, y_), alpha = 0.1, size = 0.01, color = "chartreuse") +
  geom_sf(data = kde_href_contours, aes(color = factor(level)), fill = NA, linewidth = 0.5) +
  scale_color_brewer("Level", palette = "Set1") +
  theme_bw() +
  coord_sf(xlim = ext(kde_href_rast)[1:2],
           ylim = ext(kde_href_rast)[3:4])







## Hpi (plug-in method); generally estimates smaller areas
h_pi <- hr_kde_pi(dat_track2, rescale = "xvar")
dat_kde_pi <- hr_kde(dat_track2, trast = trast, h = h_pi, levels = c(0.5, 0.95))
dat_kde_pi
plot(dat_kde_pi, col = c("red", "blue"))

# Extract contours and raster for estimated UD
kde_hpi_contours <- hr_isopleths(dat_kde_pi)  #pulls the levels supplied to `hr_kde`
kde_hpi_rast <- dat_kde_pi$ud


# Plot raster layer directly (via tidyterra::geom_spatraster)
ggplot() +
  geom_sf(data = st_transform(africa, crs = 32736)) +
  geom_spatraster(data = kde_hpi_rast) +
  scale_fill_viridis_c(option = "rocket") +
  # geom_point(data = dat_track2, aes(x_, y_), alpha = 0.1, size = 0.01, color = "chartreuse") +
  geom_sf(data = kde_hpi_contours, aes(color = factor(level)), fill = NA, linewidth = 0.5) +
  scale_color_brewer("Level", palette = "Set1") +
  theme_bw() +
  coord_sf(xlim = ext(kde_hpi_rast)[1:2],
           ylim = ext(kde_hpi_rast)[3:4])


# Plot contours by themselves
ggplot() +
  geom_sf(data = st_transform(africa, crs = 32736)) +
  # geom_point(data = dat_track2, aes(x_, y_), alpha = 0.1, size = 0.01, color = "chartreuse") +
  geom_sf(data = kde_hpi_contours, aes(color = factor(level)), fill = NA, linewidth = 0.5) +
  scale_color_brewer("Level", palette = "Set1") +
  theme_bw() +
  coord_sf(xlim = ext(kde_hpi_rast)[1:2],
           ylim = ext(kde_hpi_rast)[3:4])






##############################
#### Calculate KDE per ID ####
##############################

# If interested in installing a package that produces sounds when your code is finished running, this is my favorite (BRRR)
# https://github.com/brooke-watson/BRRR
# devtools::install_github("brooke-watson/BRRR")


## Href

# Need to turn data.frame into list to map hr_kde() function and then recombine
tic()
dat.id.kde.href <- dat.track |>
  split(.$id) |>  #split into list by ID
  map(~hr_kde(.x,  #map hr_kde() onto each list element (i.e., ID)
              trast = make_trast(.x, res = 0.5),  #define spatial res of raster
              h = hr_kde_ref(.x),  #define bandwidth for KDE
              levels = c(0.5, 0.95))
      ) |>
  map(hr_isopleths) |>  #extract contours from raster layers
  do.call(rbind, .)  #merge all contours into single `sf` object
toc()  #takes 28 sec to run
BRRR::skrrrahh('ross1')  #let me know that it's done running!

dat.id.kde.href <- dat.id.kde.href |>
  mutate(id = rownames(.), .before = level) |>
  mutate(id = str_replace(id, "\\..$", ""))  #remove decimal and extra number


ggplot() +
  geom_sf(data = brazil) +
  geom_path(data = dat, aes(lon, lat, group = id), alpha = 0.25, size = 0.3) +
  geom_sf(data = dat.id.kde.href, aes(color = factor(level)), fill = 'transparent', size = 0.75) +
  scale_color_met_d('Egypt') +
  theme_bw() +
  coord_sf(xlim = c(-44, -30), ylim = c(-9, 0)) +
  facet_wrap(~ id)





## Hpi

# Need to turn data.frame into list to map hr_kde() function and then recombine
dat.id.list <- dat.track |>
  split(.$id) |>
  map(~{.x |>
      mutate(pattern = ifelse(min(lon) > -34, 'Resident', 'Migratory'))
    })

id.pattern <- dat.id.list |>
  map(pluck, 'pattern') |>
  map(~{.x[1]})

res.ind <- which(id.pattern == 'Resident')
mig.ind <- which(id.pattern == 'Migratory')


# Run for Migratory IDs
tic()
dat.mig.kde.hpi <- dat.id.list[mig.ind] |>
  map(~hr_kde(.x,
              trast = make_trast(.x, res = 1),
              h = hr_kde_pi(.x, rescale = 'xvar'),
              levels = c(0.5, 0.95))
  ) |>
  map(hr_isopleths) |>
  do.call(rbind, .)
toc()  #takes 3 sec to run
BRRR::skrrrahh('khaled3')


# Run for Resident IDs
dat.res.kde.hpi <- dat.id.list[res.ind] |>
  map(~hr_kde(.x,
              trast = make_trast(.x, res = 0.1),
              h = hr_kde_pi(.x, rescale = 'xvar'),
              levels = c(0.5, 0.95))
  ) |>
  map(hr_isopleths) |>
  do.call(rbind, .)
BRRR::skrrrahh('liljon')

dat.id.kde.hpi <- rbind(dat.mig.kde.hpi, dat.res.kde.hpi) |>
  mutate(id = rownames(.), .before = level) |>
  mutate(id = str_replace(id, "\\..$", ""))  #remove decimal and extra number


ggplot() +
  geom_sf(data = brazil) +
  geom_path(data = dat, aes(lon, lat, group = id), alpha = 0.25, size = 0.3) +
  geom_sf(data = dat.id.kde.hpi, aes(color = factor(level)), fill = 'transparent', size = 0.75) +
  scale_color_met_d('Egypt') +
  theme_bw() +
  coord_sf(xlim = c(-42, -32), ylim = c(-8, -2)) +
  facet_wrap(~ id)







#### Compare KDE estimates between methods for bandwidth estimation ####

# Collate both ID-level data.frames together
dat.id.kde.href$method <- 'href'
dat.id.kde.hpi$method <- 'hpi'

dat.id.kde <- rbind(dat.id.kde.href, dat.id.kde.hpi)

ggplot(dat.id.kde, aes(factor(level), area, color = method)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, position = position_dodge(0.55)) +
  geom_point(alpha = 0.7, position = position_jitterdodge(jitter.width = 0.25, dodge.width = 0.55)) +
  scale_color_met_d("Hokusai3") +
  theme_bw()






#### Export datasets for easy loading ####

save(dat.id.kde.href, dat.id.kde.hpi, file = "Processed_data/KDE_fits.RData")
