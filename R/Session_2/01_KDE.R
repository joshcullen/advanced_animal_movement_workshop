
### Calculate kernel density estimates (KDE) ###

library(tidyverse)
library(amt)
library(sf)
library(terra)
library(tidyterra)
library(rnaturalearth)
library(tictoc)
library(units)


###################
#### Load data ####
###################

dat <- read_csv('processed_data/cleaned_tracks.csv') |> 
  mutate(id = as.character(id))  #convert to char. to fix potential handling problems

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
trast <- make_trast(dat_track2, res = 1000)  #resolution is 1 x 1 km




######################################
#### Calculate KDE across all IDs ####
######################################

#-- This may be of interest when generating pseudo-absences for an RSF --#
#-- Also could be helpful for measuring space use over short time intervals --#


## Href (reference bandwidth method); generally estimates larger areas
h_ref <- hr_kde_ref(dat_track2)
dat_kde_href <- hr_kde(dat_track2, trast = trast, h = h_ref, levels = c(0.5, 0.95))
dat_kde_href
plot(dat_kde_href, col = c("red", "blue"))

# Check the estimated smoothing bandwidth
dat_kde_href$h  #both values equal; denotes isotropic model

# Extract contours and raster for estimated UD
kde_href_contours <- hr_isopleths(dat_kde_href)  #pulls the levels supplied to `hr_kde`
kde_href_rast <- dat_kde_href$ud
kde_href_df <- as.data.frame(kde_href_rast, xy = TRUE)  #convert SpatRaster to data.frame

africa <- ne_countries(scale = 10, continent = c("Africa"), returnclass = "sf")


# Plot returned raster layer and contours
ggplot() +
  geom_sf(data = st_transform(africa, crs = 32736)) +
  geom_raster(data = kde_href_df, aes(x, y, fill = lyr.1), alpha = 0.7) +
  scale_fill_viridis_c("Density", option = "rocket") +
  # geom_point(data = dat_track2, aes(x_, y_), alpha = 0.1, size = 0.01, color = "chartreuse") +
  geom_sf(data = kde_href_contours, aes(color = factor(level)), fill = NA, linewidth = 0.5) +
  scale_color_brewer("Level", palette = "Set1") +
  theme_bw() +
  coord_sf(xlim = ext(kde_href_rast)[1:2],
           ylim = ext(kde_href_rast)[3:4])


# Plot raster layer directly (via tidyterra::geom_spatraster)
ggplot() +
  geom_sf(data = st_transform(africa, crs = 32736)) +
  geom_spatraster(data = kde_href_rast, alpha = 0.7) +  #replacing geom_raster()
  scale_fill_viridis_c("Density", option = "rocket") +
  # geom_point(data = dat_track2, aes(x_, y_), alpha = 0.1, size = 0.01, color = "chartreuse") +
  geom_sf(data = kde_href_contours, aes(color = factor(level)), fill = NA, linewidth = 0.5) +
  scale_color_brewer("Level", palette = "Set1") +
  theme_bw() +
  coord_sf(xlim = ext(kde_href_rast)[1:2],
           ylim = ext(kde_href_rast)[3:4])


# Plot contours by themselves
ggplot() +
  geom_sf(data = st_transform(africa, crs = 32736)) +
  geom_point(data = dat_track2, aes(x_, y_), alpha = 0.1, size = 0.01, color = "grey20") +
  geom_sf(data = kde_href_contours, aes(color = factor(level)), fill = NA, linewidth = c(0.5, 1)) +
  scale_color_brewer("Level", palette = "Set1") +
  theme_bw() +
  coord_sf(xlim = ext(kde_href_rast)[1:2],
           ylim = ext(kde_href_rast)[3:4]) +
  guides(color = guide_legend(override.aes = list(linewidth = 0.75)))







## Hpi (plug-in method); generally estimates smaller areas
h_pi <- hr_kde_pi(dat_track2, rescale = "xvar")
dat_kde_hpi <- hr_kde(dat_track2, trast = trast, h = h_pi, levels = c(0.5, 0.95))
dat_kde_hpi
plot(dat_kde_hpi, col = c("red", "blue"))

# Extract contours and raster for estimated UD
kde_hpi_contours <- hr_isopleths(dat_kde_hpi)  #pulls the levels supplied to `hr_kde`
kde_hpi_rast <- dat_kde_hpi$ud


# Plot raster layer directly (via tidyterra::geom_spatraster)
ggplot() +
  geom_sf(data = st_transform(africa, crs = 32736)) +
  geom_spatraster(data = kde_hpi_rast, alpha = 0.75) +
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
  geom_point(data = dat_track2, aes(x_, y_), alpha = 0.1, size = 0.01, color = "grey20") +
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
# pak::pak("brooke-watson/BRRR")


## Href

# Need to turn data.frame into list to map hr_kde() function and then recombine
# Follows split-apply-combine "tidy" data science principle

tic()
dat_id_kde_href <- dat_track2 |>
  split(~id) |>  #split into list by ID
  map(~hr_kde(.x,  #map hr_kde() onto each list element (i.e., ID)
              # trast = make_trast(.x, res = 1000),  #define spatial res of raster (1000 m)
              trast = trast,
              h = hr_kde_ref(.x),  #define reference bandwidth for KDE
              levels = c(0.5, 0.95))
      )
toc()  #takes 1 sec to run
BRRR::skrrrahh('ross1')  #let me know that it's done running!

# Extract contours
dat_id_kde_href2 <- dat_id_kde_href |>
  map(hr_isopleths) |>  #extract contours from raster layers
  bind_rows(.id = "id")  #merge all contours into single `sf` object



# Plotted separately by ID
ggplot() +
  geom_sf(data = africa) +
  geom_path(data = dat, aes(lon, lat, group = id), alpha = 0.25, linewidth = 0.3) +
  geom_sf(data = st_transform(dat_id_kde_href2, crs = 4326),
          aes(color = factor(level)), fill = NA, linewidth = 0.75) +
  scale_color_brewer("Level", palette = "Set1") +
  theme_bw() +
  coord_sf(xlim = range(dat$lon),
           ylim = range(dat$lat)) +
  facet_wrap(~ id)


# Plotted together
ggplot() +
  geom_sf(data = africa) +
  geom_path(data = dat, aes(lon, lat, group = id, color = factor(id)), alpha = 0.5, linewidth = 0.3) +
  # Plot 95% isopleths
  geom_sf(data = dat_id_kde_href2 |> 
            st_transform(crs = 4326) |> 
            filter(level == 0.95),
          aes(color = factor(id)), fill = NA, linewidth = 0.75) +
  # Plot 50% isopleths
  geom_sf(data = dat_id_kde_href2 |> 
            st_transform(crs = 4326) |> 
            filter(level == 0.5),
          aes(color = factor(id), fill = factor(id)), linewidth = 0.75, alpha = 0.7) +
  scale_color_brewer("ID", palette = "Dark2") +
  scale_fill_brewer("ID", palette = "Dark2") +
  theme_bw() +
  coord_sf(xlim = range(dat$lon),
           ylim = range(dat$lat) + c(0,0.25))  #expand upper ylim only to fully viz contours



## Hpi

tic()
dat_id_kde_hpi <- dat_track2 |>
  split(~id) |> 
  map(~hr_kde(.x,
              trast = make_trast(.x, res = 100),
              h = hr_kde_pi(.x, rescale = 'xvar'),
              levels = c(0.5, 0.95))
  ) 
toc()  #takes 37 sec to run
BRRR::skrrrahh('khaled3')

# Extract contours
dat_id_kde_hpi2 <- dat_id_kde_hpi|>
  map(hr_isopleths) |>
  bind_rows(.id = "id")





# Plotted separately by ID
ggplot() +
  geom_sf(data = africa) +
  geom_path(data = dat, aes(lon, lat, group = id), alpha = 0.25, linewidth = 0.3) +
  geom_sf(data = st_transform(dat_id_kde_hpi2, crs = 4326),
          aes(color = factor(level)), fill = NA, linewidth = 0.75) +
  scale_color_brewer("Level", palette = "Set1") +
  theme_bw() +
  coord_sf(xlim = range(dat$lon),
           ylim = range(dat$lat)) +
  facet_wrap(~ id)


# Plotted together
ggplot() +
  geom_sf(data = africa) +
  geom_path(data = dat, aes(lon, lat, group = id, color = factor(id)), alpha = 0.5, linewidth = 0.3) +
  # Plot 95% isopleths
  geom_sf(data = dat_id_kde_hpi2 |> 
            st_transform(crs = 4326) |> 
            filter(level == 0.95),
          aes(color = factor(id)), fill = NA, linewidth = 0.75) +
  # Plot 50% isopleths
  geom_sf(data = dat_id_kde_hpi2 |> 
            st_transform(crs = 4326) |> 
            filter(level == 0.5),
          aes(color = factor(id), fill = factor(id)), linewidth = 0.75, alpha = 0.7) +
  scale_color_brewer("ID", palette = "Dark2") +
  scale_fill_brewer("ID", palette = "Dark2") +
  theme_bw() +
  coord_sf(xlim = range(dat$lon),
           ylim = range(dat$lat))






######################################################################
### Compare KDE estimates between methods for bandwidth estimation ###
######################################################################

# Collate both ID-level data.frames together
dat_id_kde_href2$method <- 'href'
dat_id_kde_hpi2$method <- 'hpi'

dat_id_kde <- rbind(dat_id_kde_href2, dat_id_kde_hpi2)

ggplot(dat_id_kde, aes(factor(level), area, color = method)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, position = position_dodge(0.55), linewidth = 1) +
  geom_point(alpha = 0.7, position = position_dodge(width = 0.55), size = 3) +
  scale_color_manual("Bandwidth method", values = viridis::cividis(n = 2, end = 0.9)) +
  theme_bw(base_size = 14)
#h_ref produces much larger estimates of space use compared to h_pi



###############################################################
### Measure space use overlap among IDs and protected areas ###
###############################################################

# Refer to https://jmsigner.r-universe.dev/articles/amt/p2_hr.html for more examples 

#-- Could be generalized to assess overlap for same ID across seasons (or other time period) --#
#-- Could be used to evaluate ontogenetic changes in space use --#
#-- Could be used to assess diel patterns --#
#-- Could be used to assess UD overlap w/ spatial feature --#


### Quantify overlap among UDs

# For single pair
hr_overlap(dat_id_kde_href$`5605`, dat_id_kde_href$`6471`, type = 'hr', conditional = TRUE)

# For all pairs
hr_overlap(x = dat_id_kde_href, type = 'hr', which = "all", conditional = TRUE)
hr_overlap(x = dat_id_kde_href, type = 'vi', which = "all", conditional = TRUE)
hr_overlap(x = dat_id_kde_href, type = 'ba', which = "all", conditional = TRUE)




# Quantify overlap with GL protected areas
gl_pa <- st_read("raw_data/gltfca_protectedAreasDetailed.shp") |> 
  st_transform(crs = 32736)

ggplot() +
  geom_sf(data = africa) +
  geom_sf(data = gl_pa |> 
            st_transform(4326), color = "black", fill = NA, linewidth = 1) +
  geom_path(data = dat, aes(lon, lat, group = id, color = factor(id)), alpha = 0.5, linewidth = 0.3) +
  # Plot 95% isopleths
  geom_sf(data = dat_id_kde_href2 |> 
            st_transform(crs = 4326) |> 
            filter(level == 0.95),
          aes(color = factor(id)), fill = NA, linewidth = 0.75) +
  # Plot 50% isopleths
  geom_sf(data = dat_id_kde_href2 |> 
            st_transform(crs = 4326) |> 
            filter(level == 0.5),
          aes(color = factor(id), fill = factor(id)), linewidth = 0.75, alpha = 0.7) +
  scale_color_brewer("ID", palette = "Dark2") +
  scale_fill_brewer("ID", palette = "Dark2") +
  theme_bw() +
  coord_sf(xlim = range(dat$lon),
           ylim = range(dat$lat))

# One UD-feature combination at a time
hr_overlap_feature(dat_id_kde_href$`5605`,
                   gl_pa |> 
                     filter(Name == 'Limpopo'),
                   direction = "hr_with_feature")


# All UDs and NPs
NPs <- gl_pa |> 
  filter(Designatio == 'National Park')

for (i in seq_along(dat_id_kde_href)) {
  message(paste("ID:", dat_id_kde_href[[i]]$data$id[1]),"\n")
  
  for (j in 1:nrow(NPs)) {
    tmp <- hr_overlap_feature(dat_id_kde_href[[i]],
                              NPs[j,],
                              direction = "hr_with_feature",
                              feature_names = NPs[j,]$Name)
    print(tmp)  
  }
}






##########################################
#### Export datasets for easy loading ####
##########################################

save(dat_id_kde_href2, dat_id_kde_hpi2, file = "processed_data/KDE_fits.RData")
