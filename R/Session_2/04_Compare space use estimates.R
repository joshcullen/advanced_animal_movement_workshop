
#################################################
### Compare space-use estimates among methods ###
#################################################

library(tidyverse)
library(lubridate)
library(amt)
library(move)
library(sf)
library(rnaturalearth)
library(MetBrewer)
library(units)

source("R/utils.R")



#### Load the model results from each method ####

load("processed_data/KDE_fits.RData")
load("processed_data/dBBMM_fits.RData")
load("processed_data/AKDE_contours.RData")

# Change object to more informative names
kde_href <- dat_id_kde_href2; rm(dat_id_kde_href2)
kde_hpi <- dat_id_kde_hpi2; rm(dat_id_kde_hpi2)
dbbmm <- contours2; rm(contours2)

dat <- read_csv('processed_data/cleaned_tracks.csv') |> 
  mutate(id = as.character(id)) |>  #convert to char. to fix potential handling problems
  add_trans_coords(coords = c('lon','lat'), proj = 4326, new_proj = 32736)



#### Wrangle model results to match up properly ####

# Rename other output
kde_href <- kde_href |>
  mutate(method = 'KDE_href') |>
  dplyr::select(id, level, method, geometry)
kde_hpi <- kde_hpi |>
  mutate(method = 'KDE_hpi') |>
  dplyr::select(id, level, method, geometry)
dbbmm <- dbbmm |>
  mutate(method = 'dBBMM') |>
  dplyr::select(id, level, method, geometry)
akde_sf <- akde_sf |> 
  mutate(method = "AKDE",
         level = ifelse(level == "50%", 0.5, 0.95)) |> 
  rename(id_burst = id,
         id = id_orig) |> 
  dplyr::select(id, id_burst, level, method, geometry)




#### Visually compare UDs among methods ####

africa <- ne_countries(scale = 10, continent = c("Africa"), returnclass = "sf") |>
  st_transform(crs = 32736) |>
  dplyr::select(-level)

ud.all <- rbind(kde_href, kde_hpi, dbbmm, akde_sf[,-2])


# Show all IDs, methods, and levels
ggplot() +
  geom_sf(data = africa) +
  geom_path(data = dat, aes(x, y, group = id), size = 0.5, alpha = 0.5) +
  geom_sf(data = ud.all, aes(color = method), fill = "transparent", size = 0.5) +
  scale_color_met_d('Egypt') +
  theme_bw() +
  coord_sf(xlim = c(min(dat$x) - 5000, max(dat$x) + 5000),
           ylim = c(min(dat$y) - 5000, max(dat$y) + 5000)) +
  facet_grid(id ~ level)
#large disparity among methods (especially those estimating 'range' vs 'occurrence' distribution)


# Zoom in on ID 6471
ggplot() +
  geom_path(data = dat %>%
              filter(id == 6471), aes(x, y), size = 0.5, alpha = 0.15) +
  geom_sf(data = ud.all %>%
            filter(id == 6471), aes(color = method), fill = "transparent", size = 0.5) +
  scale_color_met_d('Egypt') +
  theme_bw() +
  facet_grid(id ~ level)
#AKDE and KDE_ref similar to each other (range distrib), whereas KDE_hpi and dBBMM more similar (occurrence distrib)






#### Compare estimated area of space-use ####

ud.all$area <- st_area(ud.all)

# Aggregate AKDE since split by bursts
ud.all2 <- ud.all |> 
  st_drop_geometry() |> 
  summarize(.by = c(id, level, method),
            area_full = sum(area))


# Compare by method and UD level
ggplot(ud.all2, aes(factor(level), set_units(area_full, km^2), color = method)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, position = position_dodge(0.55)) +
  geom_point(alpha = 0.7, position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.55)) +
  scale_color_met_d('Egypt') +
  labs(x = "Level", y = "Area") +
  theme_bw()



#-- While this is a bit of a simplification, this confirms the spatial patterns we saw on the map. However, this doesn't account for how the data are treated by each of the methods, for which AKDE is much more principled at estimating home ranges (assuming an animal is range-resident). To understand the relatively "tight" area traversed by an animal, especially during dispersal or migratory behavior, dBBMM provides a more principled approach to estimating what are essentially "confidence intervals" of trajectories (per a metaphor from John Fieberg). The method you choose should depend on you research questions and the properties of your data, rather than selecting a method you're familiar with to quickly measure space use. As shown here, the results can largely differ in the total area covered and the specific spatial regions covered. --# 
