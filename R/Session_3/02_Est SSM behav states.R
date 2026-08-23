
### Estimate behavioral state from SSM (continuous states) ###

library(aniMotum)
library(tidyverse)
library(bayesmove)
library(rnaturalearth)
library(sf)
library(tictoc)



###################
#### Load data ####
###################

### Load tracks

# Load tracks w/ bursts
dat <- read_csv("processed_data/Session_3/track_bursts.csv")

# Load regularized tracks (mean est.)
ssm_tracks <- read_csv("processed_data/Session_3/regularized_tracks.csv") |> 
  arrange(id, date) |>  #make sure data is properly sorted
  mutate(id = as.character(id))  #IDs are better handled as 'character'

# Load fitted SSM objects
load(file = "processed_data/Session_3/ssm_fits.RData")

# Load multiple imputations from fitted SSM
mi_tracks <- read_csv(file = "processed_data/Session_3/mi_tracks.csv")



### Load spatial layers
africa <- ne_countries(scale = 10, continent = c("Africa"), returnclass = "sf")
gl_pa <- st_read("raw_data/gltfca_protectedAreasDetailed.shp")






###############################################################
### Estimate behavioral states (as move persistence; gamma) ###
###############################################################

### From average tracks

## Individual move persistence model ('mpm') estimates behavioral states separately across IDs
tic()
rw_fit_g <- fit_mpm(rw_fit, what = "predicted", model = "mpm",
                    control = mpm_control(verbose = 1))
toc()  #took 1 min to fit
#if issues w/ model not converging, try changing time step of SSM (or optimizer for fit_mpm) and re-run

print(rw_fit_g)  #all models converged
plot(rw_fit_g)
# not really clear what's happening at this fine temporal scale


## Try fitting SSM at coarser time scale and re-evaluating
dat_sf <- dat |> 
  relocate(lon, .before = lat) |>  #fix column order for aniMotum::fit_ssm()
  mutate(lc = 'G', .after = date) |>  #need to specify "location class" (G = GPS)
  st_as_sf(coords = c('lon','lat'), crs = 4326, remove = FALSE) |>  #convert to spatial object
  st_transform(32736)  #convert to a projected ref system (e.g., UTM)

tic()
rw_fit2 <- fit_ssm(dat_sf, 
                   model = "rw", 
                   time.step = 8,  #8 hr
                   spdf = FALSE,
                   control = ssm_control(verbose = 1, tdist = "norm"))
toc()  #took 1 min to fit


print(rw_fit2)  #all indiv. models converged
summary(rw_fit2)

# Viz the fitted vs observed pts; blue = observed; gold = fitted
plot(rw_fit2, what = "predicted", type = 1, ask = TRUE)  #plot time series of coords
plot(rw_fit2, what = "predicted", type = 2, alpha = 0.1, ask = TRUE)  #plot maps of tracks

# Estimate move persistence at 8 hr scale
tic()
rw_fit_g8h <- fit_mpm(rw_fit2, what = "predicted", model = "mpm",
                      control = mpm_control(verbose = 1))
toc()  #took 8 sec to fit

print(rw_fit_g8h)  #all models converged
plot(rw_fit_g8h)





### Using previously fitted 'move persistence' model
plot(mp_fit, type = 3, normalise = TRUE)
plot(mp_fit, type = 4, normalise = TRUE)
# more variability than other 2-step approach

# Try re-fitting at coarser time scale
tic()
mp_fit2 <- fit_ssm(dat_sf, 
                   model = "mp", 
                   time.step = 8,  #8 hr
                   spdf = FALSE,
                   control = ssm_control(verbose = 1, tdist = "norm"))
toc()  #took 2 min to fit


print(mp_fit2)  #all indiv. models converged
summary(mp_fit2)

# Viz the fitted vs observed pts; blue = observed; gold = fitted
plot(mp_fit2, what = "predicted", type = 1, ask = TRUE)  #plot time series of coords
plot(mp_fit2, what = "predicted", type = 2, alpha = 0.1, ask = TRUE)  #plot maps of tracks

plot(mp_fit2, type = 3, normalise = TRUE)
plot(mp_fit2, type = 4, normalise = TRUE)
#behavioral patterns match up much better now at this coarser scale





###################################################
### Extract annotated tracks from model objects ###
###################################################

# Via 2-stage method
rw_behav <- join(ssm = rw_fit,
                 mpm = rw_fit_g,
                 what.ssm = "predicted",
                 normalise = TRUE)

# Via single stage 'mp' model
mp_behav <- grab(mp_fit2, what = "predicted", normalise = TRUE)



### Filter out large temporal gaps (using bursts)

burst_windows <- dat |>
  group_by(id, burst_id) |>
  summarize(
    start_time = min(date),
    end_time = max(date),
    .groups = "drop"
  ) |> 
  ungroup() |> 
  mutate(id = as.character(id))  #needs to match class for other df

# Filter fitted tracks (i.e., remove interpolated section during long gaps)
rw_behav2 <- rw_behav |> 
  inner_join(burst_windows,
             by = join_by(id, between(date, start_time, end_time)))

mp_behav2 <- mp_behav |> 
  inner_join(burst_windows,
             by = join_by(id, between(date, start_time, end_time)))




### Viz mapped tracks

# RW model
ggplot() +
  geom_sf(data = africa) +
  geom_path(data = rw_behav, aes(lon, lat, group = id), alpha = 0.5, linewidth = 0.25) +
  geom_point(data = rw_behav2, aes(lon, lat, color = g), alpha = 0.5, size = 0.5) +
  geom_sf(data = gl_pa, color = "black", fill = NA, linewidth = 0.25) +
  scale_color_viridis_c("Move\nPersistence", option = "viridis", limits = c(0,1)) +
  theme_bw() +
  coord_sf(xlim = range(dat$lon),
           ylim = range(dat$lat))

# MP model
ggplot() +
  geom_sf(data = africa) +
  geom_path(data = mp_behav, aes(lon, lat, group = id), alpha = 0.5, linewidth = 0.25) +
  geom_point(data = mp_behav2, aes(lon, lat, color = g), alpha = 0.5, size = 1) +
  geom_sf(data = gl_pa, color = "black", fill = NA, linewidth = 0.25) +
  scale_color_viridis_c("Move\nPersistence", option = "viridis", limits = c(0,1)) +
  theme_bw() +
  coord_sf(xlim = range(dat$lon),
           ylim = range(dat$lat))


# Interactive mapping
bayesmove::shiny_tracks(data = mp_behav2,
                        epsg = "+proj=utm +zone=36 +ellps=WGS84 +units=km +no_defs +south")


#-- Not explored in detail here, but time scale has large impact on ecological inferences. Care should be taken to address research questions when selecting a time scale and behavioral state method. For example, the behavioral stat estimates at the 1 hr time step did not appear to provide informative behavioral states from this model, but DID perform better at coarser scales for coarser behavioral patterns. However, custom Bayesian versions of this SSM may produce more informative results, such as through the inclusion of covariates. --#


###############################
### Export annotated tracks ###
###############################

write_csv(rw_behav2, "processed_data/Session_3/rw_behav_1h.csv")
write_csv(mp_behav2, "processed_data/Session_3/mp_behav_8h.csv")
