
### Prepare and regularize tracks ###

#track regularization isn't necessary if tracks are 1) highly regular, 2) only have several large time gaps, or 3) a continuous-time behavioral state model is being used

library(aniMotum)
library(tidyverse)
library(bayesmove)
library(rnaturalearth)
library(sf)
library(tictoc)

source("R/utils.R")



###################
#### Load data ####
###################

dat <- read_csv('processed_data/cleaned_tracks.csv') |> 
  arrange(id, date) |>  #make sure data is properly sorted
  mutate(id = as.character(id))  #IDs are better handled as 'character'

summary(dat)
glimpse(dat)



############################
#### Inspect time steps ####
############################

dat <- dat |>
  split(~id) |>  #split into list by ID
  purrr::map(~mutate(.x,
                     dt = difftime(c(date[-1], NA),
                                   date,
                                   units = "hours") |>  #calc time step (in hours)
                       as.numeric())
  ) |>
  bind_rows()  #combine back to a df


# Viz dt patterns
ggplot(dat, aes(date, dt)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~id, scales = "free_x")
#dt shows irregularity; will need to use multiple imputation or regularize tracks (unless using continuous-time model)
#ID 6469 has several very large gaps (>24 hrs) in early 2025

dat |> 
  pull(dt) |> 
  summary()
table(dat$dt) |> 
  sort(decreasing = TRUE) |> 
  head()
#1 hr is mean/median/mode
#we'll need to regularize all tracks




###########################################################
### Check that tracks have sufficient data for analysis ###
###########################################################

# Sample size
dat |>
  group_by(id) |>
  count()
#min of 8796 obs for ID 5605

# Duration
dat |>
  group_by(id) |>
  summarize(duration = last(date) - first(date))
#min is 376 days for ID 5605

#All tracks have sufficient data




##################################
### Segment tracks into bursts ###
##################################

# If "large" gaps exist in the tracks, you'll likely need to segment them into bursts to avoid excessive interpolation by the model when an animal wasn't observed

# Let's explore if we'd have enough data per burst if analyzing separately when dt > 24 hrs


# Define time threshold on which to split tracks into bursts
dt_thresh <- 24

# Define a new burst ID every time the threshold is exceeded
dat2 <- dat |>
  group_by(id) |>
  arrange(date, .by_group = TRUE) |>
  mutate(
    # lag(dt) looks at the gap BEFORE the current row
    # coalesce(..., FALSE) prevents NA values from breaking the cumsum
    # Adding 1 ensures burst IDs start at 1 (instead of 0)
    burst_id = 1 + cumsum(coalesce(lag(dt) > dt_thresh, FALSE)) # Increment burst ID whenever dt > threshold
  ) |>  
  ungroup() |> 
  mutate(burst_id = paste(id, burst_id, sep = "_"))  #create unique burst IDs for model fitting


# View burst summary 
dat2 |>
  group_by(burst_id) |>
  summarize(
    start_time = min(date),
    end_time = max(date),
    n = n(),
    .groups = "drop"
  )
#ID 6469 has most bursts (n = 5), and 3 of these have small N (N < 30)




# Let's increase the threshold to 36 hrs to limit the number of short bursts and avoid tossing out data

# Define time threshold on which to split tracks into bursts
dt_thresh <- 36

# Define a new burst ID every time the threshold is exceeded
dat2 <- dat |>
  group_by(id) |>
  arrange(date, .by_group = TRUE) |>
  mutate(burst_id = 1 + cumsum(coalesce(lag(dt) > dt_thresh, FALSE))) |> # Increment burst ID when dt > dt_thresh
  ungroup() |> 
  mutate(burst_id = paste(id, burst_id, sep = "_"))  #create unique burst IDs for model fitting


# View burst summary 
dat2 |>
  group_by(burst_id) |>
  summarize(
    start_time = min(date),
    end_time = max(date),
    n = n(),
    .groups = "drop"
  )


# Remove bursts w/ few points (e.g., n < 30)
dat3 <- dat2 |>
  group_by(id, burst_id) |>
  filter(n() >= 30) |>
  group_by(id) |>
  mutate(burst_id = dense_rank(burst_id)) |>  #ensure all bursts are consecutive and begin at 1
  ungroup() |>
  mutate(burst_id = paste(id, burst_id, sep = "_"))  #create unique burst IDs for model fitting










###############
### Fit SSM ###
###############

# These models use continuous-time formulations, which are well-suited to irregular time series
# The similar {bsam} package relied on discrete-time SSMs to estimate track locations
# Similar continuous-time models also available from {ctmm} and {crawl}

# 3 models to choose from:
#  - "rw": Random walk model is better suited to tracks w/ little error (e.g., GPS), but can overfit to noisy data, especially when time steps are small compared to specified time step
#  - "crw": Correlated random walk model is most commonly used since it accounts for autocorrelation/persistence in movements and therefore better handles these short/medium time steps. However, it can often produce nonsensical "looping" artefacts during long gaps or when animals are stationary for extended periods, as well as "move" fitted tracks away from accurate locations
#  - "mp": Move persistence model can possibly handle these time gaps better than the other models and also accounts for a behavioral movement process when estimating these locations (via a time-varying parameter). Like the CRW, may not always keep fitted tracks close to locations, even when using GPS data


# Prep data for analysis
dat_sf <- dat3 |> 
  rename(id_orig = id, id = burst_id) |>  #treat 'burst_id' as primary ID for model fitting
  relocate(lon, .before = lat) |>  #fix column order for aniMotum::fit_ssm()
  mutate(lc = 'G', .after = date) |>  #need to specify "location class" (G = GPS)
  st_as_sf(coords = c('lon','lat'), crs = 4326, remove = FALSE) |>  #convert to spatial object
  st_transform(32736)  #convert to a projected ref system (e.g., UTM)
#leaving in lon/lat will result in default projection being used (global Mercator; epsg 3395)



### Random walk model (RW)
tic()
rw_fit <- fit_ssm(dat_sf, 
                  model = "rw", 
                  time.step = 1,  #1 hr; set to NA if wanting locations at observed timestamps ONLY
                  spdf = FALSE,  #turn off pre-filtering obs
                  map = list(rho_o = factor(NA)),  #turn off est. of the obs. error corr. param; helps for GPS data
                  control = ssm_control(verbose = 1,  #setting to positive integer prints real-time param est.
                                        tdist = "norm"),  #Normal (Gaussian) err. distr. more appropriate for GPS
                  )
toc()  #took 2 min to fit


print(rw_fit)  #all indiv. models converged
summary(rw_fit)

# Viz the fitted vs observed pts; blue = observed; gold = fitted
plot(rw_fit, what = "predicted", type = 1, ask = TRUE)  #plot time series of coords
plot(rw_fit, what = "predicted", type = 2, alpha = 0.1, ask = TRUE)  #plot maps of tracks




### Correlated random walk model (CRW)
tic()
crw_fit <- fit_ssm(dat_sf, 
                  model = "crw", 
                  time.step = 1,  #1 hr
                  spdf = FALSE,  #turn off pre-filtering obs
                  map = list(rho_o = factor(NA)),  #turn off obs error corr
                  control = ssm_control(verbose = 1, tdist = "norm"))
toc()  #took 2.5 min to fit


print(crw_fit)  #all indiv. models converged
summary(crw_fit)

# Viz the fitted vs observed pts; blue = observed; gold = fitted
plot(crw_fit, what = "predicted", type = 1, ask = TRUE)  #plot time series of coords
plot(crw_fit, what = "predicted", type = 2, alpha = 0.1, ask = TRUE)  #plot maps of tracks




### Move persistence model (MP)
tic()
mp_fit <- fit_ssm(dat_sf, 
                  model = "mp", 
                  time.step = 1,  #1 hr
                  spdf = FALSE,  #turn off pre-filtering obs
                  map = list(rho_o = factor(NA)),  #turn off obs error corr
                  control = ssm_control(verbose = 1, tdist = "norm"))
toc()  #took 5 min to fit


print(mp_fit)  #all indiv. models converged
summary(mp_fit)

# Viz the fitted vs observed pts; blue = observed; gold = fitted
plot(mp_fit, what = "predicted", type = 1, ask = TRUE)  #plot time series of coords
plot(mp_fit, what = "predicted", type = 2, alpha = 0.1, ask = TRUE)  #plot maps of tracks



## Since we're using GPS data here and want fitted tracks to more-or-less stay close to these points, the Random Walk model may be the way to go. But let's compare by AIC

model_selection <- rbind(summary(rw_fit)[[1]],
                         summary(crw_fit)[[1]],
                         summary(mp_fit)[[1]]) |> 
  data.frame() |> 
  mutate(AICc = as.numeric(AICc))

model_selection |> 
  split(~Animal.id) |> 
  map(~{
    .x |> 
      mutate(dAICc = AICc - min(AICc)) |> 
      arrange(dAICc)
  })
#Using AIC, it appears that "mp" model is best for all tracks, but no real hard-and-fast rules
#Per Ian Jonsen: "Note that AIC statistics can be misleading for time-series models and should not be used as the sole criterion for preferring one model fit over another."
#Since the "rw" and "mp" models seem to most closely track the observed locs for our accurate GPS data, we'll use the "mp" model results; but possible to mix and match like what was done for the AKDE workflow


### Check goodness-of-fit (GOF) via one-step-ahead residuals
#WARNING: especially for large datasets, this can take a LONG time to run

# Calculate residuals
tic()
res_mp <- osar(mp_fit)
toc()  #took x min

# Viz plots of residuals
plot(res_mp, type = "ts")  #resids over time
plot(res_mp, type = "qq")  # QQ plot
plot(res_mp, type = "acf")  #autocorr plot




#####################################
### Extract fitted tracks and viz ###
#####################################

### Grab results and create data.frame
ssm_res <- grab(mp_fit, what = "predicted")  #if what = "fitted", returns locations at observed timestamps

# Add back the "original" ID (in addition to burst ID)
ssm_res <- ssm_res |> 
  mutate(id_orig = as.vector(str_match(id, "[0-9]+")),
         .after = id)

### Viz fitted tracks
africa <- ne_countries(scale = 10, continent = c("Africa"), returnclass = "sf")

# Plot fitted tracks only
ggplot() +
  geom_sf(data = africa) +
  #Plot trajectories by burst ("group"), but color trajectories by animal ("color")
  geom_path(data = ssm_res, aes(lon, lat, color = id_orig, group = id), linewidth = 0.25) +
  scale_color_brewer(palette = 'Dark2') +
  theme_bw() +
  coord_sf(xlim = range(ssm_res$lon),
           ylim = range(ssm_res$lat))

# Compare original and fitted tracks
ggplot() +
  geom_path(data = dat3 |> 
              rename(id_orig = id, id = burst_id),
            aes(lon, lat, group = id), color = "black", linewidth = 0.25) +  #observed locs
  geom_point(data = ssm_res, aes(lon, lat), color = "red", size = 0.1) +  #predicted locs @ regularized intervals
  theme_bw() +
  facet_wrap(~id_orig, scales = "free")
#fitted locs show high fidelity to original tracks


# For visualizing *many* tracks, use extension function from {ggforce}
for (i in 1:n_distinct(ssm_res$id_orig)) {
  print(
    ggplot() +
      geom_path(data = dat3 |> 
                  rename(id_orig = id, id = burst_id),
                aes(lon, lat, group = id), color = "black", linewidth = 0.25) +  #observed locs
      geom_point(data = ssm_res, aes(lon, lat), color = "red", size = 0.1) +  #predicted locs @ regularized intervals
      theme_bw() +
      ggforce::facet_wrap_paginate(~id_orig, scales = "free", nrow = 1, ncol = 1, page = i)
  )
}




# Viz on interactive map
ssm_res |> 
  select(-c(x, y)) |> 
  rename(x = lon, y = lat) |> 
  shiny_tracks(epsg = 4326)





#########################################
### Simulate tracks from fitted model ###
#########################################

# Useful if wanting to account for uncertainty in tracks due to location error and time series irregularity
# Essentially performs *multiple imputation* from posterior estimates of model parameters
# In this case, the X number of imputed tracks would be used directly for estimating behavioral states and/or habitat selection in place of the average location estimates from the SSM

# Conduct multiple imputation
mi_tracks <- sim_post(mp_fit, what = "predicted", reps = 50)

plot(mi_tracks[1,], type = "lines", alpha = 0.05, ortho = FALSE)
plot(mi_tracks[2,], type = "lines", alpha = 0.05)
plot(mi_tracks[3,], type = "lines", alpha = 0.05)
plot(mi_tracks[4,], type = "lines", alpha = 0.05)


# Convert to data.frame
mi_tracks2 <- mi_tracks |> 
  unnest(cols = psims) |> 
  select(-c(lon, lat)) |>  #remove incorrect coords
  add_trans_coords(coords = c('x','y'),
                   proj = "+proj=utm +zone=36 +ellps=WGS84 +units=km +no_defs +south",
                   new_proj = 4326) |> 
  mutate(id_orig = as.vector(str_match(id, "[0-9]+")),
         .after = id)


# Create custom viz
ggplot() +
  #plot imputed tracks (rep #s 1-50)
  geom_path(data = mi_tracks2 |> 
              filter(id_orig == 6469, rep > 0), aes(lon, lat, group = interaction(id,rep)), color = "dodgerblue",
            alpha = 0.2, linewidth = 0.1) +
  #plot predicted track (rep = 0)
  geom_path(data = mi_tracks2 |> 
              filter(id_orig == 6469, rep == 0), aes(lon, lat, group = id), color = "firebrick",
            linewidth = 0.5) +
  theme_bw() +
  coord_equal()






#######################################
### Export fitted models and tracks ###
#######################################

save(rw_fit, mp_fit, file = "processed_data/Session_3/ssm_fits.RData")  #model fit objects
write_csv(ssm_res, file = "processed_data/Session_3/regularized_tracks.csv")
write_csv(mi_tracks2, file = "processed_data/Session_3/mi_tracks.csv")
write_csv(dat2, file = "processed_data/Session_3/track_bursts.csv")
