
### Calculate autocorrelated kernel density estimates (AKDE) ###

library(tidyverse)
library(ctmm)
library(sf)
library(rnaturalearth)
library(tictoc)
library(terra)
library(progressr)
library(furrr)

source("R/utils.R")  #load in custom functions


###################
#### Load data ####
###################

dat <- read_csv('processed_data/cleaned_tracks.csv')

glimpse(dat)
summary(dat)




###############################
#### Wrangle and prep data ####
###############################

# Switch to Movebank naming conventions
dat2 <- dat |> 
  rename(tag.local.identifier = id,
         timestamp = date,
         location.long = lon,
         location.lat = lat) |> 
  mutate(GPS.Horizontal.Error = 30)  #add column to store presumed location error (30 m)

# Convert to 'telemetry' class
dat.telem <- as.telemetry(dat2, projection = "EPSG:32736")

dat.telem
summary(dat.telem)
plot(dat.telem, col = rainbow(length(dat.telem)))
#in this case, all IDs seem to be range-resident (i.e., not migrating/dispersing)



### Explore time intervals

# For one ID
dt.plot(dat.telem$`6471`)
#mostly at 1 hr 

# For all tracks
dt.plot(dat.telem)
#mostly at 1 hr, but also 2 and 4 hrs


# 1, 2, 4 hour sampling intervals
dt <- c(1,2,4) %#% "hour"



#########################
#### Plot variograms ####
#########################

# Generate variogram per track (using primary observed sampling intervals)
svf <- map(dat.telem, variogram, dt = dt)
level <- c(0.5, 0.95) # 50% and 95% CIs
xlim <- c(0,12) %#% "hour" # 0-12 hour window

### Explore individually
plot(svf$`5605`, xlim = xlim, level = level)
title("zoomed in")

plot(svf$`5605`, fraction = 0.65, level = level)
title("zoomed out")



### Explore together
par(mfrow = c(1,3))

# Zoomed-in
map(svf, \(x) plot(x, xlim = xlim, level = level))

# Zoomed-out
map(svf, \(x) plot(x, fraction = 0.65, level = level))

par(mfrow = c(1,1))
#all variograms look good
#if dispersal/migratory movements are present, tracks would need to be segmented to only analyze the range-resident components


### Viz population-level variogram (if similar across tracks)
svf_mean <- mean(svf)
plot(svf_mean, fraction = 0.65, level = level)






##################################################
### Fit continuous-time movement models (CTMM) ###
##################################################

# Guess at the variogram and model to be used per track
ctmm_guess <- map2(.x = dat.telem,
                   .y = svf,
                   .f = ~ctmm.guess(data = .x, variogram = .y, interactive = FALSE)
                   )
map(ctmm_guess, plot)


### Fit multiple CTMMs and select the best-fitting model

# Example for 1 track/burst
tic()
ctmm_5605 <- ctmm.select(data = dat.telem$`5605`, CTMM = ctmm_guess$`5605`, IC = 'AICc',
                         verbose = TRUE, trace = 1, cores = 10, method = 'pHREML')
toc()  #took 2.5 min

summary(ctmm_5605)
#shows model comparison by AICc (when `verbose = TRUE`), otherwise returns best model ONLY




# Example for multiple tracks/bursts
plan(multisession, workers = 3)  #run segments in parallel
handlers(handler_progress(format = ":spin :current/:total [:bar] :percent in :elapsed",
                                                width = 100,
                                                clear = FALSE))


with_progress({
  
  p <- progressor(along = dat.telem)  #create progress bar
  
  ctmm_fit <- future_map2(.x = dat.telem,
                          .y = ctmm_guess,
                          .f = ~{
                            tmp <- ctmm.select(data = .x, CTMM = .y, verbose = TRUE,
                                               method = 'pHREML', cores = 4)
                            p()  #print progress
                            tmp  #print results
                          },
                          .options = furrr_options(seed = 2026)
  )
  
})

plan(sequential)  #return to single core
# takes 7 min to run w/ 3 cores

map(ctmm_fit, summary)
#OUF anisotropic model best for all 3 tracks (of the 4 models compared)
#DOF (effective degrees of freedom) generally reflects number of range crossings; this should ideally be >3 for use in AKDE
#To compare only 'dispersive' models, specify `CTMM = ctmm(range = FALSE)` within ctmm::ctmm.guess()
#Can use ctmm.fit() to simply fit an IID model (essentially the same as KDE w/ href bandwidth estimator)




### Test fit of IID model
m.iid <- map(list(600, 60, 1200),  #values derived from asymptotes of variograms
             ~ctmm(sigma = .x %#% "km^2"))
ctmm_iid <- map2(.x = dat.telem,
                 .y = m.iid,
                ~ctmm.fit(data = .x, CTMM = .y))

# Append IID results per track
ctmm_all <- map2(.x = ctmm_fit,
                 .y = ctmm_iid,
                 .f = ~append(.x, list(IID = .y)))

# Re-do model comparison w/ IID included
map(ctmm_all, summary)
#IID models perform much worse (i.e., autocorrelation is informative!)


# Pull the best fitting models per track (listed first)
ctmm_fit_best <- map(ctmm_fit, pluck, 1)  #"plucks" first list element per track




#######################################
### Generate and viz AKDE per track ###
#######################################

### Fit standard AKDE model for 1 track
tic()
akde_5605 <- akde(data = dat.telem$`5605`, CTMM = ctmm_fit_best$`5605`,
             grid = list(dr = 1000, align.to.origin = TRUE, dr.fn = max))
toc()  #took 16 sec
#option 'dr' for arg `grid` specifies output spatial resolution (in meters)

summary(akde_5605)
plot(dat.telem$`5605`, UD = akde_5605)



### Fit weighted AKDE model for 1 track (to account for sampling bias from irregular time series)
dt1 <- 0.5 %#% "hour"  #setting to half of median time step

tic()
akde_5605w <- akde(data = dat.telem$`5605`, CTMM = ctmm_fit_best$`5605`, weights = TRUE, dt = dt1,
                   grid = list(dr = 1000, align.to.origin = TRUE, dr.fn = max))
toc()  #took 3.5 min

summary(akde_5605w)
plot(dat.telem$`5605`, UD = akde_5605w)
#very similar results



### Fit weighted AKDE across all tracks
tic()
akde <- akde(data = dat.telem, CTMM = ctmm_fit_best, weights = TRUE, dt = dt1,
             grid = list(dr = 1000, align.to.origin = TRUE, dr.fn = max))
toc()  #took 4 min

map(akde, summary)  #all DOFs > 4
plot(dat.telem, UD = akde, col = rainbow(length(dat.telem)))




### Extract 50 and 95% contours (core area, home range)
akde_sf <- map(akde,
               ~as.sf(.x, level.UD = c(0.5, 0.95))) |> 
  bind_rows() |> 
  filter(str_detect(name, "est")) |>  #keep only the mean prediction
  separate_wider_delim(cols = name, delim = " ", names = c("id", "level", NA)) |>  #split messy 'name' column
  st_sf(crs = 'epsg:32736')




### Viz map of estimates
africa <- ne_countries(scale = 10, continent = c("Africa"), returnclass = "sf") |>
  st_transform(crs = 32736)

# Plot separately by ID
ggplot() +
  geom_sf(data = africa |> 
            dplyr::select(-level)) +
  geom_sf(data = akde_sf, aes(color = level), fill = NA, linewidth = 0.5) +
  scale_color_brewer(palette = 'Set1') +
  theme_bw() +
  coord_sf(xlim = st_bbox(akde_sf)[c(1,3)],
           ylim = st_bbox(akde_sf)[c(2,4)]) +
  facet_grid(~id)


# Plot all tracks per UD level
ggplot() +
  geom_sf(data = africa |> 
            dplyr::select(-level)) +
  geom_path(data = dat |> 
              add_trans_coords(coords = c('lon','lat'), proj = 4326, new_proj = 32736),
            aes(x, y, color = factor(id), group = id), linewidth = 0.25, alpha = 0.5) +
  geom_sf(data = akde_sf, aes(color = factor(id)), fill = NA, linewidth = 1) +
  scale_color_brewer("ID", palette = "Dark2") +
  theme_bw() +
  coord_sf(xlim = st_bbox(akde_sf)[c(1,3)],
           ylim = st_bbox(akde_sf)[c(2,4)]) +
  facet_wrap(~ level, ncol = 2)


# Plotted all tracks and UD levels together
ggplot() +
  geom_sf(data = africa |> 
            st_transform(4326)) +
  geom_path(data = dat, aes(lon, lat, group = id, color = factor(id)), alpha = 0.5, linewidth = 0.25) +
  # Plot 95% isopleths
  geom_sf(data = akde_sf |> 
            st_transform(crs = 4326) |> 
            filter(level == "95%"),
          aes(color = factor(id)), fill = NA, linewidth = 0.75) +
  # Plot 50% isopleths
  geom_sf(data = akde_sf |> 
            st_transform(crs = 4326) |> 
            filter(level == "50%"),
          aes(color = factor(id), fill = factor(id)), linewidth = 0.75, alpha = 0.7) +
  scale_color_brewer("ID", palette = "Dark2") +
  scale_fill_brewer("ID", palette = "Dark2") +
  theme_bw() +
  coord_sf(xlim = range(dat$lon),
           ylim = range(dat$lat))







###########################################
### Measure space use overlap among IDs ###
###########################################

# Quantify Bhattacharyya's Affinity among UDs (w/ uncertainty)
overlap(akde)

# Show only the average estimated overlap
overlap(akde)$CI[,,"est"]






####################################################
### Assess pairwise interactions between animals ###
####################################################

### Calculate conditional distribution of encounters (CDE)
#a spatially defined PDF that describes the long-term encounter location probabilities for movement within home ranges

# Calculate for each pair
ctmm_cde <- combn(names(akde), 2, simplify = FALSE) |>  #create list of unique combos
  map(function(pair) {
    id1 <- pair[1]
    id2 <- pair[2]
    
    cde(akde[c(id1, id2)]) 
  }) |> 
  (\(x) set_names(x, map_chr(x, ~{.x@info$identity})))()  #set pair names per element

map(ctmm_cde, summary)

# Viz regions of potential encounters
plot(dat.telem[c('5605','6470')], UD = ctmm_cde$`5605 6470`, col = c('red','green'))
plot(dat.telem[c('5605','6471')], UD = ctmm_cde$`5605 6471`, col = c('red','blue'))
plot(dat.telem[c('6470','6471')], UD = ctmm_cde$`6470 6471`, col = c('green','blue'))




### Calculate pairwise distances over time

# Check if tracks co-occur in time
dat |> 
  summarize(.by = id,
            start = first(date),
            end = last(date))
#5605 and 6470 don't overlap in time

ctmm_dist <- combn(names(ctmm_fit_best), 2, simplify = FALSE) |>  #create list of unique combos
  discard_at(1) |> 
  map(function(pair) {
    id1 <- pair[1]
    id2 <- pair[2]
    
    distances(data = dat.telem[c(id1, id2)], CTMM = ctmm_fit_best[c(id1, id2)]) 
  }) |> 
  set_names(combn(names(ctmm_fit_best), 2, simplify = FALSE)[-1] |> 
              map(~{paste(.x[1],.x[2],sep = "_")}) |> 
              unlist())  #set pair names per element


ctmm_dist_df <- ctmm_dist |> 
  map(data.frame) |> 
  bind_rows(.id = "pair")


ggplot(ctmm_dist_df) +
  geom_path(aes(timestamp, est, color = pair)) +
  theme_bw() +
  labs(x = "Time", y = "Separation Distance (meters)")






########################################################
### Generate population-level inference of space use ###
########################################################

# mean() estimates population range; however, doesn't extrapolate to rest of pop.
#assumption that the sample of distributions is the population of distributions - such as for combining summer and winter ranges
pop_mean <- mean(akde)

summary(pop_mean)
plot(pop_mean)
plot(dat.telem, col = rainbow(length(dat.telem)), add = T)



# pkde() estimates range for entire pop. (via extrapolation) by accounting for inter-individual variability
#for when you have a small sample of a larger population - such as for a herd or colony
tic()
pkde <- pkde(data = dat.telem, UD = akde, kernel = "individual", ref = "Gaussian", weights = TRUE,
             dt = dt1, population = 3)
toc()  #took 16 min

summary(pkde)
plot(dat.telem, UD = pkde, col = rainbow(length(dat.telem)))
#in this case, results nearly identical to mean(); different results and shorter run-time when kernel = "population"


# meta() averages home range areas (in square meters) to get the population average and coefficient of variation
pop_meta <- ctmm::meta(akde, level.UD = 0.95, sort = TRUE)
#need to include 'ctmm' namespace explicitly due to conflict w/ another pkg

pop_meta






################################################################
### Make predictions of regularized tracks from fitted CTMMs ###
################################################################

### Predict from CTMMs
tic()
dat_pred <- map2(.x = ctmm_fit_best,
                 .y = dat.telem,
                 .f = ~ctmm::predict(object = .x, data = .y, dt = 3600, complete = TRUE))
toc()  # takes 8 s to run

# Convert from list to df
dat_pred_df <- dat_pred |> 
  map(data.frame) |> 
  bind_rows(.id = "id") |> 
  dplyr::select(id, timestamp, longitude, latitude, x, y, vx, vy)


# Viz regularized tracks
ggplot() +
  geom_point(data = dat_pred_df, aes(x, y, color = factor(id)), size = 0.1) +
  scale_color_brewer("ID", palette = "Dark2") +
  coord_equal() +
  theme_bw()




### Simulate from CTMMs
tic()
sim_5605 <- ctmm::simulate(object = ctmm_fit_best$`5605`, nsim = 5, seed = 2026, data = dat.telem$`5605`,
                           dt = 3600, complete = TRUE)
toc()  #took 7 sec

sim_5605_df <- sim_5605 |> 
  map(data.frame) |> 
  bind_rows(.id = "sim")


# Viz simulated tracks for 5605
ggplot() +
  geom_path(data = sim_5605_df, aes(x, y, color = factor(sim)), linewidth = 0.1, alpha = 0.5) +
  scale_color_brewer("Simulation #", palette = "Set1") +
  guides(color = guide_legend(override.aes = list(linewidth = 1, alpha = 1))) +
  coord_equal() +
  theme_bw()
#they're essentially all identical in this case





#####################################
### Perform workflow in Shiny app ###
#####################################

#remotes::install_github("ctmm-initiative/ctmmweb")
ctmmweb::app(dat.telem)





##########################################
#### Export datasets for easy loading ####
##########################################

save(akde_sf, file = "processed_data/AKDE_fits.RData")  #contours
save(ctmm_fit_best, file = "processed_data/CTMM_fits.RData")  #fitted CTMMs
save(pkde, file = "processed_data/PKDE_fits.RData")  #fitted PKDE
