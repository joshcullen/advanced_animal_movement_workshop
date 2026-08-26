
### Estimate behavioral states from HMM (discrete states) ###

library(momentuHMM)
library(tidyverse)
# library(bayesmove)
library(rnaturalearth)
library(sf)
library(tictoc)
library(furrr)

source("R/utils.R")

### Additional notes to discuss:
# * {momentuHMM} can estimate multiple discrete behavioral states using any number of covariates; this includes ancillary biologging data and evironmental variables
# * Options for BRW and CRW models, including activity centers that have attractive or repulsive forces
# * Option to perform multiple imputation via wrapper function for {crawl} methods, which can provide options to analyze tracks w/ measurement error and/or irregular time series indirectly
# * Can include random effects by ID on TPM
# * Can accommodate data streams/behaviors at multiple time scales (hierarchical HMM)



###################
#### Load data ####
###################

### Load tracks

# Load tracks w/ bursts
dat <- read_csv("processed_data/Session_3/track_bursts.csv")

# Load regularized tracks (mean est.)
ssm_tracks <- read_csv("processed_data/Session_3/regularized_tracks.csv") |> 
  arrange(id, date) |>  #make sure data is properly sorted
  mutate(id_orig = as.character(id_orig))  #IDs are better handled as 'character'

# Load multiple imputations from fitted SSM
mi_tracks <- read_csv(file = "processed_data/Session_3/mi_tracks.csv") |> 
  mutate(id_orig = as.character(id_orig))  #IDs are better handled as 'character'



### Load spatial layers
africa <- ne_countries(scale = 10, continent = c("Africa"), returnclass = "sf")
gl_pa <- st_read("raw_data/gltfca_protectedAreasDetailed.shp")





##############################
### Prep data for modeling ###
##############################

### If using irregular tracks and wanting to perform multiple imputation via {momentuHMM}/{crawl} ###

## Fit CTCRW via {crawl}

# format isotropic error ellipse
# dat2 <- cbind(dat,
#               crawl::argosDiag2Cov(Major = rep(30, nrow(dat)),  #specify GPS error in meters
#                                    Minor = rep(30, nrow(dat)),
#                                    Orientation = 0)
#               )

# Define how to handle error CTCRW model
# err.model <- list(x =  ~ ln.sd.x - 1,
#                   y =  ~ ln.sd.y - 1,
#                   rho =  ~ error.corr)

# Define prior for model param
ln.prior <- function(theta) dnorm(theta[2],-4,2,log=TRUE)


# Fit CTCRW
tic()
crw <- dat |> 
  rename(ID = id) |>  #needs to be capitalized for crawlWrap()
  add_trans_coords(coords = c('lon','lat'), proj = 4326, new_proj = 32736) |> 
  crawlWrap(timeStep = "1 hour",
            Time.name = "date",
            coord = c("x","y"),
            prior = ln.prior,
            theta = c(5,3))
toc()  #took 12 sec

# Calc SL and TA from fitted model
crw_dat <- prepData(crw)
plot(crw_dat)



### If using avg track locs from SSM regularization ###

# Change name of 'id' column
ssm_tracks2 <- ssm_tracks |>
  # add_trans_coords(coords = c('lon','lat'), proj = 4326, new_proj = 32736) |>  #preferred to have coords in m or km
  rename(ID = id) |>  #prepData() requires name as 'ID'; use bursts!
  data.frame()  #CAN'T be a 'tibble'

# Calc SL and TA (also creates as "momentuHMMData" class)
# to fit HMM, object needs to be of this class, making 'step' and 'angle' reserved colnames
ssm_tracks3 <- prepData(data = ssm_tracks2, type = 'UTM', coordNames = c('x','y'))  #can also be done using lat/long
plot(ssm_tracks3)



### If using multiple imputations generated w/ {aniMotum} ###

# Change name of 'id' column
mi_tracks2 <- mi_tracks |>
  filter(rep > 0) |>  #remove avg track fit
  mutate(ID = paste(id, rep, sep = "_")) |>  #create new ID that includes sim #
  data.frame()  #CAN'T be a 'tibble'

# Calc SL and TA (also creates as "momentuHMMData" class)
# to fit HMM, object needs to be of this class, making 'step' and 'angle' reserved colnames
mi_tracks3 <- prepData(data = mi_tracks2, type = 'UTM', coordNames = c('x','y'))  #can also be done using lat/long







#####################################################
### Filter out large temporal gaps (using bursts) ###
#####################################################

#-- Function moveHMM::splitAtGaps() can do this for you if only fitting HMM and wanting more concise approach --#

burst_windows <- dat |>
  group_by(id, burst_id) |>
  summarize(
    start_time = min(date),
    end_time = max(date),
    .groups = "drop"
  ) |> 
  ungroup() |> 
  mutate(ID = as.character(id))  #needs to match class for other df



# Filter tracks (i.e., remove interpolated section during long gaps)
crw_dat2 <- crw_dat |> 
  inner_join(burst_windows,
             by = join_by(ID, between(date, start_time, end_time))) |> 
  group_by(id, burst_id) |>
  filter(n() >= 6) |>  #remove bursts shorter than 6 hrs
  ungroup()

# ssm_tracks4 <- ssm_tracks3 |> 
#   inner_join(burst_windows,
#              by = join_by(ID, between(date, start_time, end_time))) |> 
#   group_by(id, burst_id) |>
#   filter(n() >= 6) |>  #remove bursts shorter than 6 hrs
#   ungroup()

mi_tracks4 <- mi_tracks3 |> 
  mutate(ID_sim = ID,
         ID = str_match(ID_sim, "[0-9]*")) |>  #need orig ID to join w/ 'burst_windows'
  inner_join(burst_windows,
             by = join_by(ID, between(date, start_time, end_time))) |> 
  group_by(id, burst_id) |>
  filter(n() >= 6) |>  #remove bursts shorter than 6 hrs
  ungroup()


# Example if using moveHMM::splitAtGaps
# tmp <- dat |> 
#   rename(ID = id, time = date) |>  #need to follow naming convention
#   moveHMM::splitAtGaps(maxGap = 8, shortestTrack = 6, units = "hours")





##########################################
### Define initial values for model(s) ###
##########################################

### Plot distributions of each movement metric (step length and turning angle)
ggplot(ssm_tracks3) +
  geom_density(aes(step), fill = "cadetblue") +
  labs(x = "Step Length (km)") +
  theme_bw()
#primarily < 2.5 km

ggplot(ssm_tracks3) +
  geom_density(aes(angle), fill = "firebrick") +
  labs(x = "Turning Angle (rad)") +
  theme_bw()
#primarily going straight (0 radians)




### Plot time series of metrics
ggplot(ssm_tracks3, aes(date, step)) +
  geom_line() +
  theme_bw(base_size = 14) +
  labs(x = "Date", y = "Step Length (km)") +
  facet_wrap(~ID, scales = "free_x")

ggplot(ssm_tracks3, aes(date, angle)) +
  geom_line() +
  theme_bw(base_size = 14) +
  labs(x = "Date", y = "Turning Angle (rad)") +
  facet_wrap(~ID, scales = "free_x", ncol = 1)
#Look for natural breaks that may define a particular state; this HIGHLY depends on the number of states (K) being estimated
#Time series of TA usually not very informative






### Pre-define states to set "good" initial values ###
## 2-state model ##

ssm_tracks4 <- ssm_tracks3 |>
  group_by(ID) |> 
  mutate(disp = sqrt((x - x[1])^2 + (y - y[1])^2)) |>  #calc net displacement
  ungroup() |> 
  data.frame()  #make sure to use data.frame for model fitting

tmp <- ssm_tracks4 |> 
  mutate(phase = case_when(step > 0.75 ~ 'Transit',
                           TRUE ~ 'ARS'))

ggplot(tmp, aes(date, step)) +
  geom_path(aes(group = ID, color = phase)) +
  theme_bw() +
  facet_wrap(~ID, scales = "free_x")

ggplot(tmp, aes(date, disp)) +
  geom_path(aes(group = ID, color = phase)) +
  theme_bw() +
  facet_wrap(~ID, scales = "free")
#looks like it does a decent job


# Check for obs w/ SL = 0
sum(tmp$step == 0, na.rm = TRUE)  #none




# Get summary stats
tmp |>
  summarize(.by = phase,
            mean.step = mean(step, na.rm = T),
            sd.step = sd(step, na.rm = T))
#ARS: mean = 0.25; SD = 0.2
#Transit: mean = 1.5; SD = 0.75
#we'll assume both states centered at 0 w/ differing concentrations





#######################
### Fit 2-state HMM ###
#######################

# initial step distribution natural scale parameters
stepPar0 <- c(0.25, 1.5, 0.2, 0.75) # (mu_1, mu_2, sd_1, sd_2)

# initial angle distribution natural scale parameters
# anglePar0 <- c(0, 0, 0.5, 0.9) # (mean_1, mean_2, concentration_1, concentration_2)
anglePar0 <- c(0.5, 0.9) # (concentration_1, concentration_2)

# Convert back to "momentuHMMData" class to fit model
class(ssm_tracks4) <- append("momentuHMMData", class(ssm_tracks4))
ssm_tracks4$ID <- as.factor(ssm_tracks4$ID)  #needs to be factor for TMB optimization


# Fit model
set.seed(2026)
#if states get flipped, try adjusting initial params or random seed number (e.g., seed 2026 gave problems w/ flipped states when using retryFits)
#alternatively, you can specify pseudo-design matrix to prevent state label switching (but doesn't work for TMB)
# K <- 2
# angleDM <- diag(K)  #pseudo-design matrix: KxK identity matrix keeping all states independent for angle conc.
# rownames(angleDM) <- paste0("concentration_", 1:K)
# colnames(angleDM) <- paste0("concentration_", 1:K, ":(Intercept)")

tic()
test_2states <- fitHMM(data = ssm_tracks4,
                       nbStates = 2,
                       dist = list(step = "gamma", angle = "wrpcauchy"),  #can use other distribs as well
                       Par0 = list(step = stepPar0, angle = anglePar0),
                       formula = ~ 1,
                       # DM = list(angle = angleDM),
                       stationary = TRUE,
                       estAngleMean = list(angle=FALSE),
                       stateNames = c('ARS', 'Transit'),
                       optMethod = "TMB",
                       # ncores = 10,
                       # retryFits = 10
)
toc()  #took 4 sec to run; 75 sec w/ retryFits via 'nlm' optim
#If problems w/ model fit or convergence, you can try 1) changing the random seed value, 2) changing the selected distributions for the movement metrics, 3) redefine new initial values for params, 4) increase number of retryFits, 5) adjust retrySD, 6) change optimization method, 7) set 'estAngleMean' to FALSE, 8) set `stationary = FALSE`, 9) specify DM, userBounds, and/or workBounds


test_2states
plot(test_2states)
plotPR(test_2states, ncores = 5)  #plot of pseudo-residuals looks pretty good; some autocorr in 'step' every 24 h




# Fit w/ random perturbations to ensure global ML found
#given the way that values are tweaked for built-in 'retryFits' arg, this method may work better sometimes
set.seed(2026)
tic()
fit_hmm_2states <- run_HMMs(data = ssm_tracks4,
                            K = 2,  #number of states
                            dist = list(step = "gamma", angle = "wrpcauchy"),
                            Par0 = list(step = stepPar0, angle = anglePar0),
                            stationary = TRUE,
                            estAngleMean = list(angle=FALSE),
                            state.names = c('ARS','Transit'),
                            optMethod = "TMB",
                            niter = 30,
                            ncores = 10)
toc()  #took 20 sec to run

fit_hmm_2states
fit_hmm_2states$mod$code  #check that it converged (code = 0)
plot(fit_hmm_2states)
plotPR(fit_hmm_2states, ncores = 5)







#######################
### Fit 3-state HMM ###
#######################

# Plot time series of SL
ggplot(ssm_tracks4, aes(date, step)) +
  geom_line() +
  theme_bw(base_size = 14) +
  labs(x = "Date", y = "Step Length (km)") +
  facet_wrap(~ID, scales = "free_x")


# Manually classify to determine initial values
tmp <- ssm_tracks4 |>
  group_by(ID) |> 
  mutate(phase = case_when(step >= 1.5 ~ 'Transit',
                           step >= 0.5 & step < 1.5 ~ 'Exploratory',
                           TRUE ~ 'Encamped')) |> 
  ungroup() |> 
  data.frame()  #make sure to use data.frame for model fitting

ggplot(tmp, aes(date, step)) +
  geom_path(aes(group = ID, color = phase)) +
  theme_bw() +
  facet_wrap(~ID, scales = "free_x")

ggplot(tmp, aes(date, disp)) +
  geom_path(aes(group = ID, color = phase)) +
  theme_bw() +
  facet_wrap(~ID, scales = "free")
#looks like it does a decent job


# Get summary stats
tmp |>
  summarize(.by = phase,
            mean.step = mean(step, na.rm = T),
            sd.step = sd(step, na.rm = T))
#Encamped: mean = 0.2; SD = 0.15
#Exploratory: mean = 0.85; SD = 0.25
#Transit: mean = 2.25; SD = 0.75




### Define inits

# initial step distribution natural scale parameters
stepPar0 <- c(0.2, 0.85, 2.25, 0.15, 0.25, 0.75) # (mu_1, mu_2, mu_3, sd_1, sd_2, sd_3)

# initial angle distribution natural scale parameters
anglePar0 <- c(0.5, 0.75, 0.9) # (conc_1, conc_2, conc_3); assuming mean fixed at 0

# Convert back to "momentuHMMData" class to fit model
# class(ssm_tracks5) <- append("momentuHMMData", class(ssm_tracks5))
# ssm_tracks5$ID <- as.factor(ssm_tracks5$ID)  #needs to be factor for TMB optimization



# Fit single model
set.seed(2026)
tic()
test_3states <- fitHMM(data = ssm_tracks4,
                       nbStates = 3,
                       dist = list(step = "gamma", angle = "wrpcauchy"),  #can use other distribs as well
                       Par0 = list(step = stepPar0, angle = anglePar0),
                       formula = ~ 1,
                       stationary = TRUE,
                       estAngleMean = list(angle=FALSE),  
                       stateNames = c('Encamped','Exploratory','Transit'),
                       optMethod = "TMB"
)
toc()  #took 7 sec to run

test_3states
plot(test_3states)
plotPR(test_3states, ncores = 5)  #looks pretty good






# Fit HMM (w/ random perturbations)
set.seed(2026)
tic()
fit_hmm_3states <- run_HMMs(data = ssm_tracks4,
                            K = 3,  #number of states
                            dist = list(step = "gamma", angle = "wrpcauchy"),
                            Par0 = list(step = stepPar0, angle = anglePar0),
                            stationary = TRUE,
                            estAngleMean = list(angle=FALSE),
                            state.names = c('Encamped','Exploratory','Transit'),
                            optMethod = "TMB",
                            niter = 30,
                            ncores = 10)
toc()  #took 30 sec to run

fit_hmm_3states
fit_hmm_3states$mod$code  #check that it converged (code = 0)
plot(fit_hmm_3states)
plotPR(fit_hmm_3states, ncores = 5)







#######################
### Fit 4-state HMM ###
#######################

# Plot time series of SL
ggplot(ssm_tracks4, aes(date, step)) +
  geom_line() +
  theme_bw(base_size = 14) +
  labs(x = "Date", y = "Step Length (km)") +
  facet_wrap(~ID, scales = "free_x")


# Manually classify to determine initial values
tmp <- ssm_tracks4 |>
  group_by(ID) |> 
  mutate(phase = case_when(step >= 3 ~ 'Ranging',
                           step >= 1.5 & step < 3 ~ 'Transit',
                           step >= 0.5 & step < 1.5 ~ 'Exploratory',
                           step < 0.5 ~ 'Encamped'),
         phase = factor(phase, levels = c('Encamped','Exploratory','Transit','Ranging'))) |> 
  ungroup() |> 
  data.frame()  #make sure to use data.frame for model fitting

ggplot(tmp, aes(date, step)) +
  geom_path(aes(group = ID, color = phase)) +
  theme_bw() +
  facet_wrap(~ID, scales = "free_x")

ggplot(tmp, aes(date, disp)) +
  geom_path(aes(group = ID, color = phase)) +
  theme_bw() +
  facet_wrap(~ID, scales = "free")
#looks like it does a decent job


# Get summary stats
tmp |>
  summarize(.by = phase,
            mean.step = mean(step, na.rm = T),
            sd.step = sd(step, na.rm = T))
#Encamped: mean = 0.2; SD = 0.15
#Exploratory: mean = 0.85; SD = 0.25
#Transit: mean = 2; SD = 0.4
#Ranging: mean = 3.5; SD = 0.7





### Define inits

# initial step distribution natural scale parameters
stepPar0 <- c(0.2, 0.85, 2, 3.5,  #(mu_1, mu_2, mu_3, mu_4)
              0.15, 0.25, 0.4, 0.7)  #(sd_1, sd_2, sd_3, sd_4)

# initial angle distribution natural scale parameters
anglePar0 <- c(0.4, 0.65, 0.8, 0.9)  #(conc_1, conc_2, conc_3, conc_4)

# Convert back to "momentuHMMData" class to fit model
# class(ssm_tracks5) <- append("momentuHMMData", class(ssm_tracks5))
# ssm_tracks5$ID <- as.factor(ssm_tracks5$ID)  #needs to be factor for TMB optimization



# Fit single model
set.seed(2026)
tic()
test_4states <- fitHMM(data = ssm_tracks4,
                       nbStates = 4,
                       dist = list(step = "gamma", angle = "wrpcauchy"),  #can use other distribs as well
                       Par0 = list(step = stepPar0, angle = anglePar0),
                       formula = ~ 1,
                       stationary = FALSE,  #changed to try to get model to converge
                       estAngleMean = list(angle=FALSE),  
                       stateNames = c('Encamped','Exploratory','Transit','Ranging'),
                       optMethod = "TMB"
)
toc()  #took 12 sec to run

test_4states
plot(test_4states)
plotPR(test_4states, ncores = 5)






# Fit HMM (w/ random perturbations)
set.seed(2026)
tic()
fit_hmm_4states <- run_HMMs(data = ssm_tracks4,
                            K = 4,  #number of states
                            dist = list(step = "gamma", angle = "wrpcauchy"),
                            Par0 = list(step = stepPar0, angle = anglePar0),
                            stationary = TRUE,
                            estAngleMean = list(angle=FALSE),
                            state.names = c('Encamped','Exploratory','Transit','Ranging'),
                            optMethod = "TMB",
                            niter = 30,
                            ncores = 10)
toc()  #took 45 sec to run

fit_hmm_4states
fit_hmm_4states$mod$code  #didn't converge
plot(fit_hmm_4states)
plotPR(fit_hmm_4states, ncores = 5)





###############################
### Perform order selection ###
###############################

#we're now going to select the "best" model of the 3 fitted
#using a combination of factors (laid out in more detail in Pohle et al., 2017)
#general recommendation is that you typically can only estimate a max of M + 1 state, where 'M' represents the number of data streams analyzed; so for only 2 data streams (SL, TA), 3 is probably our max

# AIC comparison
AIC(fit_hmm_2states, fit_hmm_3states, fit_hmm_4states) |> 
  mutate(dAIC = AIC - min(AIC))
#4-state model "best" per AIC; but IC often overfits, so take w/ grain of salt

# BIC comparison
BIC(fit_hmm_2states, fit_hmm_3states, fit_hmm_4states) |> 
  mutate(dBIC = BIC - min(BIC))
#4-state model "best" per BIC; but IC often overfits, so take w/ grain of salt


# Compare state-dependent density distributions and annotated tracks
plot(fit_hmm_2states)
plot(fit_hmm_3states)
plot(fit_hmm_4states)
#all look reasonable, but 4-state model probs overfits and model didn't converge
#so we'll select the 3-state model here as the best-fitting of the three





#############################################
### Explore 3-state model more thoroughly ###
#############################################

# Explore time series of state probabilities
plotStates(fit_hmm_3states)

# Annotate tracks w/ model results
ssm_tracks5 <- ssm_tracks4 |> 
  mutate(state_vit = viterbi(fit_hmm_3states),  #use Viterbi algorithm (most likely sequence)
         state_vit = case_when(state_vit == 1 ~ 'Encamped',
                               state_vit == 2 ~ 'Exploratory',
                               state_vit == 3 ~ 'Transit')) |> 
  cbind(stateProbs(fit_hmm_3states)) |>  #uses forward-backward algorithm (probabilities of each state)
  mutate(state_fb = case_when(  #using est. from forward-backward algorithm, assign states when likely (Pr > 50%)
    Encamped > 0.5 ~ "Encamped",
    Exploratory > 0.5 ~ "Exploratory",
    Transit > 0.5 ~ "Transit",
    TRUE ~ 'Unclassified'  # Assigns NA if no value is > 0.5
  ))

# Any discrepencies between the 2 methods for assigning states?
all.equal(ssm_tracks5$state_vit, ssm_tracks5$state_fb)
#looks like there are some differences (1833 to be exact)

# Any 'Unclassified' obs?
table(ssm_tracks5$state_fb)  #5 obs


#-- In general, I like to use the approach of assigning states based on the confidence in the state estimates. But plenty of people use and publish the Viterbi seq. of states too --#




############################
### Viz annotated tracks ###
############################

# Plot all annotated tracks together
ggplot() +
  geom_sf(data = africa) +
  geom_path(data = ssm_tracks5, aes(lon, lat, group = ID), alpha = 0.5, linewidth = 0.25) +
  geom_point(data = ssm_tracks5, aes(lon, lat, color = state_fb), alpha = 0.5, size = 1) +
  geom_sf(data = gl_pa, color = "black", fill = NA, linewidth = 0.25) +
  scale_color_manual("State", values = c(RColorBrewer::brewer.pal(n = 3, "Dark2"), "grey")) +
  theme_bw() +
  coord_sf(xlim = range(dat$lon),
           ylim = range(dat$lat))


# Facet by ID
ggplot() +
  geom_sf(data = africa) +
  geom_path(data = ssm_tracks5, aes(lon, lat, group = ID), alpha = 0.5, linewidth = 0.25) +
  geom_point(data = ssm_tracks5, aes(lon, lat, color = state_fb), alpha = 0.5, size = 1) +
  geom_sf(data = gl_pa, color = "black", fill = NA, linewidth = 0.25) +
  scale_color_manual("State", values = c(RColorBrewer::brewer.pal(n = 3, "Dark2"), "grey")) +
  theme_bw() +
  coord_sf(xlim = range(dat$lon),
           ylim = range(dat$lat)) +
  facet_wrap(~id_orig, ncol = 2)


for (i in 1:n_distinct(ssm_tracks5$id_orig)) {
  print(
    ggplot() +
      geom_sf(data = africa) +
      geom_path(data = ssm_tracks5, aes(lon, lat, group = ID), alpha = 0.5, linewidth = 0.25) +
      geom_point(data = ssm_tracks5, aes(lon, lat, color = state_fb), alpha = 0.5, size = 1) +
      geom_sf(data = gl_pa, color = "black", fill = NA, linewidth = 0.25) +
      scale_color_manual("State", values = c(RColorBrewer::brewer.pal(n = 3, "Dark2"), "grey")) +
      theme_bw() +
      coord_sf(xlim = range(dat$lon),
               ylim = range(dat$lat)) +
      ggforce::facet_wrap_paginate(~id_orig, ncol = 1, nrow = 1, page = i)
  )
}


# Separate by state
for (i in 1:n_distinct(ssm_tracks5$id_orig)) {
  print(
    ggplot() +
      geom_sf(data = africa) +
      geom_path(data = ssm_tracks5 |> 
                  select(-state_fb), aes(lon, lat, group = ID), alpha = 0.5, linewidth = 0.25) +
      geom_point(data = ssm_tracks5, aes(lon, lat, color = state_fb), alpha = 0.5, size = 1) +
      geom_sf(data = gl_pa, color = "black", fill = NA, linewidth = 0.25) +
      scale_color_manual("State", values = c(RColorBrewer::brewer.pal(n = 3, "Dark2"), "grey")) +
      theme_bw() +
      coord_sf(xlim = range(dat$lon),
               ylim = range(dat$lat)) +
      ggforce::facet_wrap_paginate(~state_fb, ncol = 1, nrow = 1, page = i)
  )
}


# Explore state probabilities (e.g., 'Transit' in ID 6469)
ggplot() +
  geom_path(data = ssm_tracks5 |> 
              filter(id_orig == 6469), aes(lon, lat, group = ID), alpha = 0.5, linewidth = 0.25) +
  geom_point(data = ssm_tracks5 |> 
               filter(id_orig == 6469), aes(lon, lat, color = Transit), alpha = 0.5, size = 1) +
  scale_color_distiller("Pr(Transit)", palette = "Spectral", direction = -1, limits = c(0,1)) +
  theme_bw() +
  coord_equal()



# Interactively explore results
ssm_tracks5 |> 
  bayesmove::shiny_tracks(epsg = "+proj=utm +zone=36 +ellps=WGS84 +units=km +no_defs +south")




###############################
### Export annotated tracks ###
###############################

write_csv(ssm_tracks5, "processed_data/Session_3/HMM_3state_simple.csv")
