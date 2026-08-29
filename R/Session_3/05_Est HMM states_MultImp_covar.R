
### Estimate behavioral states from HMM (discrete states) ###
## Multiple Imputation with Covariates ##

library(momentuHMM)
library(tidyverse)
library(rnaturalearth)
library(sf)
library(tictoc)
library(furrr)

source("R/utils.R")



###################
#### Load data ####
###################

### Load tracks

# Load tracks w/ bursts
dat <- read_csv("processed_data/Session_3/track_bursts.csv")



### Load spatial layers
africa <- ne_countries(scale = 10, continent = c("Africa"), returnclass = "sf")
gl_pa <- st_read("raw_data/gltfca_protectedAreasDetailed.shp")





##############################
### Prep data for modeling ###
##############################

### If using irregular tracks and wanting to perform multiple imputation via {momentuHMM}/{crawl} ###

## Fit CTCRW via {crawl}

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
            theta = c(6.855, -0.007),
            ncores = 3
  )
toc()  #took 15 sec


# Viz what these CTMM predictions look like
ggplot() +
  # Plot lines using observed locs
  geom_path(data = crw$crwPredict |> 
              data.frame() |> 
              filter(locType == 'o'), aes(x, y, color = factor(ID), group = ID), linewidth = 0.1) +
  # Plot points using predicted locs
  geom_point(data = crw$crwPredict |> 
               data.frame() |> 
               filter(locType == 'p'), aes(mu.x, mu.y, color = factor(ID)), size = 0.2, alpha = 0.5) +
  theme_bw() +
  coord_equal()



# Calc SL and TA from fitted model
crw_dat <- prepData(crw)
plot(crw_dat)

head(crw_dat)  #we now have hourly (regularized) data





##########################################
### Interpolate values for temperature ###
##########################################

# Perform linear interpolation to add temperature to predicted locs
cov_interp <- crw$crwPredict |>
  group_by(ID) |>
  mutate(
    # Interpolate temperature based on the numeric time vector
    temperature = approx(
      x = date, 
      y = temp, 
      xout = date, 
      rule = 2 # Prevents NAs if a 'p' step occurs just before the first 'o' step
    )$y,
    # Extract hour
    hour = hour(date)
  ) |>
  ungroup() 

# Visually explore that these interpolations make sense
ggplot() +
  geom_point(data = cov_interp |> 
               group_by(ID) |> 
               slice_head(n = 500) |> 
               ungroup() |> 
               filter(locType == 'p'), aes(date, temperature), color = "red") +  #predicted
  geom_point(data = cov_interp |> 
               group_by(ID) |> 
               slice_head(n = 500) |> 
               ungroup() |> 
               filter(locType == 'o'), aes(date, temperature), color = "black", alpha = 0.6) +  #observed
  theme_bw() +
  facet_wrap(~ID, scales = "free_x", ncol = 1)
#looks like it does a good job

ggplot() +
  geom_point(data = cov_interp |> 
               group_by(ID) |> 
               slice_head(n = 500) |> 
               ungroup() |> 
               filter(locType == 'p'), aes(date, hour), color = "red") +  #predicted
  geom_point(data = cov_interp |> 
               group_by(ID) |> 
               slice_head(n = 500) |> 
               ungroup() |> 
               filter(locType == 'o'), aes(date, hour), color = "black", alpha = 0.6) +  #observed
  theme_bw() +
  facet_wrap(~ID, scales = "free_x", ncol = 1)
#looks like it does a good job


# Subset and convert to df
cov_interp2 <- cov_interp|> 
  data.frame() |>  #can't be a tibble
  select(ID, date, temperature, hour)  #only keep necessary cols

# Merge the interpolated covariates into the 'crw' object
crw_merged <- crawlMerge(crw, cov_interp2, Time.name = "date")





##########################################
### Define initial values for model(s) ###
##########################################

### Plot distributions of each movement metric (step length and turning angle)
ggplot(crw_dat) +
  geom_density(aes(step), fill = "cadetblue") +
  labs(x = "Step Length (m)") +
  theme_bw()
#primarily < 2500 km

ggplot(crw_dat) +
  geom_density(aes(angle), fill = "firebrick") +
  labs(x = "Turning Angle (rad)") +
  theme_bw()
#primarily going straight (0 radians)




### Plot time series of metrics
ggplot(crw_dat, aes(date, step)) +
  geom_line() +
  theme_bw(base_size = 14) +
  labs(x = "Date", y = "Step Length (m)") +
  facet_wrap(~ID, scales = "free_x")

ggplot(crw_dat, aes(date, angle)) +
  geom_line() +
  theme_bw(base_size = 14) +
  labs(x = "Date", y = "Turning Angle (rad)") +
  facet_wrap(~ID, scales = "free_x", ncol = 1)
#Look for natural breaks that may define a particular state; this HIGHLY depends on the number of states (K) being estimated
#Time series of TA usually not very informative






### Pre-define states to set "good" initial values ###
## 2-state model ##

crw_dat2 <- crw_dat |>
  group_by(ID) |> 
  mutate(disp = sqrt((x - x[1])^2 + (y - y[1])^2)) |>  #calc net displacement
  ungroup() |> 
  data.frame()  #make sure to use data.frame for model fitting

tmp <- crw_dat2 |> 
  mutate(phase = case_when(step > 750 ~ 'Transit',
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
#ARS: mean = 250; SD = 200
#Transit: mean = 1500; SD = 750
#we'll assume both states centered at 0 w/ differing concentrations





#######################
### Fit 2-state HMM ###
#######################

# initial step distribution natural scale parameters
stepPar0 <- c(250, 1500, 200, 750) # (mu_1, mu_2, sd_1, sd_2)

# initial angle distribution natural scale parameters
anglePar0 <- c(0.5, 0.9) # (concentration_1, concentration_2)



# Fit model w/o covars first
set.seed(2026)
tic()
test_2states <- MIfitHMM(miData = crw_merged,
                         nSims = 3,  #start w/ 3 imputations as a test
                         nbStates = 2,
                         ncores = 3,
                         dist = list(step = "gamma", angle = "wrpcauchy"),  #can use other distribs as well
                         Par0 = list(step = stepPar0, angle = anglePar0),
                         formula = ~ 1,
                         covNames = c("temperature", "hour"),  # THIS PREVENTS STRIPPING FROM OBJECT
                         stationary = TRUE,
                         estAngleMean = list(angle=FALSE),
                         stateNames = c('ARS', 'Transit'),
                         optMethod = "TMB"
)
toc()  #took 1 min to run for 3 sims

test_2states



## Now re-fit w/ covars using fitted model params

# Define formula on tpm
form <- ~ temperature * cosinor(hour, period = 24)

# initial parameters (obtained from 'test_2states')
Par0_2states <- getPar0(model = test_2states, formula = form)


# Fit model w/ covars
set.seed(2026)
tic()
fit_hmm_2states <- MIfitHMM(miData = crw_merged,
                            nSims = 10,  #10 imputations
                            nbStates = 2,
                            ncores = 10,
                            dist = list(step = "gamma", angle = "wrpcauchy"),
                            
                            #Using initial vals here from simpler model for step, angle, and betas
                            Par0 = list(step = Par0_2states$Par$step, angle = Par0_2states$Par$angle),
                            beta0 = Par0_2states$beta,
                            
                            formula = ~ temperature * cosinor(hour, period = 24),
                            covNames = c("temperature", "hour"),  # THIS PREVENTS STRIPPING FROM OBJECT
                            stationary = FALSE,  #needs to be FALSE when using covars
                            estAngleMean = list(angle=FALSE),
                            stateNames = c('ARS', 'Transit'),
                            optMethod = "TMB"
)
toc()  #took 2.5 min to run for 10 sims

fit_hmm_2states  #most imputations didn't converge
plot(fit_hmm_2states, plotCI = TRUE)
plot(fit_hmm_2states, plotCI = TRUE, plotTracks = FALSE,
     covs = data.frame(hour = 16, temperature = 45))  #manually set fixed values

plotStationary(fit_hmm_2states, plotCI = TRUE)
plotStationary(fit_hmm_2states, plotCI = TRUE,
               covs = data.frame(hour = 16, temperature = 45))  #manually set fixed values


plotPR(fit_hmm_2states, ncores = 5)  #looks pretty good







#######################
### Fit 3-state HMM ###
#######################

# Plot time series of SL
ggplot(crw_dat2, aes(date, step)) +
  geom_line() +
  theme_bw(base_size = 14) +
  labs(x = "Date", y = "Step Length (m)") +
  facet_wrap(~ID, scales = "free_x")


# Manually classify to determine initial values
tmp <- crw_dat2 |>
  group_by(ID) |> 
  mutate(phase = case_when(step >= 1500 ~ 'Transit',
                           step >= 500 & step < 1500 ~ 'Exploratory',
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
#Encamped: mean = 200; SD = 150
#Exploratory: mean = 850; SD = 250
#Transit: mean = 2250; SD = 750




### Define inits

# initial step distribution natural scale parameters
stepPar0 <- c(250, 850, 2250, 150, 250, 750) # (mu_1, mu_2, mu_3, sd_1, sd_2, sd_3)

# initial angle distribution natural scale parameters
anglePar0 <- c(0.25, 0.5, 0.9) # (concentration_1, concentration_2, concentration_3)



### Fit model
set.seed(2026)
tic()
test_3states <- MIfitHMM(miData = crw_merged,
                         nSims = 3,  #start w/ 3 imputations as a test
                         nbStates = 3,
                         ncores = 3,
                         dist = list(step = "gamma", angle = "wrpcauchy"),  #can use other distribs as well
                         Par0 = list(step = stepPar0, angle = anglePar0),
                         formula = ~ 1,
                         covNames = c("temperature", "hour"),  # THIS PREVENTS STRIPPING FROM OBJECT
                         stationary = TRUE,
                         estAngleMean = list(angle=FALSE),
                         stateNames = c('Encamped','Exploratory','Transit'),
                         optMethod = "TMB"
)
toc()  #took 1 min to run for 3 sims

test_3states  #looks pretty good



## Now re-fit w/ covars using fitted model params

# initial parameters (obtained from 'test_3states')
Par0_3states <- getPar0(model = test_3states, formula = form)


# Fit model w/ covars
set.seed(2026)
tic()
fit_hmm_3states <- MIfitHMM(miData = crw_merged,
                            nSims = 10,  #10 imputations
                            nbStates = 3,
                            ncores = 10,
                            na.rm = TRUE,  #remove imputations w/ NA params or SEs from pooling
                            dist = list(step = "gamma", angle = "wrpcauchy"),  #can use other distribs as well
                            
                            #Using initial vals here from simpler model for step, angle, and betas
                            Par0 = list(step = Par0_3states$Par$step, angle = Par0_3states$Par$angle),
                            beta0 = Par0_3states$beta,
                            
                            formula = ~ temperature * cosinor(hour, period = 24),
                            covNames = c("temperature", "hour"),  # THIS PREVENTS STRIPPING FROM OBJECT
                            stationary = FALSE,  #needs to be FALSE when using covars
                            estAngleMean = list(angle=FALSE),
                            stateNames = c('Encamped','Exploratory','Transit'),
                            optMethod = "TMB"
)
toc()  #took 2 min to run for 10 sims

fit_hmm_3states  #most imputations didn't converge and/or have missing uncertainty est.
plot(fit_hmm_3states, plotCI = TRUE)
plot(fit_hmm_3states, plotCI = TRUE, plotTracks = FALSE,
     covs = data.frame(hour = 16, temperature = 45))  #manually set fixed values

plotStationary(fit_hmm_3states, plotCI = TRUE)
plotStationary(fit_hmm_3states, plotCI = TRUE,
               covs = data.frame(hour = 16, temperature = 45))  #manually set fixed values

plotPR(fit_hmm_3states, ncores = 5)  #looks fine







#######################
### Fit 4-state HMM ###
#######################

# Plot time series of SL
ggplot(crw_dat2, aes(date, step)) +
  geom_line() +
  theme_bw(base_size = 14) +
  labs(x = "Date", y = "Step Length (m)") +
  facet_wrap(~ID, scales = "free_x")


# Manually classify to determine initial values
tmp <- crw_dat2 |>
  group_by(ID) |> 
  mutate(phase = case_when(step >= 3000 ~ 'Ranging',
                           step >= 1500 & step < 3000 ~ 'Transit',
                           step >= 500 & step < 1500 ~ 'Exploratory',
                           step < 500 ~ 'Encamped'),
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
#Encamped: mean = 200; SD = 150
#Exploratory: mean = 850; SD = 250
#Transit: mean = 2000; SD = 400
#Ranging: mean = 3500; SD = 700




### Define inits

# initial step distribution natural scale parameters
stepPar0 <- c(200, 850, 2000, 3500, 150, 250, 400, 700)  #(mu_1, mu_2, mu_3, mu_4, sd_1, sd_2, sd_3, sd_4)

# initial angle distribution natural scale parameters
anglePar0 <- c(0.1, 0.5, 0.7, 0.9)  #(conc_1, conc_2, conc_3, conc_4)



### Fit model
set.seed(2026)
tic()
test_4states <- MIfitHMM(miData = crw_merged,
                         nSims = 3,  #start w/ 3 imputations as a test
                         nbStates = 4,
                         ncores = 3,
                         dist = list(step = "gamma", angle = "wrpcauchy"),  #can use other distribs as well
                         Par0 = list(step = stepPar0, angle = anglePar0),
                         formula = ~ 1,
                         covNames = c("temperature", "hour"),  # THIS PREVENTS STRIPPING FROM OBJECT
                         stationary = TRUE,
                         estAngleMean = list(angle=FALSE),
                         stateNames = c('Encamped','Exploratory','Transit','Ranging'),
                         optMethod = "TMB"
)
toc()  #took 1 min to run for 3 sims

test_4states  #looks pretty good



## Now re-fit w/ covars using fitted model params

# initial parameters (obtained from 'test_4states')
Par0_4states <- getPar0(model = test_4states, formula = form)


# Fit model w/ covars
set.seed(2026)
tic()
fit_hmm_4states <- MIfitHMM(miData = crw_merged,
                            nSims = 10,  #10 imputations
                            nbStates = 4,
                            ncores = 10,
                            na.rm = TRUE,  #remove imputations w/ NA params or SEs from pooling
                            dist = list(step = "gamma", angle = "wrpcauchy"),  #can use other distribs as well
                            
                            #Using initial vals here from simpler model for step, angle, and betas
                            # Par0 = list(step = Par0_4states$Par$step, angle = Par0_4states$Par$angle),
                            # beta0 = Par0_4states$beta,
                            
                            formula = ~ temperature * cosinor(hour, period = 24),
                            covNames = c("temperature", "hour"),  # THIS PREVENTS STRIPPING FROM OBJECT
                            stationary = FALSE,  #needs to be FALSE when using covars
                            estAngleMean = list(angle=FALSE),
                            stateNames = c('Encamped','Exploratory','Transit','Ranging'),
                            optMethod = "TMB"
)
toc()  #took 2 min to run for 10 sims

fit_hmm_4states  #no models converged
plot(fit_hmm_4states, plotCI = TRUE)
plot(fit_hmm_4states, plotCI = TRUE, plotTracks = FALSE,
     covs = data.frame(hour = 16, temperature = 45))  #manually set fixed values

plotStationary(fit_hmm_4states, plotCI = TRUE)
plotStationary(fit_hmm_4states, plotCI = TRUE,
               covs = data.frame(hour = 16, temperature = 45))  #manually set fixed values

plotPR(fit_hmm_4states, ncores = 5)  #looks fine





###############################
### Perform order selection ###
###############################

#we're now going to select the "best" model of the 3 fitted
#using a combination of factors (laid out in more detail in Pohle et al., 2017)
#general recommendation is that you typically can only estimate a max of M + 1 state, where 'M' represents the number of data streams analyzed; so for only 2 data streams (SL, TA), 3 is probably our max

# No clear way to calculate AIC/BIC since now we're dealing w/ multiple results per track. This isn't the best metric for order selection anyway, so we'll rely on the GOF diagnostics, state-dependent density distribs, and annotated tracks


# Compare state-dependent density distributions and annotated tracks
plot(fit_hmm_2states)
plot(fit_hmm_3states)
plot(fit_hmm_4states)
#2-state model probably over-simplifies, 4-state model probs overfits and this model had the most convergence issues
#so we'll select the 3-state model here as the best-fitting of the three

#NOTE: these are examples showing how to fit HMMs w/ covariates while also using multiple imputation to account for time series irregularity and location error. Most of the imputed model results didn't converge, which should ultimately be fixed when performing more than a simple exploratory analysis. Likewise, it is also possible to fit an HMM w/ covariates more simply by doing this for a single regularized track instead of multiple imputations.





#############################################
### Explore 3-state model more thoroughly ###
#############################################

# Explore time series of state probabilities
plotStates(fit_hmm_3states)

# Get estimates of activity budget
timeInStates(fit_hmm_3states)
#Encamped: 55%; Exploratory: 36%; Transit: 9%

# Annotate tracks w/ model results
crw_dat3 <- crw_dat2 |> 
  mutate(state_vit = fit_hmm_3states$miSum$Par$states,  #used Viterbi algorithm (most likely sequence)
         state_vit = case_when(state_vit == 1 ~ 'Encamped',
                               state_vit == 2 ~ 'Exploratory',
                               state_vit == 3 ~ 'Transit')) |> 
  cbind(fit_hmm_3states$miSum$Par$stateProbs$est) |>  #used forward-backward algorithm (probabilities of each state)
  mutate(state_fb = case_when(  #using est. from forward-backward algorithm, assign states when likely (Pr > 50%)
    Encamped > 0.5 ~ "Encamped",
    Exploratory > 0.5 ~ "Exploratory",
    Transit > 0.5 ~ "Transit",
    TRUE ~ 'Unclassified'  # Assigns NA if no value is > 0.5
  ))

# Any discrepencies between the 2 methods for assigning states?
all.equal(crw_dat3$state_vit, crw_dat3$state_fb)
#looks like there are some differences (2064 to be exact)

# Any 'Unclassified' obs?
table(crw_dat3$state_fb)  #33 obs


#-- In general, I like to use the approach of assigning states based on the confidence in the state estimates. But plenty of people use and publish the Viterbi seq. of states too --#




############################
### Viz annotated tracks ###
############################

# Plot all annotated tracks together
ggplot() +
  geom_sf(data = africa |> 
            st_transform(32736)) +
  geom_path(data = crw_dat3, aes(x, y, group = ID), alpha = 0.5, linewidth = 0.25) +
  geom_point(data = crw_dat3, aes(x, y, color = state_fb), alpha = 0.5, size = 1) +
  geom_sf(data = gl_pa |> 
            st_transform(32736), color = "black", fill = NA, linewidth = 0.25) +
  scale_color_manual("State", values = c(RColorBrewer::brewer.pal(n = 3, "Dark2"), "grey")) +
  theme_bw() +
  coord_sf(xlim = range(crw_dat3$x),
           ylim = range(crw_dat3$y))


# Facet by ID
ggplot() +
  geom_sf(data = africa |> 
            st_transform(32736)) +
  geom_path(data = crw_dat3, aes(x, y, group = ID), alpha = 0.5, linewidth = 0.25) +
  geom_point(data = crw_dat3, aes(x, y, color = state_fb), alpha = 0.5, size = 1) +
  geom_sf(data = gl_pa |> 
            st_transform(32736), color = "black", fill = NA, linewidth = 0.25) +
  scale_color_manual("State", values = c(RColorBrewer::brewer.pal(n = 3, "Dark2"), "grey")) +
  theme_bw() +
  coord_sf(xlim = range(crw_dat3$x),
           ylim = range(crw_dat3$y)) +
  facet_wrap(~ID, ncol = 2)


for (i in 1:n_distinct(crw_dat3$ID)) {
  print(
    ggplot() +
      geom_sf(data = africa |> 
                st_transform(32736)) +
      geom_path(data = crw_dat3, aes(x, y, group = ID), alpha = 0.5, linewidth = 0.25) +
      geom_point(data = crw_dat3, aes(x, y, color = state_fb), alpha = 0.5, size = 1) +
      geom_sf(data = gl_pa |> 
                st_transform(32736), color = "black", fill = NA, linewidth = 0.25) +
      scale_color_manual("State", values = c(RColorBrewer::brewer.pal(n = 3, "Dark2"), "grey")) +
      theme_bw() +
      coord_sf(xlim = range(crw_dat3$x),
               ylim = range(crw_dat3$y)) +
      ggforce::facet_wrap_paginate(~ID, ncol = 1, nrow = 1, page = i)
  )
}


# Separate by state
for (i in 1:n_distinct(crw_dat3$state_fb)) {
  print(
    ggplot() +
      geom_sf(data = africa |> 
                st_transform(32736)) +
      # geom_path(data = crw_dat3, aes(x, y, group = ID), alpha = 0.5, linewidth = 0.25) +
      geom_point(data = crw_dat3, aes(x, y, color = state_fb), alpha = 0.5, size = 1) +
      geom_sf(data = gl_pa |> 
                st_transform(32736), color = "black", fill = NA, linewidth = 0.25) +
      scale_color_manual("State", values = c(RColorBrewer::brewer.pal(n = 3, "Dark2"), "grey")) +
      theme_bw() +
      coord_sf(xlim = range(crw_dat3$x),
               ylim = range(crw_dat3$y)) +
      ggforce::facet_wrap_paginate(~state_fb, ncol = 1, nrow = 1, page = i)
  )
}


# Explore state probabilities (e.g., 'Transit' in ID 6469)
ggplot() +
  geom_path(data = crw_dat3 |> 
              filter(ID == 6469), aes(x, y, group = ID), alpha = 0.5, linewidth = 0.25) +
  geom_point(data = crw_dat3 |> 
               filter(ID == 6469), aes(x, y, color = Transit), alpha = 0.5, size = 1) +
  scale_color_distiller("Pr(Transit)", palette = "Spectral", direction = -1, limits = c(0,1)) +
  theme_bw() +
  coord_equal()



# Interactively explore results
crw_dat3 |> 
  rename(id = ID) |> 
  bayesmove::shiny_tracks(epsg = 32736)




###########################################################
### Create custom state-dependent density distrib plots ###
###########################################################

### Step lengths

# Generate the plotting data for step lengths
df_steps <- get_hmm_densities(fit_hmm_3states, metric = "step")

# Define colors (add enough for K states + Black for total)
pal <- c("#E69F00", "#56B4E9", "#009E73", "#000000")

ggplot() +
  # Add empirical histogram
  geom_histogram(data = crw_dat, aes(x = step, y = after_stat(density)), fill = "grey85",
                 color = "grey65", linewidth = 0.25) +
  geom_line(data = df_steps, aes(x = x, y = dens, color = state, linetype = state), linewidth = 1) +
  scale_colour_manual(values = pal) +
  # Assign "solid" to all states, and "dashed" to the final "Total" line
  scale_linetype_manual(values = c(rep("solid", 3), "dashed")) +
  theme_bw(base_size = 14) +
  theme(legend.title = element_blank()) +
  labs(x = "Step Length (m)", y = "Density")



### Turning angles

# Generate the plotting data for turning angles
df_angle <- get_hmm_densities(fit_hmm_3states, metric = "angle")

ggplot() +
  # Add empirical histogram
  geom_histogram(data = crw_dat, aes(x = angle, y = after_stat(density)), fill = "grey85",
                 color = "grey65", linewidth = 0.25) +
  geom_line(data = df_angle, aes(x = x, y = dens, color = state, linetype = state), linewidth = 1) +  
  scale_colour_manual(values = pal) +
  # Assign "solid" to all states, and "dashed" to the final "Total" line
  scale_linetype_manual(values = c(rep("solid", 3), "dashed")) +
  scale_x_continuous(
    limits = c(-pi, pi),
    breaks = seq(-pi, pi, by = pi/2),
    labels = expression(-pi, -frac(pi, 2), 0, frac(pi, 2), pi)
  ) +
  theme_bw(base_size = 14) +
  theme(legend.title = element_blank()) +
  labs(x = "Turning Angle (rad)", y = "Density")






########################################################
### Create custom stationary state probability plots ###
########################################################

# Generate predictions (where other covar held at its mean)
stat_probs <- plotStationary(fit_hmm_3states, plotCI = T, return = T,
                             # covs = data.frame(temperature = 25, hour = 16)  #if wanting to use diff values
                             )

# Reformat into long df
stat_probs2 <- stat_probs |> 
  map(bind_rows, .id = "state") |> 
  bind_rows() |> 
  pivot_longer(cols = c(temperature, hour), names_to = "covar", values_to = "x") |> 
  drop_na(x)

ggplot(stat_probs2, aes(x = x, y = est, color = state, fill = state)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2, color = NA) +
  geom_line(size = 1) +
  scale_color_manual("", values = pal[-4]) +
  scale_fill_manual("", values = pal[-4]) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "",
       y = "Stationary State Probability") +
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.4),
        strip.text = element_text(size = 15),
        strip.background = element_blank(),
        strip.placement = "outside"  # Moves strip text below the x-axis
        ) +
  facet_wrap(~ covar, scales = "free_x", strip.position = "bottom",
             labeller = labeller(
               covar = c("temperature" = "Temperature (°C)",
                         "hour" = "Hour of Day")
               )
             )




###############################
### Export annotated tracks ###
###############################

write_csv(crw_dat3, "processed_data/Session_3/HMM_3state_MultImp_covar.csv")
