
### Estimate behavioral states from HMM (discrete states) ###
## Multiple Imputation Version ##

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
            # theta = c(5,3)
            theta = c(6.855, -0.007)
            )
toc()  #took 14 sec


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



# Fit model
set.seed(2026)
#if states get flipped, try adjusting initial params or random seed number (e.g., seed 2026 gave problems w/ flipped states when using retryFits)
#alternatively, you can specify pseudo-design matrix to prevent state label switching (but doesn't work for TMB)

tic()
test_2states <- MIfitHMM(miData = crw,
                         nSims = 1,  #start w/ single imputation as a test
                         nbStates = 2,
                         dist = list(step = "gamma", angle = "wrpcauchy"),  #can use other distribs as well
                         Par0 = list(step = stepPar0, angle = anglePar0),
                         formula = ~ 1,
                         stationary = TRUE,
                         estAngleMean = list(angle=FALSE),
                         stateNames = c('ARS', 'Transit'),
                         optMethod = "TMB"
)
toc()  #took 1 min to run for 1 sim
#If problems w/ model fit or convergence, you can try 1) changing the random seed value, 2) changing the selected distributions for the movement metrics, 3) redefine new initial values for params, 4) increase number of retryFits, 5) adjust retrySD, 6) change optimization method, 7) set 'estAngleMean' to FALSE, 8) set `stationary = FALSE`, 9) specify DM, userBounds, and/or workBounds


test_2states  #looks pretty good




# Fit w/ multiple imputations and random perturbations to ensure global ML
set.seed(2026)
tic()
fit_hmm_2states <- MIfitHMM(miData = crw,
                            nSims = 10,  #fit for 10 imputations
                            ncores = 10,  #number of cores over which to process
                            nbStates = 2,
                            dist = list(step = "gamma", angle = "wrpcauchy"),  #can use other distribs as well
                            Par0 = list(step = stepPar0, angle = anglePar0),
                            formula = ~ 1,
                            stationary = TRUE,
                            estAngleMean = list(angle=FALSE),
                            stateNames = c('ARS', 'Transit'),
                            optMethod = "TMB",
                            retryFits = 10  #typically want this to be a little larger
)
toc()  #took 3 min to run for 10 sim and 10 retryFits

fit_hmm_2states
plot(fit_hmm_2states)  #error ellipses now shown for tracks, particularly at large time gaps
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
test_3states <- MIfitHMM(miData = crw,
                         nSims = 1,  #start w/ single imputation as a test
                         nbStates = 3,
                         dist = list(step = "gamma", angle = "wrpcauchy"),  #can use other distribs as well
                         Par0 = list(step = stepPar0, angle = anglePar0),
                         formula = ~ 1,
                         stationary = TRUE,
                         estAngleMean = list(angle=FALSE),
                         stateNames = c('Encamped','Exploratory','Transit'),
                         optMethod = "TMB"
)
toc()  #took 1 min to run for 1 sim

test_3states  #looks pretty good




# Fit w/ multiple imputations and random perturbations to ensure global ML
set.seed(2026)
tic()
fit_hmm_3states <- MIfitHMM(miData = crw,
                            nSims = 10,  #start w/ single imputation as a test
                            ncores = 10,  #number of cores over which to process
                            nbStates = 3,
                            dist = list(step = "gamma", angle = "wrpcauchy"),  #can use other distribs as well
                            Par0 = list(step = stepPar0, angle = anglePar0),
                            formula = ~ 1,
                            stationary = TRUE,
                            estAngleMean = list(angle=FALSE),
                            stateNames = c('Encamped','Exploratory','Transit'),
                            optMethod = "TMB",
                            retryFits = 10
)
toc()  #took 3.5 min to run for 10 sim and 10 retryFits

fit_hmm_3states
plot(fit_hmm_3states)  #error ellipses now shown for tracks, particularly at large time gaps
plotPR(fit_hmm_3states, ncores = 10)  #looks pretty good; a bit better than the 2-state model for ACF







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
test_4states <- MIfitHMM(miData = crw,
                         nSims = 1,  #start w/ single imputation as a test
                         nbStates = 4,
                         dist = list(step = "gamma", angle = "wrpcauchy"),  #can use other distribs as well
                         Par0 = list(step = stepPar0, angle = anglePar0),
                         formula = ~ 1,
                         stationary = TRUE,
                         estAngleMean = list(angle=FALSE),
                         stateNames = c('Encamped','Exploratory','Transit','Ranging'),
                         optMethod = "TMB"
)
toc()  #took 1 min to run for 1 sim

test_4states  #looks okay, but didn't converge




# Fit w/ multiple imputations and random perturbations to ensure global ML
set.seed(2026)
tic()
fit_hmm_4states <- MIfitHMM(miData = crw,
                            nSims = 10,  #start w/ single imputation as a test
                            ncores = 10,  #number of cores over which to process
                            nbStates = 4,
                            dist = list(step = "gamma", angle = "wrpcauchy"),  #can use other distribs as well
                            Par0 = list(step = stepPar0, angle = anglePar0),
                            formula = ~ 1,
                            stationary = TRUE,
                            estAngleMean = list(angle=FALSE),
                            stateNames = c('Encamped','Exploratory','Transit','Ranging'),
                            optMethod = "TMB",
                            retryFits = 10
)
toc()  #took 4.5 min to run for 10 sim and 10 retryFits

fit_hmm_4states  #looks decent
map(fit_hmm_4states$HMMfits, ~{.x$mod$code})  #8 of the 10 imputations didn't quite converge
plot(fit_hmm_4states)
plotPR(fit_hmm_4states, ncores = 10)  #also looks a bit better than the 2-state model for ACF





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
#2-state model probably over-simplifies, 4-state model probs overfits and this model didn't converge
#so we'll select the 3-state model here as the best-fitting of the three





#############################################
### Explore 3-state model more thoroughly ###
#############################################

# Explore time series of state probabilities
plotStates(fit_hmm_3states)

# Get estimates of activity budget
timeInStates(fit_hmm_3states)
#Encamped: 54%; Exploratory: 34%; Transit: 12%

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
#looks like there are some differences (2480 to be exact)

# Any 'Unclassified' obs?
table(crw_dat3$state_fb)  #1219 obs


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
my_colors <- c("#E69F00", "#56B4E9", "#009E73", "#000000")

ggplot() +
  geom_line(data = df_steps, aes(x = x, y = dens, color = state, linetype = state), linewidth = 1) +
  scale_colour_manual(values = my_colors) +
  # Assign "solid" to all states, and "dashed" to the final "Total" line
  scale_linetype_manual(values = c(rep("solid", 3), "dashed")) +
  theme_bw(base_size = 14) +
  theme(legend.title = element_blank()) +
  labs(x = "Step Length (m)", y = "Density")


### Turning angles

# Generate the plotting data for turning angles
df_angle <- get_hmm_densities(fit_hmm_3states, metric = "angle")

ggplot() +
  geom_line(data = df_angle, aes(x = x, y = dens, color = state, linetype = state), linewidth = 1) +  
  scale_colour_manual(values = my_colors) +
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




###############################
### Export annotated tracks ###
###############################

write_csv(crw_dat3, "processed_data/Session_3/HMM_3state_MultImp.csv")
