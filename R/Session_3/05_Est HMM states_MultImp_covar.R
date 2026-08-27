
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
  ungroup() |> 
  data.frame() |>  #can't be a tibble
  select(ID, date, temperature, hour)  #only keep necessary cols

# Merge the interpolated covariates into the 'crw' object
crw_merged <- crawlMerge(crw, cov_interp, Time.name = "date")

# foo <- crw$crwPredict |> data.frame() |> filter(ID == 5605) |> select(TimeNum, date, temp)
# z <- approx(x = foo$date, y = foo$temp, xout = foo$date, rule = 2)$y


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
                         nSims = 3,  #start w/ single imputation as a test
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
                            dist = list(step = "gamma", angle = "wrpcauchy"),  #can use other distribs as well
                            Par0 = list(step = stepPar0, angle = anglePar0),
                            formula = ~ temperature * cosinor(hour, period = 24),
                            covNames = c("temperature", "hour"),  # THIS PREVENTS STRIPPING FROM OBJECT
                            stationary = FALSE,  #needs to be FALSE when using covars
                            estAngleMean = list(angle=FALSE),
                            stateNames = c('ARS', 'Transit'),
                            optMethod = "TMB"
)
toc()  #took 2.5 min to run for 10 sims

fit_hmm_2states
plot(fit_hmm_2states, plotCI = TRUE)
plotStationary(fit_hmm_2states, plotCI = TRUE)
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
                         nSims = 3,  #start w/ single imputation as a test
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

# initial parameters (obtained from 'test_2states')
Par0_3states <- getPar0(model = test_3states, formula = form)


# Fit model w/ covars
set.seed(2026)
tic()
fit_hmm_3states <- MIfitHMM(miData = crw_merged,
                            nSims = 10,  #10 imputations
                            nbStates = 3,
                            ncores = 10,
                            dist = list(step = "gamma", angle = "wrpcauchy"),  #can use other distribs as well
                            Par0 = list(step = stepPar0, angle = anglePar0),
                            formula = ~ temperature * cosinor(hour, period = 24),
                            covNames = c("temperature", "hour"),  # THIS PREVENTS STRIPPING FROM OBJECT
                            stationary = FALSE,  #needs to be FALSE when using covars
                            estAngleMean = list(angle=FALSE),
                            stateNames = c('Encamped','Exploratory','Transit'),
                            optMethod = "TMB"
)
toc()  #took 2.5 min to run for 10 sims

fit_hmm_3states
plot(fit_hmm_3states, plotCI = FALSE)
plotStationary(fit_hmm_3states, plotCI = TRUE)
plotPR(fit_hmm_3states, ncores = 5)  #looks fine
