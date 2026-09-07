
### Approximate Bayesian computation w/ INLA/inlabru ###

library(tidyverse)
library(inlabru)



#################
### Load data ###
#################

# Load tracks (used previously for RSF)
dat <- read_csv("processed_data/Session_5/dat_presabs_covars.csv") |> 
  mutate(id = as.character(id)) |> 
  filter(obs == 1)  #only keep "used" points for this analysis

# Calc step lengths
dat2 <- bayesmove::prep_data(dat, coord.names = c("x","y"), id = "id")



###########################
### Prep data for model ###
###########################

# Remove any obs w/ NAs for covars
dat3 <- dat2 |> 
  drop_na(dist2pop:ndvi)

# Z-score transform covars (i.e., center and scale by SD)
#Helps with model stability during parameter estimation
#Also easier to directly compare effect sizes across covars
dat3 <- dat3 |> 
  mutate(across(.cols = c(dist2pop:ndvi),
                .fns = ~{scale(.x) |>  #transform w/ scale() function
                    c()},
                .names = "{.col}_s"))  #create new cols for scaled covars



# Remove all records where `step = 0` (needs to be positive for Gamma distrib)
dat_bayes <- dat3 |> 
  filter(step > 0)



#####################
### Fit the model ###
#####################

#We'll use the same priors as we did in {brms}

# Define model components (these are predictors from the data)
cmp <- ~ Intercept(1) + ndvi_s + dist2pop_s + dist2water_s

# Fit model via inlabru
fit_inla <- bru(
  components = cmp,
  formula    = step ~ Intercept + ndvi_s + dist2pop_s + dist2water_s,
  family     = "gamma",  #uses log-link by default
  data       = dat_bayes,
  options    = list(
    # Set fixed-effect priors using precision (prec = 1 / SD^2)
    control.fixed = list(
      mean = 0, 
      prec = 4.0,           # b ~ normal(0, 0.5)  -> prec = 1 / (0.5^2) = 4.0
      mean.intercept = 6, 
      prec.intercept = 1.0  # Intercept ~ normal(6, 1) -> prec = 1 / (1^2) = 1.0
    ),
    # Set gamma(2, 0.1) prior on the shape parameter
    control.family = list(
      hyper = list(
        prec = list(prior = "loggamma", param = c(2, 0.1))  #needs to be loggamma for these prior values
      )
    )
  )
)


# View model results
summary(fit_inla)
#Almost identical to full Bayesian model!

#estimates from {brms} model were:
## Intercept = 6.4
## ndvi_s = 0.05
## dist2pop_s = -0.00
## dist2water = -0.07





######################################
### Explore results from posterior ###
######################################

# Draw 100 posterior samples of the linear predictor matrix (N_obs x N_samples)
eta_draws <- generate(
  fit_inla, 
  newdata = dat_bayes, 
  formula = ~ Intercept + ndvi_s + dist2pop_s + dist2water_s, 
  n.samples = 100
)

# Convert from log scale to expectation scale: mu = exp(eta)
mu_draws <- exp(eta_draws)

# Extract shape parameter estimate from INLA summary
shape_est <- fit_inla$summary.hyperpar[["mean"]]

# Simulate step lengths (y_rep) for each draw using the Gamma distribution
# In R: rgamma rate = shape / mu
y_rep <- apply(mu_draws, 2, function(mu) {
  rgamma(n = length(mu), shape = shape_est, rate = shape_est / mu)
})

# Format simulated matrix for plotting
y_rep_df <- as.data.frame(y_rep) |> 
  pivot_longer(cols = everything(), names_to = "draw", values_to = "y_rep")



# Plot Posterior Predictive Overlay (matching bayesplot style)
ggplot() +
  # Simulated draws (y_rep) in blue
  geom_density(data = y_rep_df, aes(x = y_rep, group = draw), color = "#2c7fb8", alpha = 0.05, linewidth = 0.3) +
  # Observed data (y) in black
  geom_density(data = dat_bayes, aes(x = step), color = "black", linewidth = 1) +
  scale_x_log10() +
  labs(
    title = "INLA Posterior Predictive Check",
    subtitle = "Black = Observed step lengths | Blue = Simulated draws (y_rep)",
    x = "Step Length (m, log scale)",
    y = "Density"
  ) +
  theme_minimal(base_size = 14)
#Looks essentially identical to fully Bayesian model
