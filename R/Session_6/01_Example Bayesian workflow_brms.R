
### Example Bayesian workflow for simple model ###

library(tidyverse)
library(brms)
library(rstan)
library(bayesplot)
library(tidybayes)

# Enable parallel processing for brms
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)



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


# Check relationships among covars
cor(dat_bayes |> select(dist2pop_s, dist2water_s, ndvi_s))  #no strong corrs





##############################
### Prior predictive check ###
##############################

# Inspect default/suggested priors by {brms}
get_prior(
  step ~ ndvi_s + dist2pop_s + dist2water_s,  #this is our model formula
  data = dat_bayes,
  family = Gamma(link = "log")  #using a Gamma distribution to model step lengths
)

# Viz the impact of different types of priors
data.frame(beta = c(rnorm(n = 1000, mean = 0, sd = 5), rnorm(n = 1000, mean = 0, sd = 0.5)),
           type = rep(c("Uninformative", "Regularizing"), each = 1000)) |> 
  
  ggplot() +
  geom_density(aes(beta, color = type)) +
  theme_bw(base_size = 14)
#While previous suggestions for prior selection often used "flat" priors, this still can have a consequence on the coefficient estimates by the model. Depending on the signal from the data, a large and unreasonable value may be estimated, or the model may have difficulty converging on the posterior distribution



# Fit PRIOR-ONLY model (using default priors)
flat_priors <- c(
  prior(normal(0, 5), class = "Intercept"),         
  prior(normal(0, 5), class = "b"),                 
  prior(gamma(0.1, 0.1), class = "shape")           
)

fit_prior1 <- brm(
  formula = step ~ ndvi_s + dist2pop_s + dist2water_s,
  data = dat_bayes,
  family = Gamma(link = "log"),
  sample_prior = "only", # IGNORES DATA; draws only from prior distributions
  prior = flat_priors,
  chains = 2,
  cores = 2,
  iter = 2000,
  seed = 2026
)
#lots of divergent transitions

summary(fit_prior1)


# Prior Predictive Check
# Evaluate if defined priors produce biologically plausible step lengths
pp_check(fit_prior1, ndraws = 50) +
  scale_x_log10() +
  labs(title = "Prior Predictive Check (Flat Priors)",
       subtitle = "Are simulated step lengths biologically plausible?")
#the data (black) overlaps a bit w/ the draws from the posterior (blue), but the left tail of the posterior is VERY long an unrealistic



# Define domain-informed, weakly informative priors (to regularize estimated params)
#Log link means coefficients represent proportional changes on log scale
log(mean(dat_bayes$step))  # ~6 on log scale

reg_priors <- c(
  prior(normal(6, 1), class = "Intercept"),           # Log-mean step length ~ exp(6) ≈ 600 m
  prior(normal(0, 0.5), class = "b"),                 # Slopes constrained to plausible effect sizes
  prior(gamma(2, 0.1), class = "shape")               # Stabilize residual variance around mean
)

# Fit PRIOR-ONLY model (Sampling from Regularizing Priors)
fit_prior2 <- brm(
  formula = step ~ ndvi_s + dist2pop_s + dist2water_s,
  data = dat_bayes,
  family = Gamma(link = "log"),
  prior = reg_priors,   #supply pre-defined priors
  sample_prior = "only", # IGNORES DATA; draws only from prior distributions
  chains = 2,
  cores = 2,
  iter = 2000,
  seed = 2026
)
#no warnings about divergent transitions

summary(fit_prior2)


# Prior Predictive Check
# Evaluate if defined priors produce biologically plausible step lengths (e.g., not 10,000 km)
pp_check(fit_prior2, ndraws = 50) +
  scale_x_log10() +
  labs(title = "Prior Predictive Check (Regularizing Priors)",
       subtitle = "Are simulated step lengths biologically plausible?")
#yep, the data (black) matches pretty well w/ the draws from the posterior (blue)
#this is before we've included any data in the model




##############################################
### Fit full model (sample from posterior) ###
##############################################

fit_bayes <- brm(
  formula = step ~ ndvi_s + dist2pop_s + dist2water_s,
  data = dat_bayes,
  family = Gamma(link = "log"),
  prior = reg_priors,
  sample_prior = "yes", # Retains prior samples for posterior comparison
  chains = 4,
  cores = 4,
  iter = 2000,
  warmup = 1000,
  seed = 2026,
  # control = list(adapt_delta = 0.95) # Helps prevent divergent transitions
)
# took 26 sec



###########################################
### Check model convergence diagnostics ###
###########################################

# some of the main ones include 1) visual inspection of parameter traceplots, 2) checking the R_hat value, and 3) checking the ESS values across params

# Check Numerical Diagnostics (Rhat < 1.01, Bulk/Tail ESS > 1000)
summary(fit_bayes)
#all values look good

# Visual Convergence - Traceplots (Fuzzy caterpillars indicate convergence)
plot(fit_bayes, nvariables = 3)

# Rank Plots (Uniform distribution across chains indicates proper mixing)
mcmc_rank_overlay(fit_bayes, pars = c("b_Intercept","b_ndvi_s","b_dist2pop_s","b_dist2water_s"))





##################################
### Posterior predictive check ###
##################################

# Density Overlay (Does model simulated data match observed data distribution?)
pp_check(fit_bayes, ndraws = 100) +
  scale_x_log10() +
  labs(title = "Posterior Predictive Check", subtitle = "Observed (y) vs. Model Simulated (yrep)")
#yes, looks very close!

# Summary Statistic Checks (e.g., checking if model captures maximum step length)
pp_check(fit_bayes, type = "stat", stat = "max", ndraws = 100)
pp_check(fit_bayes, type = "stat", stat = "median", ndraws = 100)
#model currently underestimates max and median step lengths
#this may require further model tweaking to get these better aligned, such as accounting for varying effects by ID





#############################################
### Compare priors to posterior estimates ###
#############################################

# Extract posterior draws
posterior_df <- as_draws_df(fit_bayes) |> 
  select(b_ndvi_s, b_dist2water_s, b_dist2pop_s, shape, b_Intercept) |> 
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") |> 
  mutate(Distribution = "Posterior")

# Extract prior draws
prior_df <- prior_draws(fit_bayes) |> 
  # Map generic 'b' class prior draws to specific covariate names if since priors weren't named individually
  rename(b_ndvi_s = b, b_Intercept = Intercept) |> 
  mutate(b_dist2water_s = b_ndvi_s, b_dist2pop_s = b_ndvi_s) |>  #all slopes have some prior
  select(b_ndvi_s, b_dist2water_s, b_dist2pop_s, shape) |> 
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") |> 
  mutate(Distribution = "Prior")

# Combine into a single tidy data.frame
comparison_df <- bind_rows(posterior_df, prior_df)

#Trim extreme prior tail draws using Posterior IQR
comparison_trimmed <- comparison_df |> 
  group_by(parameter) |> 
  filter({
    post_vals <- value[Distribution == "Posterior"]
    p_med     <- median(post_vals)
    p_iqr     <- IQR(post_vals)
    
    # Keep values within a reasonable window around the posterior (20 x IQR)
    value >= (p_med - 20 * p_iqr) & value <= (p_med + 20 * p_iqr)
  }) |> 
  ungroup()

# Plot using "ridges" for clean vertical separation
ggplot(comparison_trimmed, aes(x = value, y = Distribution, fill = Distribution, color = Distribution)) +
  ggridges::geom_density_ridges(alpha = 0.5, scale = 0.9, rel_min_height = 0.01) +
  facet_wrap(~ parameter, scales = "free", strip.position = "bottom") +
  scale_fill_manual(values = c("Prior" = "grey70", "Posterior" = "#2c7fb8")) +
  scale_color_manual(values = c("Prior" = "grey40", "Posterior" = "#1c4d6f")) +
  labs(
    title = "Prior vs. Posterior Distributions (Focused View)",
    subtitle = "x-axis scaled to 20x Posterior IQR to reveal density curvature",
    x = "Parameter Value",
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.placement = "outside",
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )
#priors look much flatter than posterior estimates
#this means that the data were very informative when fitting the model





##############################################################
### Visualize the posterior and make ecological inferences ###
##############################################################

### Plotting Half-Eye Density Interval Distributions of Parameters

# On log (link) scale
fit_bayes |> 
  gather_draws(b_ndvi_s, b_dist2pop_s, b_dist2water_s) |> 
  
  ggplot(aes(x = .value, y = .variable, fill = .variable)) +
  stat_halfeye(point_interval = median_hdi, .width = c(0.89, 0.95)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkred") +
  labs(
    title = "Posterior Probability Distributions of Fixed Effects",
    x = "Effect Size (log-scale)",
    y = "Parameter"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")
#effect of NDVI and dist2water are strong (in opposing directions) on step lengths


# On odds (response) scale
fit_bayes |> 
  gather_draws(b_ndvi_s, b_dist2pop_s, b_dist2water_s) |> 
  mutate(.value = exp(.value)) |>  #transform from link- to response-scale
  
  ggplot(aes(x = .value, y = .variable, fill = .variable)) +
  stat_halfeye(point_interval = median_hdi, .width = c(0.89, 0.95)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "darkred") +
  labs(
    title = "Posterior Probability Distributions of Fixed Effects",
    x = "Effect Size (response-scale)",
    y = "Parameter"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")
#Suggests that mean step length decreases ~7% w/ each 1 SD (~11 km) increase in distance to water, and that mean step length increases by ~6% w/ each 1 SD (0.11) increase in NDVI
#Essentially no effect of dist2pop on mean step length



### Extracting Posterior Draws for Hypothesis Testing
posterior_draws <- as_draws_df(fit_bayes)

# Calculate Probability of Direction (p_d): Proportion of posterior > or < 0
pd_ndvi <- mean(posterior_draws$b_ndvi_s > 0)
message(sprintf("Probability of NDVI positive effect on step length: %.2f%%", pd_ndvi * 100))

pd_dist2pop <- mean(posterior_draws$b_dist2pop_s < 0)
message(sprintf("Probability of Dist to Population negative effect on step length: %.2f%%", pd_dist2pop * 100))

pd_dist2water <- mean(posterior_draws$b_dist2water_s < 0)
message(sprintf("Probability of Dist to Water negative effect on step length: %.2f%%", pd_dist2water * 100))





### Simple Conditional Marginal Effects Plots
conditional_effects(fit_bayes, effects = c("ndvi_s", "dist2water_s", "dist2pop_s"))



### Custom Conditional Marginal Effects Plots

# Generate prediction list
ce <- conditional_effects(
  fit_bayes, 
  effects = c("ndvi_s", "dist2water_s", "dist2pop_s"), 
  plot = FALSE
)

# Un-scale x-values and combine into a single tidy data.frame
effects_df <- map_dfr(names(ce), function(eff) {
  raw_var <- sub("_s$", "", eff)
  
  var_mean <- mean(dat3[[raw_var]], na.rm = TRUE)
  var_sd   <- sd(dat3[[raw_var]], na.rm = TRUE)
  
  ce[[eff]] |> 
    mutate(
      covariate = raw_var,
      x_natural = (.data[[eff]] * var_sd) + var_mean
    )
})

# Plot all marginal effects
ggplot(effects_df, aes(x = x_natural, y = estimate__)) +
  geom_line(color = "#2c7fb8", linewidth = 1) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), fill = "#2c7fb8", alpha = 0.2) +
  facet_wrap(~ covariate, scales = "free_x", strip.position = "bottom") +
  labs(
    title = "Bayesian Marginal Effects (Natural Scale)",
    y = "Expected Step Length (m)",
    x = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.placement = "outside", 
    strip.text = element_text(face = "bold")
  )
#we see essentially no change in mean step length w/ dist2pop
#avg step length ~45% lower at max recorded dist2water than at water
#avg step length ~40% higher at max recorded NDVI compared to lowest NDVI
##these results are quite similar to those of the SSF, which reflects that both models seem to be coming to the same conclusion (despite this one only evaluating step length and not step selection)

