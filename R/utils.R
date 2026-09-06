
# Function to add transformed coordinates

add_trans_coords <- function(data, coords, proj, new_proj) {
  ## data = a data.frame with coordinates to be re-projected
  ## coords = char vector; The names of the coordinates to be tranformed
  ## proj = either an EPSG code or proj string defining the projection of `coords`
  ## new_proj = an EPSG code or proj string defining the projection of the new coords to be added
  
  tmp <- sf::st_as_sf(data, coords = coords, crs = proj, remove = FALSE)
  
  if (sum(coords %in% c('x','X','y','Y')) > 0) {
    tmp |> 
      sf::st_transform(4326) |> 
      mutate(lon = unlist(map(geometry, 1)),
             lat = unlist(map(geometry, 2))) |>  #add cols from new CRS
      sf::st_drop_geometry()
  } else if (sum(coords %in% c('lon','long','longitude','Lon','Long','Longitude',
                               'lat','latitude','Lat','Latitude')) > 0) {
    tmp |> 
      sf::st_transform(new_proj) |> 
      mutate(x = unlist(map(geometry, 1)),
             y = unlist(map(geometry, 2))) |>  #add cols from new CRS
      sf::st_drop_geometry()
  } else {
    stop("Coordinate names must either be x/y or lon/lat.")
  }
  
}

#-----------------------------------

# Function to create common SpatRaster grid across IDs

create_shared_grid <- function(data, rasters, crs, res, ext = 0.3) {
  ## data = data.frame with coordinates labeled 'x' and 'y'
  ## rasters = a list object containing the SpatRasters of interest; list or raster layers should be named by ID for clarity
  ## crs = projection to be used for raster; could be EPSG code or proj4 string; For example, UTM Zone 36S could be specified as "EPSG:32736"
  ## res = numeric. Value of spatial resolution in units of specified CRS
  ## ext = the proportion by which to extend the bounding box around the data
  
  dat_bbox <- terra::ext(min(data$x), max(data$x), min(data$y), max(data$y))
  
  # Extend bbox
  dat_bbox2 <- dat_bbox * (1 + ext)
  
  # Create shared template grid
  r_template <- terra::rast(dat_bbox2, resolution = res, crs = crs)
  
  # Resample to template
  r <- lapply(rasters, \(x) terra::resample(x, r_template, method = "mean", threads = TRUE))
  
  # Replace NAs with 0 (areas where the animal didn't go should have 0 probability)
  r2 <- lapply(r, \(x) terra::subst(x, from = NA, to = 0)) |> 
    terra::rast()  #create stack
  
  # Add names
  if(!is.null(names(rasters))) {
    names(r2) <- names(rasters)
  } else if (lapply(rasters, names) |> unlist() |> unique() |> length() > 1) {
    names(r2) <- lapply(rasters, names) |> unlist()
  }
  
  return(r2)
}

#-----------------------------------

# Function to calculate space use overlap indices from SpatRasters
#works on whatever raster layers provided (e.g., full surface, 95% UD, 50% UD, etc)

calc_ud_overlap <- function(r, index, feature = NULL) {
  ## r = a SpatRaster object where the layers represent different UDs on shared grid
  ## index = character. "vi" (Volume of Intersection), "ba" (Bhattacharyya's Affinity), or "feature" (Overlap with spatial feature).
  ## feature = an sf or SpatVector polygon (required if index == "feature").
  
  
  # Normalize the SpatRaster so each layer sums to 1
  layer_sums <- terra::global(r, "sum", na.rm = TRUE)$sum
  r_norm <- r / layer_sums
  
  
  # ---------------------------------------------------------
  # Feature Overlap (Probability Mass inside a Polygon)
  # ---------------------------------------------------------
  if (index == 'feature') {
    if (is.null(feature)) {
      stop("You must provide an 'sf' or 'SpatVector' object to the 'feature' argument.")
    }
    
    # Convert sf to terra's SpatVector format
    if (inherits(feature, "sf")) {
      feature <- terra::vect(feature)
    }
    
    # Ensure CRS matches between raster and feature
    if (terra::crs(feature) != terra::crs(r_norm)) {
      feature <- terra::project(feature, terra::crs(r_norm))
    }
    
    # Mask the raster to the feature (sets all cells outside the polygon to NA)
    r_masked <- terra::mask(r_norm, feature)
    
    # Sum the remaining density values inside the polygon
    overlap_sums <- terra::global(r_masked, "sum", na.rm = TRUE)$sum
    
    overlap <- tibble(
      id = names(r_norm),
      overlap = overlap_sums
    )
    
    # ---------------------------------------------------------
    # Pairwise Overlap (VI or BA)
    # ---------------------------------------------------------
  } else if (index %in% c('vi', 'ba')) {
    
    # Extract spatial data into a list of numeric vectors for fast pairwise math
    vals_list <- as.list(r_norm) |> 
      map(~ as.numeric(terra::values(.x))) |> 
      set_names(names(r_norm))
    
    overlap <- combn(names(vals_list), 2, simplify = FALSE) |> 
      map_dfr(function(pair) {
        id1 <- pair[1]
        id2 <- pair[2]
        
        if (index == 'vi') {
          # Calculate Volume of Intersection
          val <- sum(pmin(vals_list[[id1]], vals_list[[id2]], na.rm = TRUE))
        } else if (index == 'ba') {
          # Calculate Bhattacharyya's Affinity
          val <- sum(sqrt(vals_list[[id1]] * vals_list[[id2]]), na.rm = TRUE)
        }
        
        tibble(
          from = id1,
          to = id2,
          overlap = val
        )
      })
    
  } else {
    stop("Index must be 'vi', 'ba', or 'feature'.")
  }
  
  return(overlap)
}

#-----------------------------------

### Function to test different sets of initial values in iterative manner
run_HMMs_internal = function(data, K, Par0, dist, state.names, optMethod, p, seed, ...) {
  
  set.seed(seed)
  
  # Step length
  stepMean0 <- runif(K,
                     min = Par0$step[1:K] / 2,
                     max = Par0$step[1:K] * 2)
  stepSD0 <- runif(K,
                   min = Par0$step[(K+1):(K*2)] / 2,
                   max = Par0$step[1:K] * 2)
  whichzero_sl <- which(data$step == 0)
  propzero_sl <- length(whichzero_sl)/nrow(data)
  zeromass0_sl <- c(propzero_sl, rep(0, K-1))        #for zero distances by state
  
  
  # Fit model
  if(propzero_sl > 0) {  #don't include zero mass if no 0s present
    stepPar0 <- c(stepMean0, stepSD0, zeromass0_sl)
  } else {
    stepPar0 <- c(stepMean0, stepSD0)
  }
  
  anglePar0 <- Par0$angle
  
  
  hmm.res <- fitHMM(data = data,
                    nbStates = K,
                    Par0 = list(step = stepPar0, angle = anglePar0),
                    dist = dist,
                    formula = ~ 1,
                    # stationary=TRUE, #stationary for a slightly better fit
                    # estAngleMean = list(angle=TRUE),
                    stateNames = state.names,
                    optMethod = optMethod,
                    ...
  )
  
  # Update progress bar
  p()
  
  return(hmm.res)
}


#-----------------------------------

# Function to run multiple HMMs w/ tweaked initial params for data streams
# Specific to SL and TA in current form
run_HMMs = function(data, K, Par0, dist, state.names, optMethod, niter, ncores, ...) {
  
  # Convert data.frame into list of identical elements to fit w/ different initial values
  hmm.list <- map(1:niter, function(x) data)
  
  # Fit model
  future::plan(future::multisession, workers = ncores)
  seed <- future.apply::future_lapply(seq_along(hmm.list), FUN = function(x) sample(1:5e3, 1),
                                       future.chunk.size = Inf, future.seed = TRUE)  #set seed per list element
  
  progressr::handlers(progressr::handler_progress(incomplete=".", complete="*", current="o", clear = FALSE))
  progressr::with_progress({
    #set up progress bar
    p <- progressr::progressor(steps = length(hmm.list))
    
    hmm.res <- furrr::future_map2(.x = hmm.list, .y = seed,
                                  function(x, s) {
                                    run_HMMs_internal(data = x,
                                                     K = K,
                                                     Par0 = Par0,
                                                     dist = dist,
                                                     state.names = state.names,
                                                     optMethod = optMethod,
                                                     p = p,
                                                     seed = s,
                                                     ...)  #add other args to fitHMM
                                  },
                                  .options = furrr::furrr_options(seed = TRUE))
  })
  
  future::plan(future::sequential)
  
  # Extract likelihoods of fitted models
  allnllk <- unlist(map(hmm.res, function(m) m$mod$minimum))
  
  # Sort model indices from lowest to highest negative log-likelihood
  sorted_indices <- order(allnllk, decreasing = FALSE)
  
  best_idx <- NULL
  
  # Iteratively evaluate models from best fit to worst fit
  for (idx in sorted_indices) {
    m <- hmm.res[[idx]]
    
    # 1. Optimizer convergence check (code == 0)
    code_ok <- !is.null(m$mod$code) && m$mod$code == 0
    
    # 2. Gradient checks (no Inf or NA)
    grad_ok <- !is.null(m$mod$gradient) && 
      !any(is.na(m$mod$gradient)) && 
      !any(is.infinite(m$mod$gradient))
    
    # 3. Hessian check (no NAs)
    hess_ok <- !is.null(m$mod$hessian) && !anyNA(m$mod$hessian)
    
    # 4. Parameter boundary check (flag exploding concentrations > 100)
    conc_ok <- TRUE
    if (!is.null(m$mle$angle)) {
      # von Mises concentrations above 100 indicate numerical boundary explosion
      if (any(m$mle$angle > 100, na.rm = TRUE)) {
        conc_ok <- FALSE
      }
    }
    
    # If all criteria are met, select this model and stop searching
    if (code_ok && grad_ok && hess_ok && conc_ok) {
      best_idx <- idx
      break
    }
  }
  
  # Fallback if no model passed all criteria
  if (is.null(best_idx)) {
    warning("No fitted model met all convergence and parameter sanity criteria. Returning the model with the lowest negative log-likelihood.")
    best_idx <- sorted_indices[1]
  }
  
  return(hmm.res[[best_idx]])
}

#-----------------------------------

#' Prepare state-dependent density distributions from a fitted momentuHMM or MIfitHMM model
#' 
#' @param model A fitted momentuHMM object, MIfitHMM list, or miSum element
#' @param metric Character string of the data stream (e.g., "step" or "angle")
#' @param n_points Number of points to generate for smooth density curves
get_hmm_densities <- function(model, metric = "step", n_points = 1000) {
  
  # 1. Determine model type and set up the base extraction object
  if ("miSum" %in% names(model)) {
    # User passed the full MIfitHMM list
    is_mi <- TRUE
    m_core <- model$miSum
  } else if ("MIcombine" %in% names(model)) {
    # User passed the miSum element directly
    is_mi <- TRUE
    m_core <- model
  } else {
    # User passed a standard fitHMM object
    is_mi <- FALSE
    m_core <- model
  }
  
  # 2. Extract shared components (names are identical in both object types)
  nbStates    <- length(m_core$stateNames)
  state_names <- m_core$stateNames
  dist_name   <- m_core$conditions$dist[[metric]]
  obs_data    <- m_core$data[[metric]]
  
  # 3. Extract parameter matrices and Viterbi sequences based on model type
  if (is_mi) {
    pars        <- m_core$Par$real[[metric]]$est
    viterbi_seq <- m_core$Par$states
  } else {
    pars        <- m_core$mle[[metric]]
    viterbi_seq <- momentuHMM::viterbi(m_core)
  }
  
  # 4. Calculate state frequencies
  state_freq <- table(factor(viterbi_seq, levels = 1:nbStates)) / length(viterbi_seq)
  
  # 5. Generate x-axis sequence based on observed data limits
  obs_data <- obs_data[!is.na(obs_data)]
  
  if (metric == "angle") {
    x_seq <- seq(-pi, pi, length.out = n_points)
  } else {
    # Pad the max step length slightly for a cleaner plot tail
    x_seq <- seq(0, max(obs_data, na.rm = TRUE) * 1.05, length.out = n_points)
  }
  
  # Initialize lists to store data frames
  df_list <- list()
  y_tot <- rep(0, n_points)
  
  # 6. Loop through states and calculate probability densities
  for (i in 1:nbStates) {
    
    # --- STEP LENGTH DISTRIBUTIONS ---
    if (dist_name == "gamma") {
      mu <- pars[1, i]
      sigma <- pars[2, i]
      y <- dgamma(x_seq, shape = (mu^2 / sigma^2), scale = (sigma^2 / mu))
      
    } else if (dist_name == "weibull") {
      shape <- pars[1, i]
      scale <- pars[2, i]
      y <- dweibull(x_seq, shape = shape, scale = scale)
      
      # --- TURNING ANGLE DISTRIBUTIONS ---
    } else if (dist_name == "vm") {
      # Handle cases where mean angle is fixed at 0 (estAngleMean = FALSE)
      mu <- ifelse(nrow(pars) == 1, 0, pars[1, i])
      kappa <- ifelse(nrow(pars) == 1, pars[1, i], pars[2, i])
      y <- (1 / (2 * pi * besselI(kappa, 0))) * exp(kappa * cos(x_seq - mu))
      
    } else if (dist_name == "wrpcauchy") {
      mu <- ifelse(nrow(pars) == 1, 0, pars[1, i])
      rho <- ifelse(nrow(pars) == 1, pars[1, i], pars[2, i])
      y <- (1 - rho^2) / (2 * pi * (1 + rho^2 - 2 * rho * cos(x_seq - mu)))
      
    } else {
      stop(paste("Helper function does not currently support the", dist_name, "distribution."))
    }
    
    # Scale density by the proportion of time spent in this state
    y_scaled <- y * as.numeric(state_freq[i])
    y_tot <- y_tot + y_scaled
    
    # Store in list
    df_list[[i]] <- data.frame(
      x = x_seq,
      dens = y_scaled,
      state = state_names[i]
    )
  }
  
  # 7. Append the overall "Total" density
  df_tot <- data.frame(
    x = x_seq,
    dens = y_tot,
    state = "Total"
  )
  
  # Combine everything into one tidy dataframe
  cmb <- do.call(rbind, c(df_list, list(df_tot)))
  
  # Lock factor levels so "Total" is always the last legend item
  cmb$state <- factor(cmb$state, levels = c(state_names, "Total"))
  
  return(cmb)
}



#--------------------------


predict_rsf_margeff <- function(fit, 
                                focal_covars, 
                                data, 
                                intercept = TRUE, 
                                level = "population", # "population" or "individual"
                                id_col = NULL,        # Name of ID grouping factor for individual level
                                length.out = 100) {
  
  # 1. Extract formula terms, coefficients, and covariance matrix
  is_tmb <- inherits(fit, "glmmTMB")
  
  if (is_tmb) {
    # Strip random effect (bar) terms to safely build fixed-effects design matrix X
    fixed_formula <- reformulas::nobars(formula(fit))
    cond_terms    <- delete.response(terms(fixed_formula))
    
    beta <- glmmTMB::fixef(fit)$cond
    V    <- glmmTMB:::vcov.glmmTMB(fit)$cond
  } else {
    cond_terms <- delete.response(terms(fit))
    beta <- coef(fit)
    V    <- vcov(fit)
  }
  
  # Remove NA coefficients if present
  valid_coefs <- !is.na(beta)
  beta <- beta[valid_coefs]
  V    <- V[names(beta), names(beta), drop = FALSE]
  
  # Identify raw variables required for fixed terms
  required_vars <- all.vars(cond_terms)
  
  # 2. Setup individual-level parameters if requested
  if (level == "individual") {
    if (!is_tmb) {
      stop("Individual-level predictions are currently designed for glmmTMB mixed models.")
    }
    
    # Auto-detect ID column if not specified
    if (is.null(id_col)) {
      ran_names <- names(glmmTMB::ranef(fit)$cond)
      if (length(ran_names) > 0) {
        id_col <- ran_names[1]
        message(sprintf("Auto-detected ID column: '%s'", id_col))
      } else {
        stop("No random effects found in the model to generate individual predictions.")
      }
    }
    
    # Extract individual-level coefficients (Fixed + Random)
    indiv_coef_df <- coef(fit)$cond[[id_col]]
    if (is.null(indiv_coef_df)) {
      stop(sprintf("Could not extract random coefficients for group '%s'.", id_col))
    }
  }
  
  # 3. Iterate over focal covariates
  pred_df <- map(focal_covars, function(focal_var) {
    
    raw_var <- sub("_s$", "", focal_var)
    
    if (!raw_var %in% names(data)) {
      stop(sprintf("Cannot find raw variable '%s' in data.", raw_var))
    }
    
    var_mean <- mean(data[[raw_var]], na.rm = TRUE)
    var_sd   <- sd(data[[raw_var]], na.rm = TRUE)
    
    focal_min <- min(data[[focal_var]], na.rm = TRUE)
    focal_max <- max(data[[focal_var]], na.rm = TRUE)
    focal_seq <- seq(focal_min, focal_max, length.out = length.out)
    
    # Build prediction dataframe with non-focal variables set to 0 (their scaled mean)
    df_new <- as.data.frame(matrix(0, nrow = length.out, ncol = length(required_vars)))
    names(df_new) <- required_vars
    df_new[[focal_var]] <- focal_seq
    
    # Generate Model Matrix X
    X <- model.matrix(cond_terms, data = df_new)
    common_cols <- intersect(colnames(X), names(beta))
    X <- X[, common_cols, drop = FALSE]
    
    # Zero out the intercept column so it doesn't inflate SEs
    if ("(Intercept)" %in% colnames(X)) {
      X[, "(Intercept)"] <- 0
    }
    
    # --- A. POPULATION-LEVEL PREDICTIONS ---
    if (level == "population") {
      beta_sub <- beta[common_cols]
      V_sub    <- V[common_cols, common_cols, drop = FALSE]
      
      # if (intercept && "(Intercept)" %in% names(beta_sub)) {
      #   int_val <- beta_sub["(Intercept)"]
      # } else {
      #   int_val <- 0
      # }
      
      # Linear predictor (X %*% beta) and Standard Errors sqrt(diag(X %*% V %*% t(X)))
      fit_link <- as.vector(X %*% beta_sub)
      se_link  <- sqrt(rowSums((X %*% V_sub) * X))
      
      res <- data.frame(
        covariate = raw_var,
        x_scaled  = focal_seq,
        log_rss   = fit_link,
        se_link   = se_link
      ) |> 
        mutate(
          log_rss_lwr = log_rss - (1.96 * se_link),
          log_rss_upr = log_rss + (1.96 * se_link),
          rss         = exp(log_rss),
          rss_lwr     = exp(log_rss_lwr),
          rss_upr     = exp(log_rss_upr),
          x_natural   = (x_scaled * var_sd) + var_mean
        )
      
      return(res)
      
      # --- B. INDIVIDUAL-LEVEL PREDICTIONS ---
    } else if (level == "individual") {
      
      id_levels <- rownames(indiv_coef_df)
      
      indiv_res_list <- map(id_levels, function(id_i) {
        # Extract individual specific coefficients (Fixed + Random)
        b_i <- as.numeric(indiv_coef_df[id_i, common_cols])
        names(b_i) <- common_cols
        
        # if (intercept && "(Intercept)" %in% common_cols) {
        #   int_val_i <- b_i["(Intercept)"]
        # } else {
        #   int_val_i <- 0
        # }
        
        # Individual linear predictor
        fit_link_i <- as.vector(X %*% b_i)
        
        data.frame(
          id        = id_i,
          covariate = raw_var,
          x_scaled  = focal_seq,
          log_rss   = fit_link_i
        ) |> 
          mutate(
            rss       = exp(log_rss),
            x_natural = (x_scaled * var_sd) + var_mean
          )
      }) |> list_rbind()
      
      return(indiv_res_list)
    }
    
  }) |> list_rbind()
  
  return(pred_df)
}


#----------------------------


predict_gam_margeff <- function(fit, data) {
  
  # 1. Extract smooth estimates on the link scale
  # This returns .estimate (log-RSS) and .se for each smooth term
  sm <- gratia::smooth_estimates(fit)
  
  # 2. Parse covariate names from the smooth term
  # e.g., "s(dist2pop_s)" becomes "dist2pop_s", and the raw becomes "dist2pop"
  sm <- sm |>
    mutate(
      covar_scaled = sub("^s\\(([^)]+)\\).*", "\\1", .smooth),
      covar_raw    = sub("_s$", "", covar_scaled)
    )
  
  # 3. Process estimates and exponentiate to get RSS
  sm_processed <- sm |>
    mutate(
      # gratia creates a column for each covariate; extract the relevant x_scaled value dynamically
      x_scaled = purrr::map2_dbl(covar_scaled, row_number(), ~ sm[[.x]][.y]),
      
      # Calculate confidence intervals on link scale
      log_rss_lwr = .estimate - (1.96 * .se),
      log_rss_upr = .estimate + (1.96 * .se),
      
      # Exponentiate to relative selection strength
      rss     = exp(.estimate),
      rss_lwr = exp(log_rss_lwr),
      rss_upr = exp(log_rss_upr)
    )
  
  # 4. Back-transform the x-axis to natural scale
  final_df <- sm_processed |>
    group_by(covar_raw) |>
    mutate(
      var_mean  = mean(data[[covar_raw[1]]], na.rm = TRUE),
      var_sd    = sd(data[[covar_raw[1]]], na.rm = TRUE),
      x_natural = (x_scaled * var_sd) + var_mean
    ) |>
    ungroup() |>
    select(covariate = covar_raw, x_scaled, x_natural, rss, rss_lwr, rss_upr)
  
  return(final_df)
}


#-----------------------------------


# Function to predict RSS from a fitted simple SSF in {amt}
predict_ssf_margeff <- function(fit, focal_covars, data, length.out = 100) {
  
  # 1. Extract model terms excluding the strata term
  # model_obj  <- if (inherits(fit, "amt_fit")) fit$model else fit
  model_obj  <- fit$model
  all_terms  <- attr(terms(model_obj), "term.labels")
  model_vars <- all_terms[!grepl("^strata\\(", all_terms)]
  
  # 2. Build x0 (reference location with all scaled covariates held at 0)
  x0 <- as.data.frame(matrix(0, nrow = 1, ncol = length(model_vars)))
  names(x0) <- model_vars
  
  # 3. Iterate over focal covariates to generate log-RSS predictions
  pred_df <- map(focal_covars, function(focal_var) {
    
    # Identify raw variable name for back-transformation
    raw_var <- sub("_s$", "", focal_var)
    if (!raw_var %in% names(data)) {
      stop(sprintf("Cannot find raw variable '%s' in data.", raw_var))
    }
    
    var_mean <- mean(data[[raw_var]], na.rm = TRUE)
    var_sd   <- sd(data[[raw_var]], na.rm = TRUE)
    
    # Create x1 grid: vary focal covariate, hold all others at 0
    x1 <- as.data.frame(matrix(0, nrow = length.out, ncol = length(model_vars)))
    names(x1) <- model_vars
    
    focal_min <- min(data[[focal_var]], na.rm = TRUE)
    focal_max <- max(data[[focal_var]], na.rm = TRUE)
    x1[[focal_var]] <- seq(focal_min, focal_max, length.out = length.out)
    
    # Calculate log-RSS via {amt}
    lr <- amt::log_rss(fit, x1 = x1, x2 = x0, ci = "se")
    
    # Process results and back-transform x-axis
    res <- lr$df |> 
      mutate(
        covariate = raw_var,
        x_scaled  = x1[[focal_var]],
        x_natural = (x_scaled * var_sd) + var_mean,
        rss       = exp(log_rss),
        rss_lwr   = exp(lwr),
        rss_upr   = exp(upr)
      ) |> 
      select(covariate, x_scaled, x_natural, log_rss, rss, rss_lwr, rss_upr)
    
    return(res)
  }) |> 
    list_rbind()
  
  return(pred_df)
}


#------------------------


# Function to predict RSS from a fitted iSSF in {amt}
predict_issf_margeff <- function(fit, focal_covars, data_steps, data_raw = NULL, length.out = 100) {
  
  # 1. Extract base model variables (excluding response and strata)
  # model_obj <- if (inherits(fit, "amt_fit")) fit$model else fit
  model_obj <- fit$model
  raw_vars  <- all.vars(formula(model_obj))
  base_vars <- setdiff(raw_vars, c("case_", "step_id_", "strata", "id"))
  
  # 2. Build x0 baseline frame with realistic defaults
  x0 <- as.data.frame(matrix(0, nrow = 1, ncol = length(base_vars)))
  names(x0) <- base_vars
  
  # Set reference defaults for movement variables if present
  if ("sl_" %in% base_vars)     x0[["sl_"]]     <- mean(data_steps[["sl_"]], na.rm = TRUE)
  if ("log_sl_" %in% base_vars) x0[["log_sl_"]] <- log(mean(data_steps[["sl_"]], na.rm = TRUE))
  if ("cos_ta_" %in% base_vars) x0[["cos_ta_"]] <- mean(data_steps[["cos_ta_"]], na.rm = TRUE)
  
  # 3. Iterate over focal covariates
  pred_df <- map(focal_covars, function(focal_var) {
    
    raw_var <- sub("_s$", "", focal_var)
    
    # Initialize x1 as a copy of x0 across all rows
    x1 <- x0[rep(1, length.out), ]
    
    # Generate gradient sequence for focal variable
    focal_min <- min(data_steps[[focal_var]], na.rm = TRUE)
    focal_max <- max(data_steps[[focal_var]], na.rm = TRUE)
    focal_seq <- seq(focal_min, focal_max, length.out = length.out)
    x1[[focal_var]] <- focal_seq
    
    # Synchronize log_sl_ if sl_ is the focal variable
    if (focal_var == "sl_" && "log_sl_" %in% base_vars) {
      x1[["log_sl_"]] <- log(focal_seq)
    }
    
    # Calculate log-RSS (amt handles interaction columns internally)
    lr <- amt::log_rss(fit, x1 = x1, x2 = x0, ci = "se")
    
    # Back-transform x-axis to natural scale if raw dataset is provided
    if (!is.null(data_raw) && raw_var %in% names(data_raw)) {
      var_mean  <- mean(data_raw[[raw_var]], na.rm = TRUE)
      var_sd    <- sd(data_raw[[raw_var]], na.rm = TRUE)
      x_natural <- (focal_seq * var_sd) + var_mean
    } else {
      x_natural <- focal_seq
    }
    
    res <- lr$df |> 
      mutate(
        covariate = raw_var,
        x_scaled  = focal_seq,
        x_natural = x_natural,
        rss       = exp(log_rss),
        rss_lwr   = exp(lwr),
        rss_upr   = exp(upr)
      ) |> 
      select(covariate, x_scaled, x_natural, log_rss, rss, rss_lwr, rss_upr)
    
    return(res)
  }) |> 
    list_rbind()
  
  return(pred_df)
}


#-----------------------------


# Function for properly handling GLM(M) when making spatial prediction w/ terra::predict()
predict_rsf_raster <- function(model, data, type = "all", ...) {
  
  # 1. Extract fixed coefficients and covariance matrix
  if (inherits(model, "glmmTMB")) {
    beta   <- glmmTMB::fixef(model)$cond
    V_full <- glmmTMB:::vcov.glmmTMB(model)$cond
  } else {
    beta   <- coef(model)
    V_full <- vcov(model)
  }
  
  # Exclude intercept to focus on relative selection strength
  beta_covars <- beta[names(beta) != "(Intercept)"]
  cov_names   <- names(beta_covars)
  
  # Extract covariance matrix for non-intercept slope terms
  V <- V_full[cov_names, cov_names, drop = FALSE]
  
  # 2. Extract raster chunk design matrix X
  X <- as.matrix(data[, cov_names, drop = FALSE])
  
  # 3. Calculate log(RSS) point estimate
  log_rss <- as.vector(X %*% beta_covars)
  
  # 4. Standard Error on the log scale: sqrt(diag(X %*% V %*% t(X)))
  se_link <- sqrt(rowSums((X %*% V) * X))
  
  # 5. Return outputs based on specified type
  if (type == "rss") {
    return(exp(log_rss))
    
  } else if (type == "se_link") {
    return(se_link) # Spatial SD/SE on the log-RSS scale
    
  } else if (type == "se_rss") {
    return(exp(log_rss) * se_link) # Delta-method SE on exponential scale
    
  } else if (type == "all") {
    # Returns 4 columns ->terra outputs a 4-layer raster
    log_lwr <- log_rss - (1.96 * se_link)
    log_upr <- log_rss + (1.96 * se_link)
    
    out <- cbind(
      rss     = exp(log_rss),
      se_link = se_link,
      rss_lwr = exp(log_lwr),
      rss_upr = exp(log_upr)
    )
    return(out)
  }
}


#-----------------------------


# Function for properly handling GAM when making spatial prediction w/ terra::predict()
predict_gam_rsf <- function(model, data, type = "all") {
  df_chunk <- as.data.frame(data)
  
  # Generate design matrix X for all smooth and parametric terms
  X <- mgcv::predict.gam(model, newdata = df_chunk, type = "lpmatrix")
  
  # Zero out the intercept column to isolate relative selection strength
  if ("(Intercept)" %in% colnames(X)) {
    X[, "(Intercept)"] <- 0
  }
  
  beta <- coef(model)
  V    <- vcov(model)
  
  log_rss <- as.vector(X %*% beta)
  se_link <- sqrt(rowSums((X %*% V) * X))
  
  if (type == "rss") return(exp(log_rss))
  if (type == "all") {
    return(cbind(
      rss     = exp(log_rss),
      se_link = se_link,
      rss_lwr = exp(log_rss - 1.96 * se_link),
      rss_upr = exp(log_rss + 1.96 * se_link)
    ))
  }
}

