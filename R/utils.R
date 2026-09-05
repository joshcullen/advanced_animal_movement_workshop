
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


predict_rsf_margeff <- function(fit, focal_covars, data, intercept, length.out = 100) {
  # Get all predictor variables used in the model (excluding the response)
  model_vars <- attr(terms(fit), "term.labels")
  
  # Extract the arbitrary intercept (tied to background sample weight/size)
  if (intercept) {  #check whether Intercept estimated in model
    intercept <- coef(fit)["(Intercept)"]
  } else {
    intercept <- 0
  }
  
  
  # Iterate over each focal covariate to generate predictions
  pred_df <- map(focal_covars, function(focal_var) {
    
    # 1. Identify the raw variable name (assuming the scaled variables end in "_s")
    raw_var <- sub("_s$", "", focal_var)
    
    if (!raw_var %in% names(data)) {
      stop(sprintf("Cannot find raw variable '%s' in data.", raw_var))
    }
    
    # 2. Calculate summary statistics from the raw data for back-transformation
    var_mean <- mean(data[[raw_var]], na.rm = TRUE)
    var_sd   <- sd(data[[raw_var]], na.rm = TRUE)
    
    # 3. Create new_data for prediction
    # Initialize all model variables at 0 (their scaled mean)
    new_data <- as.data.frame(matrix(0, nrow = length.out, ncol = length(model_vars)))
    names(new_data) <- model_vars
    
    # Vary only the focal variable from its scaled minimum to maximum
    focal_min <- min(data[[focal_var]], na.rm = TRUE)
    focal_max <- max(data[[focal_var]], na.rm = TRUE)
    new_data[[focal_var]] <- seq(focal_min, focal_max, length.out = length.out)
    
    # 4. Generate predictions on the link (log) scale
    preds <- predict(fit, newdata = new_data, type = "link", se.fit = TRUE)
    
    # 5. Process results into a tidy format
    res <- data.frame(
      covariate = raw_var,
      x_scaled  = new_data[[focal_var]],
      log_rss   = preds$fit - intercept
    ) |> 
      mutate(
        log_rss_lwr = log_rss - (1.96 * preds$se.fit),
        log_rss_upr = log_rss + (1.96 * preds$se.fit),
        rss         = exp(log_rss),
        rss_lwr     = exp(log_rss_lwr),
        rss_upr     = exp(log_rss_upr),
        x_natural   = (x_scaled * var_sd) + var_mean
      )
    
    return(res)
  }) |> 
    list_rbind() # Combine the list of dataframes into a single tidy dataframe
  
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
