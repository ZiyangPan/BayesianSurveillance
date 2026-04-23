# ============================================================
# File: helper/utils_real_data.R
# Purpose:
#   Helper functions for the profile-based real-data analysis
#   under Bayesian latent factor modeling and bias correction.
# ============================================================

# ------------------------------------------------------------
# Function: extract_profile_data
# Purpose:
#   Validate and extract the profile-level object used in the
#   real-data analysis.
# Input:
#   profile_dat - profile-level data.frame.
# Output:
#   Returns a validated data.frame.
# ------------------------------------------------------------
extract_profile_data <- function(profile_dat) {
  required_cols <- c(
    "exposure_id", "outcome_id", "period_id",
    "point", "value"
  )
  
  assert_required_columns(profile_dat, required_cols, object_name = "profile-level object")
  profile_dat
}

# ------------------------------------------------------------
# Function: create_loglikelihood_function
# Purpose:
#   Create a profile log-likelihood function using spline
#   interpolation within range and quadratic extrapolation
#   outside range, matching the original real-data code logic.
# Input:
#   lik - data.frame with columns 'point' and 'value'.
# Output:
#   Function that takes numeric beta and returns a named list
#   with elements 'value' and 'gradient'.
# ------------------------------------------------------------
create_loglikelihood_function <- function(lik) {
  
  # --- 1. Input Validation ---
  if (!is.data.frame(lik) || !all(c("point", "value") %in% names(lik))) {
    stop("Input 'lik' must be a data frame with 'point' and 'value' columns.")
  }
  if (nrow(lik) < 4) {
    stop("Input 'lik' must have at least 4 rows for stable interpolation and extrapolation.")
  }
  
  # --- 2. Prepare Data ---
  lik_sorted <- lik[order(lik$point), ]
  
  min_beta <- min(lik_sorted$point)
  max_beta <- max(lik_sorted$point)
  
  # --- 3. Set up Interpolation ---
  spline_interpolator <- splinefun(lik_sorted$point, lik_sorted$value)
  
  # --- 4. Set up Extrapolation ---
  extrap_data_left <- head(lik_sorted, 3)
  model_left <- lm(value ~ point + I(point^2), data = extrap_data_left)
  coefs_left <- coef(model_left)
  
  extrap_data_right <- tail(lik_sorted, 3)
  model_right <- lm(value ~ point + I(point^2), data = extrap_data_right)
  coefs_right <- coef(model_right)
  
  # --- 5. Define and Return the Main Function ---
  loglikelihood <- function(beta) {
    
    if (!is.numeric(beta)) {
      stop("Input 'beta' must be a numeric value or vector.")
    }
    
    n <- length(beta)
    ll_value <- numeric(n)
    ll_gradient <- numeric(n)
    
    idx_interp <- beta >= min_beta & beta <= max_beta
    idx_extrap_left <- beta < min_beta
    idx_extrap_right <- beta > max_beta
    
    # Interpolation
    if (any(idx_interp)) {
      beta_interp <- beta[idx_interp]
      ll_value[idx_interp] <- spline_interpolator(beta_interp)
      ll_gradient[idx_interp] <- spline_interpolator(beta_interp, deriv = 1)
    }
    
    # Left extrapolation
    if (any(idx_extrap_left)) {
      beta_extrap_left <- beta[idx_extrap_left]
      new_data_left <- data.frame(point = beta_extrap_left)
      ll_value[idx_extrap_left] <- predict(model_left, newdata = new_data_left)
      ll_gradient[idx_extrap_left] <- coefs_left[2] + 2 * coefs_left[3] * beta_extrap_left
    }
    
    # Right extrapolation
    if (any(idx_extrap_right)) {
      beta_extrap_right <- beta[idx_extrap_right]
      new_data_right <- data.frame(point = beta_extrap_right)
      ll_value[idx_extrap_right] <- predict(model_right, newdata = new_data_right)
      ll_gradient[idx_extrap_right] <- coefs_right[2] + 2 * coefs_right[3] * beta_extrap_right
    }
    
    return(list(value = ll_value, gradient = ll_gradient))
  }
  
  return(loglikelihood)
}

# ------------------------------------------------------------
# Function: create_loglikelihood_batch_function
# Purpose:
#   Create a batch profile log-likelihood function from a data
#   frame containing one row per profile, matching the original
#   batch implementation logic.
# Input:
#   lik_batch_df - data.frame containing profile identifiers
#                  and columns 'point' and 'value'.
# Output:
#   Function that applies all profile interpolators to beta.
#   The returned function stores 'interpolator_list' in its
#   environment.
# ------------------------------------------------------------
create_loglikelihood_batch_function <- function(lik_batch_df) {
  
  # --- 1. Input Validation ---
  potential_id_cols <- c("exposure_id", "outcome_id", "period_id")
  id_cols_present <- intersect(potential_id_cols, names(lik_batch_df))
  
  if (length(id_cols_present) == 0) {
    stop(paste(
      "Input data frame must contain at least one of:",
      paste(potential_id_cols, collapse = ", ")
    ))
  }
  
  required_data_cols <- c("point", "value")
  if (!all(required_data_cols %in% names(lik_batch_df))) {
    stop("Input data frame is missing 'point' and/or 'value' columns.")
  }
  
  # --- 2. Pre-process and Create All Interpolator Functions ---
  keys <- do.call(paste, c(lik_batch_df[id_cols_present], sep = "_"))
  
  points_list <- lapply(strsplit(lik_batch_df$point, ";"), as.numeric)
  values_list <- lapply(strsplit(lik_batch_df$value, ";"), as.numeric)
  
  interpolator_list <- mapply(
    function(p, v, key) {
      lik_single <- data.frame(point = p, value = v)
      
      if (any(is.na(p)) || any(is.na(v)) || nrow(lik_single) < 4) {
        warning(paste("Skipping profile", key, "due to invalid data or insufficient points."))
        return(NULL)
      }
      
      return(create_loglikelihood_function(lik_single))
    },
    p = points_list,
    v = values_list,
    key = keys,
    SIMPLIFY = FALSE
  )
  
  names(interpolator_list) <- keys
  interpolator_list <- interpolator_list[!sapply(interpolator_list, is.null)]
  
  # --- 3. Define and Return the Main Batch Function ---
  loglikelihoodBatch <- function(beta) {
    results_list <- lapply(interpolator_list, function(func) {
      func(beta)
    })
    return(results_list)
  }
  
  return(loglikelihoodBatch)
}

# ------------------------------------------------------------
# Function: select_fixed_pairs
# Purpose:
#   Construct fixed exposures and fixed exposure-outcome pairs
#   over a user-specified period range.
# Input:
#   profile_dat - validated profile-level data.frame.
#   period_range - integer vector of requested periods.
# Output:
#   List containing exposure_ids_fixed, outcomes_fixed,
#   fixed_pairs, outcome_counts, and filtered_profile_dat.
# ------------------------------------------------------------
select_fixed_pairs <- function(
    profile_dat,
    period_range
) {
  profiles <- profile_dat[
    profile_dat$period_id %in% period_range,
  ]
  
  if (nrow(profiles) == 0) {
    stop("No rows left after filtering. Check the input profile object and period_range.")
  }
  
  # 1) fixed exposures appearing in all periods
  exp_by_period <- unique(profiles[, c("period_id", "exposure_id")])
  exp_counts <- aggregate(
    period_id ~ exposure_id,
    data = exp_by_period,
    FUN = function(x) length(unique(x))
  )
  names(exp_counts)[2] <- "n_period"
  
  exposure_ids_fixed <- exp_counts$exposure_id[
    exp_counts$n_period == length(unique(period_range))
  ]
  exposure_ids_fixed <- sort(exposure_ids_fixed)
  
  if (length(exposure_ids_fixed) == 0) {
    stop(
      "No exposure_id appears in all periods ",
      paste(period_range, collapse = ","),
      "."
    )
  }
  
  # 2) for each period, outcomes common to all fixed exposures
  common_outcomes_list <- lapply(period_range, function(p) {
    tmp <- profiles[
      profiles$period_id == p &
        profiles$exposure_id %in% exposure_ids_fixed,
      c("exposure_id", "outcome_id")
    ]
    tmp <- unique(tmp)
    
    outcomes_by_exp <- split(tmp$outcome_id, tmp$exposure_id)
    
    if (length(outcomes_by_exp) < length(exposure_ids_fixed)) {
      return(integer(0))
    }
    
    Reduce(intersect, outcomes_by_exp)
  })
  names(common_outcomes_list) <- as.character(period_range)
  
  common_long_list <- lapply(seq_along(period_range), function(i) {
    p <- period_range[i]
    outs <- common_outcomes_list[[i]]
    if (length(outs) == 0) return(NULL)
    data.frame(period_id = p, outcome_id = outs, stringsAsFactors = FALSE)
  })
  common_long <- do.call(rbind, common_long_list)
  
  if (is.null(common_long) || nrow(common_long) == 0) {
    stop("No outcomes satisfy the strict common-across-exposures condition in any period.")
  }
  
  # 3) count coverage across periods
  common_long <- unique(common_long)
  outcome_counts <- aggregate(
    period_id ~ outcome_id,
    data = common_long,
    FUN = function(x) length(unique(x))
  )
  names(outcome_counts)[2] <- "n_period"
  outcome_counts <- outcome_counts[order(-outcome_counts$n_period, outcome_counts$outcome_id), ]
  
  outcomes_fixed <- outcome_counts$outcome_id[
    outcome_counts$n_period == length(unique(period_range))
  ]
  outcomes_fixed <- sort(outcomes_fixed)
  
  if (length(outcomes_fixed) == 0) {
    stop(
      "No outcome_id has full coverage across all periods ",
      paste(period_range, collapse = ","),
      " for all fixed exposures."
    )
  }
  
  # 4) build all fixed pairs
  fixed_pairs <- expand.grid(
    exposure_id = exposure_ids_fixed,
    outcome_id = outcomes_fixed,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  fixed_pairs <- fixed_pairs[order(fixed_pairs$exposure_id, fixed_pairs$outcome_id), ]
  
  # 5) check
  pair_period_check <- profiles[
    profiles$exposure_id %in% exposure_ids_fixed &
      profiles$outcome_id %in% outcomes_fixed,
    c("period_id", "exposure_id", "outcome_id")
  ]
  pair_period_check <- unique(pair_period_check)
  
  pair_counts <- aggregate(
    period_id ~ exposure_id + outcome_id,
    data = pair_period_check,
    FUN = function(x) length(unique(x))
  )
  names(pair_counts)[3] <- "n_period"
  
  bad_pairs <- pair_counts[pair_counts$n_period < length(unique(period_range)), ]
  
  if (nrow(bad_pairs) > 0) {
    stop(
      "Some fixed (exposure, outcome) pairs do not appear in all periods. ",
      "This should not happen under the strict common-outcome construction."
    )
  }
  
  list(
    exposure_ids_fixed = exposure_ids_fixed,
    outcomes_fixed = outcomes_fixed,
    fixed_pairs = fixed_pairs,
    outcome_counts = outcome_counts,
    filtered_profile_dat = profiles
  )
}

# ------------------------------------------------------------
# Function: prepare_period_data
# Purpose:
#   Prepare the fixed-exposure/common-outcome analysis set for
#   one user-specified period.
# Input:
#   profile_dat - validated profile-level data.frame.
#   period_id - scalar analysis period.
#   exposure_ids_fixed - fixed exposures.
#   ooi_ids - user-specified OOI ids.
# Output:
#   List containing period-specific data, outcome sets,
#   jo_to_key matrices, and profile lookup objects.
# ------------------------------------------------------------
prepare_period_data <- function(
    profile_dat,
    period_id,
    exposure_ids_fixed,
    ooi_ids
) {
  profiles <- profile_dat[
    profile_dat$period_id == period_id,
  ]
  
  profiles <- profiles[profiles$exposure_id %in% exposure_ids_fixed, ]
  
  if (nrow(profiles) == 0L) {
    stop(sprintf("No profile records are available for period %s after filtering.", period_id))
  }
  
  outcomes_by_exposure <- split(profiles$outcome_id, profiles$exposure_id)
  common_outcomes <- Reduce(intersect, outcomes_by_exposure)
  
  if (length(common_outcomes) == 0L) {
    stop(sprintf(
      "No common outcomes remain across the selected fixed exposures in period %s.",
      period_id
    ))
  }
  
  check_ooi_ids(ooi_ids, common_outcomes)
  
  outcome_ooi <- sort(unique(ooi_ids))
  outcome_nc <- sort(setdiff(common_outcomes, outcome_ooi))
  
  if (length(outcome_nc) == 0L) {
    stop(sprintf(
      "No negative control outcomes remain after removing OOI in period %s.",
      period_id
    ))
  }
  
  outcome_all <- c(outcome_nc, outcome_ooi)
  O <- length(outcome_all)
  J <- length(exposure_ids_fixed)
  
  outcome_index <- setNames(seq_len(O), outcome_all)
  
  profiles_indexed <- profiles[
    profiles$outcome_id %in% outcome_all,
  ]
  profiles_indexed$j <- match(profiles_indexed$exposure_id, exposure_ids_fixed)
  profiles_indexed$o <- outcome_index[as.character(profiles_indexed$outcome_id)]
  profiles_indexed <- profiles_indexed[order(profiles_indexed$j, profiles_indexed$o), ]
  
  dup_check <- aggregate(
    rep(1, nrow(profiles_indexed)) ~ j + o,
    data = profiles_indexed,
    FUN = length
  )
  names(dup_check)[3] <- "n"
  
  if (any(dup_check$n > 1)) {
    stop("Multi mapping detected. Check again.")
  }
  
  batch_data <- profiles_indexed[, c("exposure_id", "outcome_id", "period_id", "point", "value")]
  
  loglik_batch <- create_loglikelihood_batch_function(batch_data)
  prof_env <- environment(loglik_batch)
  prof_fun_list <- prof_env$interpolator_list
  
  batch_data$key <- do.call(
    paste,
    c(batch_data[, c("exposure_id", "outcome_id", "period_id")], sep = "_")
  )
  
  jo_to_key_all <- matrix(NA_character_, nrow = J, ncol = O)
  for (r in seq_len(nrow(batch_data))) {
    j_r <- profiles_indexed$j[r]
    o_r <- profiles_indexed$o[r]
    jo_to_key_all[j_r, o_r] <- batch_data$key[r]
  }
  
  if (any(is.na(jo_to_key_all))) {
    stop(
      "Missing profiles detected in jo_to_key at period_id = ",
      period_id,
      ". This violates fixed exposure/OOI design."
    )
  }
  
  jo_to_key_nc <- jo_to_key_all[, seq_len(length(outcome_nc)), drop = FALSE]
  jo_to_key_ooi <- jo_to_key_all[, (length(outcome_nc) + 1L):O, drop = FALSE]
  
  list(
    period_profile_dat = profiles_indexed,
    exposure_ids = exposure_ids_fixed,
    outcome_nc = outcome_nc,
    outcome_ooi = outcome_ooi,
    outcome_all = outcome_all,
    jo_to_key_nc = jo_to_key_nc,
    jo_to_key_ooi = jo_to_key_ooi,
    jo_to_key_all = jo_to_key_all,
    loglik_batch = loglik_batch,
    prof_fun_list = prof_fun_list
  )
}

# ------------------------------------------------------------
# Function: stage1_profile_loglik_one_outcome
# Purpose:
#   Compute the NC profile log-likelihood contribution for one
#   outcome given the J-vector beta_col.
# Input:
#   beta_col - numeric vector of length J.
#   key_col - character vector of length J.
#   prof_fun_list - named list of profile likelihood functions.
# Output:
#   Numeric scalar log-likelihood.
# ------------------------------------------------------------
stage1_profile_loglik_one_outcome <- function(beta_col, key_col, prof_fun_list) {
  ll <- 0
  for (j in seq_along(beta_col)) {
    key <- key_col[j]
    if (!key %in% names(prof_fun_list)) {
      stop(sprintf("Missing profile likelihood for key '%s'.", key))
    }
    out <- prof_fun_list[[key]](beta_col[j])
    ll <- ll + as.numeric(out$value)
  }
  ll
}

# ------------------------------------------------------------
# Function: run_stage1_realdata
# Purpose:
#   Run Stage 1 MCMC for negative control profile likelihoods.
# Input:
#   jo_to_key_nc - J x m key matrix.
#   prof_fun_list - named list of profile likelihood functions.
#   n_iter, burn_in - MCMC settings.
#   sigma_b - proposal SD for b_{j,o}.
#   sd_mu0 - prior SD for mu_j.
#   a0, b0 - inverse-gamma prior parameters for tau_j^2.
# Output:
#   List containing full samples, posterior summaries based on
#   post-burn-in iterations, and acceptance rates.
# ------------------------------------------------------------
run_stage1_realdata <- function(
    jo_to_key_nc,
    prof_fun_list,
    n_iter = 5000,
    burn_in = 1000,
    sigma_b = 0.63,
    sd_mu0 = 0.20,
    a0 = 2.5,
    b0 = 0.06
) {
  J <- nrow(jo_to_key_nc)
  m <- ncol(jo_to_key_nc)
  
  mu <- stats::rnorm(J, 0, sd_mu0)
  tau2 <- 1 / stats::rgamma(J, shape = a0, rate = b0)
  
  b_mat <- matrix(0, nrow = J, ncol = m)
  for (j in seq_len(J)) {
    b_mat[j, ] <- stats::rnorm(m, mean = mu[j], sd = sqrt(tau2[j]))
  }
  
  mu_samples <- matrix(NA_real_, nrow = n_iter, ncol = J)
  tau2_samples <- matrix(NA_real_, nrow = n_iter, ncol = J)
  tau_samples <- matrix(NA_real_, nrow = n_iter, ncol = J)
  b_samples <- array(NA_real_, dim = c(n_iter, J, m))
  
  accept_b <- matrix(0, nrow = J, ncol = m)
  
  for (iter in seq_len(n_iter)) {
    for (j in seq_len(J)) {
      prec <- m / tau2[j] + 1 / (sd_mu0^2)
      var_j <- 1 / prec
      mean_j <- var_j * (sum(b_mat[j, ]) / tau2[j])
      mu[j] <- stats::rnorm(1, mean_j, sqrt(var_j))
    }
    
    for (j in seq_len(J)) {
      ss <- sum((b_mat[j, ] - mu[j])^2)
      a_post <- a0 + m / 2
      b_post <- b0 + 0.5 * ss
      tau2[j] <- 1 / stats::rgamma(1, shape = a_post, rate = b_post)
    }
    
    for (o in seq_len(m)) {
      key_col <- jo_to_key_nc[, o]
      
      for (j in seq_len(J)) {
        current <- b_mat[j, o]
        proposal <- stats::rnorm(1, mean = current, sd = sigma_b)
        
        beta_col_curr <- b_mat[, o]
        beta_col_prop <- beta_col_curr
        beta_col_prop[j] <- proposal
        
        ll_curr <- stage1_profile_loglik_one_outcome(beta_col_curr, key_col, prof_fun_list)
        ll_prop <- stage1_profile_loglik_one_outcome(beta_col_prop, key_col, prof_fun_list)
        
        lp_curr <- stats::dnorm(current, mean = mu[j], sd = sqrt(tau2[j]), log = TRUE)
        lp_prop <- stats::dnorm(proposal, mean = mu[j], sd = sqrt(tau2[j]), log = TRUE)
        
        log_acc <- (ll_prop + lp_prop) - (ll_curr + lp_curr)
        
        if (log(stats::runif(1)) < log_acc) {
          b_mat[j, o] <- proposal
          accept_b[j, o] <- accept_b[j, o] + 1L
        }
      }
    }
    
    mu_samples[iter, ] <- mu
    tau2_samples[iter, ] <- tau2
    tau_samples[iter, ] <- sqrt(tau2)
    b_samples[iter, , ] <- b_mat
  }
  
  colnames(mu_samples) <- paste0("mu_", seq_len(J))
  colnames(tau2_samples) <- paste0("tau2_", seq_len(J))
  colnames(tau_samples) <- paste0("tau_", seq_len(J))
  
  keep_idx <- (burn_in + 1L):n_iter
  
  list(
    mu_samples = mu_samples,
    tau2_samples = tau2_samples,
    tau_samples = tau_samples,
    b_samples = b_samples,
    mu_summary = summarize_matrix_samples(mu_samples[keep_idx, , drop = FALSE]),
    tau_summary = summarize_matrix_samples(tau_samples[keep_idx, , drop = FALSE]),
    b_summary = summarize_array_samples_3d(b_samples[keep_idx, , , drop = FALSE], prefix = "b"),
    accept_b = accept_b / n_iter
  )
}

# ------------------------------------------------------------
# Function: stage2_profile_loglik
# Purpose:
#   Compute the OOI profile log-likelihood under the current
#   latent factor matrices and OOI bias matrix.
# Input:
#   L_D - J x LF matrix.
#   L_O - LF x n matrix.
#   b_ooi - J x n OOI bias matrix.
#   jo_to_key_ooi - J x n key matrix.
#   prof_fun_list - named list of profile likelihood functions.
# Output:
#   Numeric scalar log-likelihood.
# ------------------------------------------------------------
stage2_profile_loglik <- function(L_D, L_O, b_ooi, jo_to_key_ooi, prof_fun_list) {
  theta_ooi <- L_D %*% L_O
  beta_ooi <- theta_ooi + b_ooi
  
  ll <- 0
  for (j in seq_len(nrow(jo_to_key_ooi))) {
    for (k in seq_len(ncol(jo_to_key_ooi))) {
      key <- jo_to_key_ooi[j, k]
      if (!key %in% names(prof_fun_list)) {
        stop(sprintf("Missing profile likelihood for key '%s'.", key))
      }
      out <- prof_fun_list[[key]](beta_ooi[j, k])
      ll <- ll + as.numeric(out$value)
    }
  }
  ll
}

# ------------------------------------------------------------
# Function: run_stage2_realdata
# Purpose:
#   Run Stage 2 MCMC on OOI profile likelihoods using latent
#   factors and bias draws generated from Stage 1 samples at
#   matching iterations.
# Input:
#   jo_to_key_ooi - J x n key matrix.
#   prof_fun_list - named list of profile likelihood functions.
#   mu_samples_stage1, tau_samples_stage1 - full Stage 1
#       samples.
#   O_total - total number of outcomes in the current analysis
#       set, including NC and OOI.
#   n_ooi - number of outcomes of interest.
#   LF - number of latent factors.
#   n_iter, burn_in - MCMC settings.
#   sigma_LD, sigma_LO - proposal SDs.
#   sd_latent_prior - prior SD for latent factors.
# Output:
#   List containing full samples, posterior summaries based on
#   post-burn-in iterations, acceptance rates, and OOI bias
#   draws used in Stage 2.
# ------------------------------------------------------------
run_stage2_realdata <- function(
    jo_to_key_ooi,
    prof_fun_list,
    mu_samples_stage1,
    tau_samples_stage1,
    O_total,
    n_ooi,
    LF = 2,
    n_iter = 5000,
    burn_in = 1000,
    sigma_LD = 0.19,
    sigma_LO = 0.33,
    sd_latent_prior = 1.0
) {
  J <- nrow(jo_to_key_ooi)
  
  if (nrow(mu_samples_stage1) < n_iter || nrow(tau_samples_stage1) < n_iter) {
    stop("Stage 1 samples must contain at least n_iter rows for Stage 2.")
  }
  
  L_D <- matrix(stats::rnorm(J * LF, mean = 0, sd = 0.1), nrow = J, ncol = LF)
  L_O <- matrix(stats::rnorm(LF * n_ooi, mean = 0, sd = 0.1), nrow = LF, ncol = n_ooi)
  
  L_D_samples <- array(NA_real_, dim = c(n_iter, J, LF))
  L_O_samples <- array(NA_real_, dim = c(n_iter, LF, n_ooi))
  theta_samples <- array(NA_real_, dim = c(n_iter, J, n_ooi))
  b_draws_ooi <- array(NA_real_, dim = c(n_iter, J, n_ooi))
  
  accept_LD <- matrix(0, nrow = J, ncol = LF)
  accept_LO <- matrix(0, nrow = LF, ncol = n_ooi)
  
  for (iter in seq_len(n_iter)) {
    b_mat <- matrix(0, nrow = J, ncol = O_total)
    for (j in seq_len(J)) {
      b_mat[j, ] <- stats::rnorm(
        O_total,
        mean = mu_samples_stage1[iter, j],
        sd = tau_samples_stage1[iter, j]
      )
    }
    
    b_ooi <- b_mat[, (O_total - n_ooi + 1L):O_total, drop = FALSE]
    
    for (j in seq_len(J)) {
      for (f in seq_len(LF)) {
        current <- L_D[j, f]
        proposal <- stats::rnorm(1, mean = current, sd = sigma_LD)
        
        ll_curr <- stage2_profile_loglik(L_D, L_O, b_ooi, jo_to_key_ooi, prof_fun_list)
        
        L_D_prop <- L_D
        L_D_prop[j, f] <- proposal
        ll_prop <- stage2_profile_loglik(L_D_prop, L_O, b_ooi, jo_to_key_ooi, prof_fun_list)
        
        lp_curr <- stats::dnorm(current, mean = 0, sd = sd_latent_prior, log = TRUE)
        lp_prop <- stats::dnorm(proposal, mean = 0, sd = sd_latent_prior, log = TRUE)
        
        log_acc <- (ll_prop + lp_prop) - (ll_curr + lp_curr)
        
        if (log(stats::runif(1)) < log_acc) {
          L_D[j, f] <- proposal
          accept_LD[j, f] <- accept_LD[j, f] + 1L
        }
      }
    }
    
    for (f in seq_len(LF)) {
      for (k in seq_len(n_ooi)) {
        current <- L_O[f, k]
        proposal <- stats::rnorm(1, mean = current, sd = sigma_LO)
        
        ll_curr <- stage2_profile_loglik(L_D, L_O, b_ooi, jo_to_key_ooi, prof_fun_list)
        
        L_O_prop <- L_O
        L_O_prop[f, k] <- proposal
        ll_prop <- stage2_profile_loglik(L_D, L_O_prop, b_ooi, jo_to_key_ooi, prof_fun_list)
        
        lp_curr <- stats::dnorm(current, mean = 0, sd = sd_latent_prior, log = TRUE)
        lp_prop <- stats::dnorm(proposal, mean = 0, sd = sd_latent_prior, log = TRUE)
        
        log_acc <- (ll_prop + lp_prop) - (ll_curr + lp_curr)
        
        if (log(stats::runif(1)) < log_acc) {
          L_O[f, k] <- proposal
          accept_LO[f, k] <- accept_LO[f, k] + 1L
        }
      }
    }
    
    theta_now <- L_D %*% L_O
    
    L_D_samples[iter, , ] <- L_D
    L_O_samples[iter, , ] <- L_O
    theta_samples[iter, , ] <- theta_now
    b_draws_ooi[iter, , ] <- b_ooi
  }
  
  keep_idx <- (burn_in + 1L):n_iter
  
  list(
    L_D_samples = L_D_samples,
    L_O_samples = L_O_samples,
    theta_samples = theta_samples,
    b_draws_ooi = b_draws_ooi,
    L_D_summary = summarize_array_samples_3d(L_D_samples[keep_idx, , , drop = FALSE], prefix = "L_D"),
    L_O_summary = summarize_array_samples_3d(L_O_samples[keep_idx, , , drop = FALSE], prefix = "L_O"),
    theta_summary = summarize_array_samples_3d(theta_samples[keep_idx, , , drop = FALSE], prefix = "theta"),
    accept_LD = accept_LD / n_iter,
    accept_LO = accept_LO / n_iter
  )
}

# ------------------------------------------------------------
# Function: run_one_realdata_period
# Purpose:
#   Run the full two-stage real-data analysis for one period.
# Input:
#   profile_dat - validated profile-level data.frame.
#   period_id - scalar analysis period.
#   exposure_ids_fixed - fixed exposure ids.
#   ooi_ids - user-specified OOI ids.
#   stage1_args - named list of Stage 1 settings.
#   stage2_args - named list of Stage 2 settings.
# Output:
#   List with Stage 1 results, Stage 2 results, and metadata.
# ------------------------------------------------------------
run_one_realdata_period <- function(
    profile_dat,
    period_id,
    exposure_ids_fixed,
    ooi_ids,
    stage1_args = list(),
    stage2_args = list()
) {
  prep <- prepare_period_data(
    profile_dat = profile_dat,
    period_id = period_id,
    exposure_ids_fixed = exposure_ids_fixed,
    ooi_ids = ooi_ids
  )
  
  stage1 <- do.call(
    run_stage1_realdata,
    c(
      list(
        jo_to_key_nc = prep$jo_to_key_nc,
        prof_fun_list = prep$prof_fun_list
      ),
      stage1_args
    )
  )
  
  stage2 <- do.call(
    run_stage2_realdata,
    c(
      list(
        jo_to_key_ooi = prep$jo_to_key_ooi,
        prof_fun_list = prep$prof_fun_list,
        mu_samples_stage1 = stage1$mu_samples,
        tau_samples_stage1 = stage1$tau_samples,
        O_total = length(prep$outcome_all),
        n_ooi = length(prep$outcome_ooi)
      ),
      stage2_args
    )
  )
  
  list(
    period_id = period_id,
    stage1 = stage1,
    stage2 = stage2,
    meta = list(
      exposure_ids = prep$exposure_ids,
      outcome_nc = prep$outcome_nc,
      outcome_ooi = prep$outcome_ooi,
      outcome_all = prep$outcome_all,
      jo_to_key_nc = prep$jo_to_key_nc,
      jo_to_key_ooi = prep$jo_to_key_ooi,
      jo_to_key_all = prep$jo_to_key_all
    )
  )
}