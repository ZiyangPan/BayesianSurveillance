# ============================================================
# File: helper/utils_simulation.R
# Purpose:
#   Helper functions for the simulation of sequential
#   analysis under Bayesian latent factor modeling and bias
#   correction.
# ============================================================

# ------------------------------------------------------------
# Function: generate_exposure_array_partial_overlap
# Purpose:
#   Generate patient-by-day-by-drug exposure trajectories with
#   short exposure windows and partial overlap allowed.
# Input:
#   N - number of patients.
#   D - number of follow-up days.
#   J - number of drugs.
#   min_run_len, max_run_len - exposure window length range.
#   min_drugs_per_patient, max_drugs_per_patient - number of
#       exposed drugs per patient.
# Output:
#   3D array of dimension N x D x J.
# ------------------------------------------------------------
generate_exposure_array_partial_overlap <- function(
    N,
    D,
    J,
    min_run_len = 3,
    max_run_len = 7,
    min_drugs_per_patient = 1,
    max_drugs_per_patient = J
) {
  exposure_array <- array(0L, dim = c(N, D, J))
  
  for (i in seq_len(N)) {
    n_drugs_i <- sample(min_drugs_per_patient:max_drugs_per_patient, size = 1)
    chosen_drugs <- sort(sample(seq_len(J), size = n_drugs_i, replace = FALSE))
    
    occupied_windows <- matrix(0L, nrow = 0, ncol = 2)
    
    for (j in chosen_drugs) {
      run_len <- sample(min_run_len:max_run_len, size = 1)
      
      success <- FALSE
      tries <- 0L
      
      while (!success && tries < 200L) {
        tries <- tries + 1L
        start_day <- sample.int(D - run_len + 1L, size = 1L)
        end_day <- start_day + run_len - 1L
        
        if (nrow(occupied_windows) == 0L) {
          success <- TRUE
        } else {
          contains_existing <- any(start_day <= occupied_windows[, 1] &
                                     end_day >= occupied_windows[, 2])
          is_contained <- any(start_day >= occupied_windows[, 1] &
                                end_day <= occupied_windows[, 2])
          success <- !(contains_existing || is_contained)
        }
        
        if (success) {
          exposure_array[i, start_day:end_day, j] <- 1L
          occupied_windows <- rbind(occupied_windows, c(start_day, end_day))
        }
      }
    }
  }
  
  exposure_array
}

# ------------------------------------------------------------
# Function: generate_true_parameters
# Purpose:
#   Generate the simulation truth for the bias parameters and
#   the latent-factor-driven OOI effects.
# Input:
#   J - number of drugs.
#   m - number of negative control outcomes.
#   n - number of outcomes of interest.
#   LF - number of latent factors.
#   sd_mu_true - SD for true mu_j.
#   tau_range - range for true tau_j.
#   sd_latent - SD for latent-factor truth.
# Output:
#   List containing mu_true, tau_true, b_mat, L_D_true,
#   L_O_true, theta_mat, beta_mat, and ooi_idx.
# ------------------------------------------------------------
generate_true_parameters <- function(
    J,
    m,
    n,
    LF,
    sd_mu_true = 0.20,
    tau_range = c(0.10, 0.30),
    sd_latent = 0.10
) {
  O <- m + n
  ooi_idx <- (m + 1L):O
  
  mu_true <- stats::rnorm(J, mean = 0, sd = sd_mu_true)
  tau_true <- stats::runif(J, min = tau_range[1], max = tau_range[2])
  
  b_mat <- matrix(0, nrow = J, ncol = O)
  for (j in seq_len(J)) {
    b_mat[j, ] <- stats::rnorm(O, mean = mu_true[j], sd = tau_true[j])
  }
  
  L_D_true <- matrix(stats::rnorm(J * LF, mean = 0, sd = sd_latent), nrow = J, ncol = LF)
  L_O_true <- matrix(stats::rnorm(LF * n, mean = 0, sd = sd_latent), nrow = LF, ncol = n)
  
  theta_mat <- matrix(0, nrow = J, ncol = O)
  theta_mat[, ooi_idx] <- L_D_true %*% L_O_true
  
  beta_mat <- theta_mat + b_mat
  
  list(
    mu_true = mu_true,
    tau_true = tau_true,
    b_mat = b_mat,
    L_D_true = L_D_true,
    L_O_true = L_O_true,
    theta_mat = theta_mat,
    beta_mat = beta_mat,
    ooi_idx = ooi_idx
  )
}

# ------------------------------------------------------------
# Function: generate_observed_counts
# Purpose:
#   Generate outcome counts under the simulation.
# Input:
#   exposure_array - N x D x J array.
#   beta_mat - J x O coefficient matrix.
# Output:
#   3D array Y_array of dimension N x D x O.
# ------------------------------------------------------------
generate_observed_counts <- function(exposure_array, beta_mat) {
  dims <- dim(exposure_array)
  N <- dims[1]
  D <- dims[2]
  J <- dims[3]
  O <- ncol(beta_mat)
  
  Y_array <- array(0L, dim = c(N, D, O))
  
  for (o in seq_len(O)) {
    beta_o <- beta_mat[, o]
    
    for (i in seq_len(N)) {
      x_mat_i <- matrix(exposure_array[i, , , drop = FALSE], nrow = D, ncol = J)
      lambda_i <- exp(as.numeric(x_mat_i %*% beta_o))
      Y_array[i, , o] <- stats::rpois(D, lambda = lambda_i)
    }
  }
  
  Y_array
}

# ------------------------------------------------------------
# Function: truncate_simulation_to_matrices
# Purpose:
#   Convert truncated patient-by-day data into the same
#   flattened structure.
# Input:
#   exposure_array_full - N x D x J array.
#   Y_array_full - N x D x O array.
#   D_trunc - truncation day.
# Output:
#   List with x_mat, y_mat_list, and sim_data.
# ------------------------------------------------------------
truncate_simulation_to_matrices <- function(
    exposure_array_full,
    Y_array_full,
    D_trunc
) {
  D_full <- dim(exposure_array_full)[2]
  if (D_trunc < 1 || D_trunc > D_full) {
    stop(sprintf("D_trunc must lie between 1 and %d.", D_full))
  }
  
  N <- dim(exposure_array_full)[1]
  J <- dim(exposure_array_full)[3]
  O <- dim(Y_array_full)[3]
  
  exposure_array <- exposure_array_full[, seq_len(D_trunc), , drop = FALSE]
  Y_array <- Y_array_full[, seq_len(D_trunc), , drop = FALSE]
  
  x_mat <- matrix(0, nrow = N * D_trunc, ncol = J)
  y_mat_list <- vector("list", O)
  
  sim_data <- data.frame(
    id = rep(seq_len(N), each = D_trunc),
    day = rep(seq_len(D_trunc), times = N)
  )
  
  row_idx <- 1L
  for (i in seq_len(N)) {
    rows_i <- row_idx:(row_idx + D_trunc - 1L)
    x_mat[rows_i, ] <- matrix(exposure_array[i, , , drop = FALSE], nrow = D_trunc, ncol = J)
    row_idx <- row_idx + D_trunc
  }
  
  colnames(x_mat) <- paste0("drug", seq_len(J))
  sim_data <- cbind(sim_data, as.data.frame(x_mat))
  
  for (o in seq_len(O)) {
    y_vec <- integer(N * D_trunc)
    row_idx <- 1L
    for (i in seq_len(N)) {
      rows_i <- row_idx:(row_idx + D_trunc - 1L)
      y_vec[rows_i] <- Y_array[i, , o]
      row_idx <- row_idx + D_trunc
    }
    y_mat_list[[o]] <- y_vec
    sim_data[[paste0("Y", o)]] <- y_vec
  }
  
  list(
    D_trunc = D_trunc,
    x_mat = x_mat,
    y_mat_list = y_mat_list,
    sim_data = sim_data
  )
}

# ------------------------------------------------------------
# Function: flatten_full_simulation_data
# Purpose:
#   Convert the full Day-60 arrays into flattened objects for
#   saving and later inspection.
# Input:
#   exposure_array_full - N x D x J array.
#   Y_array_full - N x D x O array.
# Output:
#   List with x_mat_full, y_mat_list_full, and sim_data_full.
# ------------------------------------------------------------
flatten_full_simulation_data <- function(
    exposure_array_full,
    Y_array_full
) {
  truncate_simulation_to_matrices(
    exposure_array_full = exposure_array_full,
    Y_array_full = Y_array_full,
    D_trunc = dim(exposure_array_full)[2]
  )
}

# ------------------------------------------------------------
# Function: loglik_one_o
# Purpose:
#   Compute the Poisson log-likelihood contribution for one
#   outcome under beta_col.
# Input:
#   beta_col - numeric vector of length J.
#   x_mat - matrix of exposures.
#   y_o - count vector for one outcome.
# Output:
#   Numeric scalar log-likelihood.
# ------------------------------------------------------------
loglik_one_o <- function(beta_col, x_mat, y_o) {
  eta <- as.numeric(x_mat %*% beta_col)
  lambda <- exp(eta)
  sum(y_o * eta - lambda)
}

# ------------------------------------------------------------
# Function: loglik_oi_latent
# Purpose:
#   Compute the OOI log-likelihood under the latent factor
#   model.
# Input:
#   L_D - J x LF matrix.
#   L_O - LF x n matrix.
#   b_mat - J x O matrix.
#   x_mat - flattened exposure matrix.
#   y_mat_list - list of outcome vectors.
#   ooi_idx - OOI indices.
# Output:
#   Numeric scalar log-likelihood.
# ------------------------------------------------------------
loglik_oi_latent <- function(L_D, L_O, b_mat, x_mat, y_mat_list, ooi_idx) {
  theta_mat <- matrix(0, nrow = nrow(b_mat), ncol = ncol(b_mat))
  theta_mat[, ooi_idx] <- L_D %*% L_O
  beta_mat <- theta_mat + b_mat
  
  ll <- 0
  for (o in ooi_idx) {
    ll <- ll + loglik_one_o(beta_mat[, o], x_mat, y_mat_list[[o]])
  }
  ll
}

# ------------------------------------------------------------
# Function: run_stage1_simulation
# Purpose:
#   Run Stage 1 MCMC using only negative control outcomes to
#   learn mu_j, tau_j^2, and b_{j,o}.
# Input:
#   x_mat - flattened exposure matrix.
#   y_mat_list - list of outcome vectors.
#   m - number of negative control outcomes.
#   n_iter, burn_in - MCMC settings.
#   sigma_b - proposal SD for b_{j,o}.
#   sd_mu0 - prior SD for mu_j.
#   a0, b0 - inverse-gamma prior parameters for tau_j^2.
# Output:
#   List containing full samples, posterior summaries based on
#   post-burn-in iterations, and acceptance rates.
# ------------------------------------------------------------
run_stage1_simulation <- function(
    x_mat,
    y_mat_list,
    m,
    n_iter = 5000,
    burn_in = 1000,
    sigma_b = 0.08,
    sd_mu0 = 0.20,
    a0 = 2.5,
    b0 = 0.06
) {
  J <- ncol(x_mat)
  
  mu <- rep(0, J)
  tau2 <- rep(0.10, J)
  b_mat <- matrix(0, nrow = J, ncol = m)
  
  mu_samples <- matrix(NA_real_, nrow = n_iter, ncol = J)
  tau2_samples <- matrix(NA_real_, nrow = n_iter, ncol = J)
  tau_samples <- matrix(NA_real_, nrow = n_iter, ncol = J)
  b_samples <- array(NA_real_, dim = c(n_iter, J, m))
  
  accept_b <- matrix(0, nrow = J, ncol = m)
  
  for (iter in seq_len(n_iter)) {
    for (j in seq_len(J)) {
      post_var <- 1 / (m / tau2[j] + 1 / sd_mu0^2)
      post_mean <- post_var * sum(b_mat[j, ]) / tau2[j]
      mu[j] <- stats::rnorm(1, mean = post_mean, sd = sqrt(post_var))
    }
    
    for (j in seq_len(J)) {
      shape_post <- a0 + m / 2
      rate_post <- b0 + 0.5 * sum((b_mat[j, ] - mu[j])^2)
      tau2[j] <- safe_inverse_gamma_sample(1, shape = shape_post, rate = rate_post)
    }
    
    for (o in seq_len(m)) {
      y_o <- y_mat_list[[o]]
      
      for (j in seq_len(J)) {
        current <- b_mat[j, o]
        proposal <- stats::rnorm(1, mean = current, sd = sigma_b)
        
        beta_col_curr <- b_mat[, o]
        beta_col_prop <- beta_col_curr
        beta_col_prop[j] <- proposal
        
        ll_curr <- loglik_one_o(beta_col_curr, x_mat, y_o)
        ll_prop <- loglik_one_o(beta_col_prop, x_mat, y_o)
        
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
# Function: run_stage2_simulation
# Purpose:
#   Run Stage 2 MCMC on OOI outcomes using latent factors and
#   bias draws generated from Stage 1.
# Input:
#   x_mat - flattened exposure matrix.
#   y_mat_list - list of outcome vectors.
#   m - number of NC outcomes.
#   ooi_idx - indices of OOI outcomes.
#   mu_samples_stage1, tau_samples_stage1 - full Stage 1
#       samples.
#   LF - number of latent factors.
#   n_iter, burn_in - MCMC settings.
#   sigma_LD, sigma_LO - proposal SDs.
#   sd_latent_prior - prior SD for latent factors.
# Output:
#   List containing full samples, posterior summaries based on
#   post-burn-in iterations, acceptance rates, and OOI bias
#   draws used in Stage 2.
# ------------------------------------------------------------
run_stage2_simulation <- function(
    x_mat,
    y_mat_list,
    m,
    ooi_idx,
    mu_samples_stage1,
    tau_samples_stage1,
    LF = 2,
    n_iter = 5000,
    burn_in = 1000,
    sigma_LD = 0.10,
    sigma_LO = 0.05,
    sd_latent_prior = 1.0
) {
  J <- ncol(x_mat)
  O <- length(y_mat_list)
  n <- length(ooi_idx)
  
  if (nrow(mu_samples_stage1) < n_iter || nrow(tau_samples_stage1) < n_iter) {
    stop("Stage 1 samples must contain at least n_iter rows for Stage 2.")
  }
  
  L_D <- matrix(stats::rnorm(J * LF, mean = 0, sd = 0.1), nrow = J, ncol = LF)
  L_O <- matrix(stats::rnorm(LF * n, mean = 0, sd = 0.1), nrow = LF, ncol = n)
  
  L_D_samples <- array(NA_real_, dim = c(n_iter, J, LF))
  L_O_samples <- array(NA_real_, dim = c(n_iter, LF, n))
  theta_samples <- array(NA_real_, dim = c(n_iter, J, n))
  b_draws_ooi <- array(NA_real_, dim = c(n_iter, J, n))
  
  accept_LD <- matrix(0, nrow = J, ncol = LF)
  accept_LO <- matrix(0, nrow = LF, ncol = n)
  
  for (iter in seq_len(n_iter)) {
    b_mat <- matrix(0, nrow = J, ncol = O)
    for (j in seq_len(J)) {
      b_mat[j, ] <- stats::rnorm(
        O,
        mean = mu_samples_stage1[iter, j],
        sd = tau_samples_stage1[iter, j]
      )
    }
    
    for (j in seq_len(J)) {
      for (f in seq_len(LF)) {
        current <- L_D[j, f]
        proposal <- stats::rnorm(1, mean = current, sd = sigma_LD)
        
        ll_curr <- loglik_oi_latent(L_D, L_O, b_mat, x_mat, y_mat_list, ooi_idx)
        
        L_D_prop <- L_D
        L_D_prop[j, f] <- proposal
        ll_prop <- loglik_oi_latent(L_D_prop, L_O, b_mat, x_mat, y_mat_list, ooi_idx)
        
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
      for (k in seq_len(n)) {
        current <- L_O[f, k]
        proposal <- stats::rnorm(1, mean = current, sd = sigma_LO)
        
        ll_curr <- loglik_oi_latent(L_D, L_O, b_mat, x_mat, y_mat_list, ooi_idx)
        
        L_O_prop <- L_O
        L_O_prop[f, k] <- proposal
        ll_prop <- loglik_oi_latent(L_D, L_O_prop, b_mat, x_mat, y_mat_list, ooi_idx)
        
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
    b_draws_ooi[iter, , ] <- b_mat[, ooi_idx, drop = FALSE]
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
# Function: run_one_simulation_trunc
# Purpose:
#   Run the full two-stage simulation analysis for one
#   truncation day.
# Input:
#   D_trunc - truncation day.
#   exposure_array_full - full N x D x J exposure array.
#   Y_array_full - full N x D x O count array.
#   truth - list returned by generate_true_parameters().
#   stage1_args - named list of Stage 1 settings.
#   stage2_args - named list of Stage 2 settings.
# Output:
#   List with truncation data, Stage 1 results, Stage 2
#   results, and lightweight metadata.
# ------------------------------------------------------------
run_one_simulation_trunc <- function(
    D_trunc,
    exposure_array_full,
    Y_array_full,
    truth,
    stage1_args = list(),
    stage2_args = list()
) {
  dat_trunc <- truncate_simulation_to_matrices(
    exposure_array_full = exposure_array_full,
    Y_array_full = Y_array_full,
    D_trunc = D_trunc
  )
  
  m <- truth$ooi_idx[1] - 1L
  
  stage1 <- do.call(
    run_stage1_simulation,
    c(
      list(
        x_mat = dat_trunc$x_mat,
        y_mat_list = dat_trunc$y_mat_list,
        m = m
      ),
      stage1_args
    )
  )
  
  stage2 <- do.call(
    run_stage2_simulation,
    c(
      list(
        x_mat = dat_trunc$x_mat,
        y_mat_list = dat_trunc$y_mat_list,
        m = m,
        ooi_idx = truth$ooi_idx,
        mu_samples_stage1 = stage1$mu_samples,
        tau_samples_stage1 = stage1$tau_samples
      ),
      stage2_args
    )
  )
  
  list(
    D_trunc = D_trunc,
    trunc_data = dat_trunc,
    stage1 = stage1,
    stage2 = stage2,
    meta = list(
      J = ncol(dat_trunc$x_mat),
      O = length(dat_trunc$y_mat_list),
      m = m,
      n = length(truth$ooi_idx)
    )
  )
}