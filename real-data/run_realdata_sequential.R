# ============================================================
# File: real-data/run_realdata_sequential.R
# Purpose:
#   Main script for the real-data analysis in:
#   "Multi-Signal Safety Surveillance with Bayesian Latent
#   Factor Modeling and Bias Correction"
# ============================================================

rm(list = ls())

source(file.path("..", "helper", "utils_general.R"))
source(file.path("..", "helper", "utils_real_data.R"))

message_rule("Real-data sequential analysis started")

# ------------------------------------------------------------
# Section 1. User settings
# ------------------------------------------------------------

# Real-data input:
# Update 'input_file' so that it points to a local .RData file.
# The .RData file should contain one profile-level object
# already prepared for this analysis workflow.
input_file <- "PATH_TO_YOUR_PROFILE_DATA.RData"
profile_object_name <- "YOUR_PROFILE_OBJECT_NAME"

# User-specified analysis periods
# Enter the periods to be analyzed.
period_range <- c()

# User-specified OOI ids
# Enter one or more outcome-of-interest ids.
ooi_ids <- c()

# Latent factor dimension
# The default setting is provided as editable examples.
LF <- 2

# MCMC settings
# These default settings are provided as editable examples.
# Increase them if you want more Monte Carlo precision.
n_iter_stage1 <- 5000
burn_in_stage1 <- 1000
n_iter_stage2 <- 5000
burn_in_stage2 <- 1000

sigma_b <- 0.63
sigma_LD <- 0.19
sigma_LO <- 0.33

# Prior settings
sd_mu0 <- 0.20
a0 <- 2.5
b0 <- 0.06
sd_latent_prior <- 1.0

# Output directory
output_dir <- "."

# ------------------------------------------------------------
# Section 2. Input checks
# ------------------------------------------------------------

if (input_file == "PATH_TO_YOUR_PROFILE_DATA.RData") {
  stop(
    paste0(
      "Please update 'input_file <- \"...\"' before running this script.\n",
      "The .RData file should contain a profile-level object."
    )
  )
}

if (profile_object_name == "YOUR_PROFILE_OBJECT_NAME") {
  stop(
    paste0(
      "Please update 'profile_object_name <- \"...\"' before running this script.\n",
      "This should be the name of the profile-level object stored in the .RData file."
    )
  )
}

if (length(period_range) == 0L) {
  stop(
    paste0(
      "Please supply one or more analysis periods in\n",
      "    period_range <- c(...)\n",
      "before running the real-data analysis."
    )
  )
}

if (!is.numeric(period_range) || any(is.na(period_range))) {
  stop("`period_range` must be a numeric vector without missing values.")
}

period_range <- sort(unique(as.integer(period_range)))

if (length(ooi_ids) == 0L) {
  stop(
    paste0(
      "Please supply one or more outcome-of-interest ids in\n",
      "    ooi_ids <- c(...)\n",
      "before running the real-data analysis."
    )
  )
}

if (!is.numeric(ooi_ids) || any(is.na(ooi_ids))) {
  stop("`ooi_ids` must be a numeric vector without missing values.")
}

ooi_ids <- unique(as.integer(ooi_ids))

# ------------------------------------------------------------
# Section 3. Load and validate profile-level input
# ------------------------------------------------------------

message("Loading the profile-level object ...")

profile_dat <- load_rdata_object(
  input_file = input_file,
  object_name = profile_object_name
)

profile_dat <- extract_profile_data(profile_dat)

# ------------------------------------------------------------
# Section 4. Select fixed exposure-outcome pairs
# ------------------------------------------------------------

message("Selecting fixed exposure-outcome pairs over the requested period range ...")

fixed_sel <- select_fixed_pairs(
  profile_dat = profile_dat,
  period_range = period_range
)

save_rdata_named(
  file = make_output_path(output_dir, "real_fixed_pairs_used.RData"),
  named_list = list(
    period_range = period_range,
    exposure_ids_fixed = fixed_sel$exposure_ids_fixed,
    outcomes_fixed = fixed_sel$outcomes_fixed,
    fixed_pairs = fixed_sel$fixed_pairs,
    outcome_counts = fixed_sel$outcome_counts
  )
)

write.csv(
  fixed_sel$fixed_pairs,
  file = make_output_path(output_dir, "real_fixed_pairs_used.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Section 5. Run sequential analyses over periods
# ------------------------------------------------------------

realdata_run_log <- vector("list", length(period_range))

for (idx in seq_along(period_range)) {
  current_period <- period_range[idx]
  message_rule(sprintf("Running real-data analysis for period %02d", current_period))
  
  res <- run_one_realdata_period(
    profile_dat = profile_dat,
    period_id = current_period,
    exposure_ids_fixed = fixed_sel$exposure_ids_fixed,
    ooi_ids = ooi_ids,
    stage1_args = list(
      n_iter = n_iter_stage1,
      burn_in = burn_in_stage1,
      sigma_b = sigma_b,
      sd_mu0 = sd_mu0,
      a0 = a0,
      b0 = b0
    ),
    stage2_args = list(
      LF = LF,
      n_iter = n_iter_stage2,
      burn_in = burn_in_stage2,
      sigma_LD = sigma_LD,
      sigma_LO = sigma_LO,
      sd_latent_prior = sd_latent_prior
    )
  )
  
  stage1_file <- make_output_path(
    output_dir,
    sprintf("real_stage1_period_%02d.RData", current_period)
  )
  
  save_rdata_named(
    file = stage1_file,
    named_list = list(
      period_id = current_period,
      exposure_ids = res$meta$exposure_ids,
      outcome_nc = res$meta$outcome_nc,
      outcome_all = res$meta$outcome_all,
      jo_to_key_nc = res$meta$jo_to_key_nc,
      mu_samples = res$stage1$mu_samples,
      tau2_samples = res$stage1$tau2_samples,
      tau_samples = res$stage1$tau_samples,
      b_samples = res$stage1$b_samples,
      mu_summary = res$stage1$mu_summary,
      tau_summary = res$stage1$tau_summary,
      b_summary = res$stage1$b_summary,
      accept_b = res$stage1$accept_b
    )
  )
  
  stage2_file <- make_output_path(
    output_dir,
    sprintf("real_stage2_period_%02d.RData", current_period)
  )
  
  save_rdata_named(
    file = stage2_file,
    named_list = list(
      period_id = current_period,
      exposure_ids = res$meta$exposure_ids,
      outcome_ooi = res$meta$outcome_ooi,
      outcome_all = res$meta$outcome_all,
      jo_to_key_ooi = res$meta$jo_to_key_ooi,
      L_D_samples = res$stage2$L_D_samples,
      L_O_samples = res$stage2$L_O_samples,
      theta_samples = res$stage2$theta_samples,
      b_draws_ooi = res$stage2$b_draws_ooi,
      L_D_summary = res$stage2$L_D_summary,
      L_O_summary = res$stage2$L_O_summary,
      theta_summary = res$stage2$theta_summary,
      accept_LD = res$stage2$accept_LD,
      accept_LO = res$stage2$accept_LO
    )
  )
  
  realdata_run_log[[idx]] <- data.frame(
    period_id = current_period,
    J = length(res$meta$exposure_ids),
    m = length(res$meta$outcome_nc),
    n = length(res$meta$outcome_ooi),
    mean_accept_b = mean(res$stage1$accept_b),
    mean_accept_LD = mean(res$stage2$accept_LD),
    mean_accept_LO = mean(res$stage2$accept_LO)
  )
}

realdata_run_log <- do.call(rbind, realdata_run_log)
write.csv(
  realdata_run_log,
  file = make_output_path(output_dir, "realdata_run_log.csv"),
  row.names = FALSE
)

message_rule("Real-data sequential analysis completed")
message("Real-data runtime depends on profile size, number of periods, and MCMC settings.")
message(
  paste0(
    "The real data are not distributed with this repository. ",
    "Please contact the authors for questions regarding data access."
  )
)