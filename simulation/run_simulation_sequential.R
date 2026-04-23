# ============================================================
# File: simulation/run_simulation_sequential.R
# Purpose:
#   Main script for the simulation study in:
#   "Multi-Signal Safety Surveillance with Bayesian Latent
#   Factor Modeling and Bias Correction"
#
# Notes:
#   This script generates one full Day-60 dataset and then
#   performs sequential truncation analyses at
#   D_trunc = 5, 10, ..., 60.
# ============================================================

rm(list = ls())
set.seed(42)

source(file.path("..", "helper", "utils_general.R"))
source(file.path("..", "helper", "utils_simulation.R"))

message_rule("Simulation sequential analysis started")

# ------------------------------------------------------------
# Section 1. User settings
# ------------------------------------------------------------

# Simulation design
N <- 600
D <- 60
J <- 4
m <- 50
n <- 5
O <- m + n
LF <- 2
trunc_days <- seq(5, 60, by = 5)

# MCMC settings
# These default settings are provided as editable examples.
# Increase them if you want more Monte Carlo precision.
n_iter_stage1 <- 5000
burn_in_stage1 <- 1000
n_iter_stage2 <- 5000
burn_in_stage2 <- 1000

sigma_b <- 0.08
sigma_LD <- 0.10
sigma_LO <- 0.05

# Prior settings
sd_mu0 <- 0.20
a0 <- 2.5
b0 <- 0.06
sd_latent_prior <- 1.0

# Output directory
output_dir <- "."

# ------------------------------------------------------------
# Section 2. Generate one full Day-60 dataset
# ------------------------------------------------------------

message("Generating one full Day-60 exposure dataset ...")
exposure_array <- generate_exposure_array_partial_overlap(
  N = N,
  D = D,
  J = J,
  min_run_len = 3,
  max_run_len = 7,
  min_drugs_per_patient = 1,
  max_drugs_per_patient = J
)

message("Generating true parameters ...")
truth <- generate_true_parameters(
  J = J,
  m = m,
  n = n,
  LF = LF,
  sd_mu_true = 0.20,
  tau_range = c(0.10, 0.30),
  sd_latent = 0.10
)

message("Generating observed outcome counts ...")
Y_array <- generate_observed_counts(
  exposure_array = exposure_array,
  beta_mat = truth$beta_mat
)

full_flat <- flatten_full_simulation_data(
  exposure_array_full = exposure_array,
  Y_array_full = Y_array
)

save_rdata_named(
  file = make_output_path(output_dir, "sim_truth_and_full_data.RData"),
  named_list = list(
    exposure_array = exposure_array,
    Y_array = Y_array,
    sim_data_full = full_flat$sim_data,
    x_mat_full = full_flat$x_mat,
    y_mat_list_full = full_flat$y_mat_list,
    mu_true = truth$mu_true,
    tau_true = truth$tau_true,
    b_mat_true = truth$b_mat,
    L_D_true = truth$L_D_true,
    L_O_true = truth$L_O_true,
    theta_mat_true = truth$theta_mat,
    beta_mat_true = truth$beta_mat,
    ooi_idx = truth$ooi_idx
  )
)

# ------------------------------------------------------------
# Section 3. Run sequential truncation analyses
# ------------------------------------------------------------

simulation_run_log <- vector("list", length(trunc_days))

for (idx in seq_along(trunc_days)) {
  D_trunc <- trunc_days[idx]
  message_rule(sprintf("Running simulation for truncation day %02d", D_trunc))
  
  res <- run_one_simulation_trunc(
    D_trunc = D_trunc,
    exposure_array_full = exposure_array,
    Y_array_full = Y_array,
    truth = truth,
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
    sprintf("sim_stage1_trunc_%02d.RData", D_trunc)
  )
  
  save_rdata_named(
    file = stage1_file,
    named_list = list(
      D_trunc = D_trunc,
      x_mat = res$trunc_data$x_mat,
      y_mat_list = res$trunc_data$y_mat_list,
      sim_data = res$trunc_data$sim_data,
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
    sprintf("sim_stage2_trunc_%02d.RData", D_trunc)
  )
  
  save_rdata_named(
    file = stage2_file,
    named_list = list(
      D_trunc = D_trunc,
      ooi_idx = truth$ooi_idx,
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
  
  simulation_run_log[[idx]] <- data.frame(
    D_trunc = D_trunc,
    J = J,
    m = m,
    n = n,
    mean_accept_b = mean(res$stage1$accept_b),
    mean_accept_LD = mean(res$stage2$accept_LD),
    mean_accept_LO = mean(res$stage2$accept_LO)
  )
}

simulation_run_log <- do.call(rbind, simulation_run_log)
write.csv(
  simulation_run_log,
  file = make_output_path(output_dir, "simulation_run_log.csv"),
  row.names = FALSE
)

message_rule("Simulation sequential analysis completed")
message("Simulation runtime depends on MCMC settings and hardware.")