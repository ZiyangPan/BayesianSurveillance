# ============================================================
# File: helper/utils_general.R
# Purpose:
#   General utility functions shared by the simulation and
#   real-data analysis pipelines.
# ============================================================

# ------------------------------------------------------------
# Function: ensure_dir
# Purpose:
#   Create a directory if it does not already exist.
# Input:
#   path - character string; directory path.
# Output:
#   Invisibly returns the input path.
# ------------------------------------------------------------
ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(path)
}

# ------------------------------------------------------------
# Function: make_output_path
# Purpose:
#   Create an output directory if needed and return the full
#   file path.
# Input:
#   dir_path - character string; output directory.
#   file_name - character string; file name.
# Output:
#   Character string; full file path.
# ------------------------------------------------------------
make_output_path <- function(dir_path, file_name) {
  ensure_dir(dir_path)
  file.path(dir_path, file_name)
}

# ------------------------------------------------------------
# Function: save_rdata_named
# Purpose:
#   Save a named list of objects to an .RData file.
# Input:
#   file - character string; output .RData path.
#   named_list - named list of objects.
# Output:
#   Saves file to disk; invisibly returns the file path.
# ------------------------------------------------------------
save_rdata_named <- function(file, named_list) {
  stopifnot(is.list(named_list), length(named_list) > 0L)
  
  if (is.null(names(named_list)) || any(names(named_list) == "")) {
    stop("All elements in 'named_list' must be named.")
  }
  
  env <- new.env(parent = emptyenv())
  for (nm in names(named_list)) {
    assign(nm, named_list[[nm]], envir = env)
  }
  
  ensure_dir(dirname(file))
  save(list = names(named_list), file = file, envir = env)
  invisible(file)
}

# ------------------------------------------------------------
# Function: load_rdata_object
# Purpose:
#   Load one object from an .RData file by object name.
# Input:
#   input_file - character string; .RData file path.
#   object_name - character string; expected object name.
# Output:
#   Returns the loaded object.
# ------------------------------------------------------------
load_rdata_object <- function(input_file, object_name) {
  if (!file.exists(input_file)) {
    stop(sprintf(
      "Input file not found: %s\nPlease update 'input_file <- \"...\"' in the script.",
      input_file
    ))
  }
  
  e <- new.env(parent = emptyenv())
  loaded_names <- load(input_file, envir = e)
  
  if (!exists(object_name, envir = e, inherits = FALSE)) {
    stop(sprintf(
      paste0(
        "Object '%s' was not found in %s.\n",
        "Objects available in the file: %s\n",
        "Please make sure the .RData file contains the requested profile-level object."
      ),
      object_name,
      input_file,
      paste(loaded_names, collapse = ", ")
    ))
  }
  
  get(object_name, envir = e, inherits = FALSE)
}

# ------------------------------------------------------------
# Function: assert_required_columns
# Purpose:
#   Check whether a data.frame contains the required columns.
# Input:
#   dat - data.frame.
#   required_cols - character vector.
#   object_name - character string used in the error message.
# Output:
#   Invisibly returns TRUE if all columns are present.
# ------------------------------------------------------------
assert_required_columns <- function(dat, required_cols, object_name = "data") {
  if (!is.data.frame(dat)) {
    stop(sprintf("The %s must be a data.frame.", object_name))
  }
  
  missing_cols <- setdiff(required_cols, colnames(dat))
  if (length(missing_cols) > 0L) {
    stop(sprintf(
      "The %s is missing required columns: %s",
      object_name,
      paste(missing_cols, collapse = ", ")
    ))
  }
  
  invisible(TRUE)
}

# ------------------------------------------------------------
# Function: summarize_vector_samples
# Purpose:
#   Summarize a numeric vector of posterior samples.
# Input:
#   x - numeric vector.
#   probs - quantiles to report.
# Output:
#   Named numeric vector of summary statistics.
# ------------------------------------------------------------
summarize_vector_samples <- function(
    x,
    probs = c(0.025, 0.25, 0.5, 0.75, 0.975)
) {
  x <- as.numeric(x)
  qs <- stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE)
  
  out <- c(
    mean = mean(x, na.rm = TRUE),
    sd = stats::sd(x, na.rm = TRUE),
    qs
  )
  
  names(out) <- c(
    "mean", "sd",
    paste0("q", formatC(probs * 100, format = "f", digits = 1))
  )
  
  out
}

# ------------------------------------------------------------
# Function: summarize_matrix_samples
# Purpose:
#   Summarize posterior samples stored in a matrix where rows
#   are MCMC iterations and columns are parameters.
# Input:
#   mat - numeric matrix.
#   param_names - optional character vector of parameter names.
# Output:
#   data.frame with one row per parameter.
# ------------------------------------------------------------
summarize_matrix_samples <- function(mat, param_names = NULL) {
  if (!is.matrix(mat)) stop("'mat' must be a matrix.")
  
  if (is.null(param_names)) {
    param_names <- colnames(mat)
  }
  if (is.null(param_names)) {
    param_names <- paste0("param_", seq_len(ncol(mat)))
  }
  
  out <- lapply(seq_len(ncol(mat)), function(j) {
    s <- summarize_vector_samples(mat[, j])
    data.frame(parameter = param_names[j], t(s), row.names = NULL, check.names = FALSE)
  })
  
  do.call(rbind, out)
}

# ------------------------------------------------------------
# Function: summarize_array_samples_3d
# Purpose:
#   Summarize posterior samples stored in a 3D array with
#   dimensions iterations x rows x cols.
# Input:
#   arr - numeric 3D array.
#   prefix - character prefix used in parameter labels.
# Output:
#   data.frame with one row per matrix entry.
# ------------------------------------------------------------
summarize_array_samples_3d <- function(arr, prefix = "theta") {
  if (length(dim(arr)) != 3L) stop("'arr' must be a 3D array.")
  
  dims <- dim(arr)
  out <- vector("list", length = dims[2] * dims[3])
  idx <- 1L
  
  for (i in seq_len(dims[2])) {
    for (j in seq_len(dims[3])) {
      s <- summarize_vector_samples(arr[, i, j])
      out[[idx]] <- data.frame(
        parameter = sprintf("%s[%d,%d]", prefix, i, j),
        t(s),
        row.names = NULL,
        check.names = FALSE
      )
      idx <- idx + 1L
    }
  }
  
  do.call(rbind, out)
}

# ------------------------------------------------------------
# Function: acceptance_rate
# Purpose:
#   Compute the empirical acceptance rate from a logical or
#   numeric indicator vector.
# Input:
#   x - logical or numeric vector.
# Output:
#   Numeric scalar between 0 and 1.
# ------------------------------------------------------------
acceptance_rate <- function(x) {
  mean(as.numeric(x), na.rm = TRUE)
}

# ------------------------------------------------------------
# Function: safe_inverse_gamma_sample
# Purpose:
#   Sample from an inverse-gamma distribution using a
#   shape/rate parameterization.
# Input:
#   n - number of samples.
#   shape - positive shape parameter.
#   rate - positive rate parameter.
# Output:
#   Numeric vector of inverse-gamma draws.
# ------------------------------------------------------------
safe_inverse_gamma_sample <- function(n, shape, rate) {
  if (shape <= 0 || rate <= 0) {
    stop("shape and rate must be positive.")
  }
  1 / stats::rgamma(n = n, shape = shape, rate = rate)
}

# ------------------------------------------------------------
# Function: check_ooi_ids
# Purpose:
#   Check that all user-specified OOI ids are available in the
#   current analysis set.
# Input:
#   ooi_ids - numeric/integer vector supplied by user.
#   available_ids - numeric/integer vector of available ids.
# Output:
#   Invisibly returns TRUE if all OOI ids are available.
# ------------------------------------------------------------
check_ooi_ids <- function(ooi_ids, available_ids) {
  missing_ids <- setdiff(ooi_ids, available_ids)
  
  if (length(missing_ids) > 0L) {
    stop(sprintf(
      paste0(
        "The following user-specified OOI ids are not available in the current analysis set: %s\n",
        "Please update 'ooi_ids <- c(...)' so that all selected outcomes are present."
      ),
      paste(missing_ids, collapse = ", ")
    ))
  }
  
  invisible(TRUE)
}

# ------------------------------------------------------------
# Function: message_rule
# Purpose:
#   Print a simple formatted section header to the console.
# Input:
#   text - character string.
# Output:
#   Prints a formatted message.
# ------------------------------------------------------------
message_rule <- function(text) {
  line <- paste(rep("=", 70), collapse = "")
  message("\n", line)
  message(text)
  message(line)
}