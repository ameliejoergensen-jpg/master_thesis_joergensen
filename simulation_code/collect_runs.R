################################################################################
################################################################################

# gather the per-run rds files of one scenario into a single rda for download from the HPC


################################################################################
################################################################################

base_dir <- path.expand(Sys.getenv("THESIS_DIR", unset = "~/qREML"))
scenario_id <- as.integer(Sys.getenv("SCENARIO_ID", unset = "1"))
sim_runs <- as.integer(Sys.getenv("SIM_RUNS", unset = "100"))

run_dir <- file.path(base_dir, "sim_results", sprintf("scenario_%d_runs", scenario_id))
out_file <- file.path(base_dir, "sim_results", sprintf("scenario_%d_all_runs.rda", scenario_id))

files <- list.files(run_dir, pattern = "^run_[0-9]+\\.rds$", full.names = TRUE)
files <- sort(files)

runs <- lapply(files, readRDS)
names(runs) <- sub("\\.rds$", "", basename(files))

present <- vapply(runs, function(x) x$run, integer(1))
missing_runs <- setdiff(seq_len(sim_runs), present)

# compact summary
mean_or_na <- function(x) if (is.null(x) || all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
ise_total <- function(x) if (is.null(x) || !is.list(x)) NA_real_ else x$total

model_names <- c("hmm", "nb", "nonp", "cv", "aic")

summary_df <- do.call(rbind, lapply(runs, function(r) {
  do.call(rbind, lapply(model_names, function(m) {
    fit <- r[[m]]
    data.frame(
      run = r$run,
      scenario_id = r$scenario_id,
      model = m,
      fitted = !is.null(fit$pars),
      accuracy = mean_or_na(fit$accuracy),
      local_accuracy = mean_or_na(fit$local_accuracy),
      ise = ise_total(fit$ise),
      elapsed = if (is.null(fit$elapsed)) NA_real_ else fit$elapsed,
      lambda = if (is.null(fit$pars$lambda)) NA_character_ else paste(signif(fit$pars$lambda, 4), collapse = "|"),
      stringsAsFactors = FALSE
    )
  }))
}))

rownames(summary_df) <- NULL

save(runs, summary_df, missing_runs, file = out_file, compress = "xz")

# anything missing?
if (length(missing_runs)) {
  cat("missing runs: ", paste(missing_runs, collapse = ", "), "\n", sep = "")
} else {
  cat("no missing runs\n")}





