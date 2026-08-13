
# exactly the same script as run array, with minimal adjustments to fit fifth dwell-scenario
# adjustment: skip the HMM-based warm init for nb-HSMM, because HMM struggles


base_dir <- path.expand(Sys.getenv("THESIS_DIR", unset = "~/qREML"))
scenario_id <- as.integer(Sys.getenv("SCENARIO_ID", unset = "1"))
run_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1"))
sim_runs <- as.integer(Sys.getenv("SIM_RUNS", unset = "100"))
n_select <- as.integer(Sys.getenv("N_SELECT", unset = "100"))

stopifnot(run_id >= 1, run_id <= sim_runs)

source(file.path(base_dir, "git_master_make_scenarios_complex.R"))
source(file.path(base_dir, "git_master_functions.R"))
scenarios <- read.csv(file.path(base_dir, "master_scenarios_complex.csv"))

suppressPackageStartupMessages({
  library(LaMa)
  library(RTMB)
})

################################################################################

sc <- scenarios[scenarios$scenario_id == scenario_id, ]
stopifnot(nrow(sc) == 1)

sojourn <- sojourn_distributions[[sc$sojourn]]
emis <- emission_distributions[[sc$emission]]
stopifnot(!is.null(sojourn), !is.null(emis))

N <- sc$N
T_len <- sc$T_len
order_diff <- sc$order_diff
seed <- sc$seed

J <- length(emis$mu1)
stopifnot(length(sojourn) == J,
          length(emis$sigma1) == J,
          length(emis$mu2) == J,
          length(emis$sigma2) == J)
stopifnot(all(diff(emis$mu1) > 0)) 
R_vec <- rep(20, J)

N_test <- 10 # does not cost much, for accuracy

RNGkind("L'Ecuyer-CMRG")
set.seed(seed)
cv_aic_runs <- sort(sample(seq_len(sim_runs), n_select))
do_cv_aic <- run_id %in% cv_aic_runs

set.seed(seed)
run_seeds <- vector("list", sim_runs)
.s <- .Random.seed
for (r in seq_len(sim_runs)) { 
  run_seeds[[r]] <- .s
  .s <- parallel::nextRNGStream(.s) }

run_dir <- file.path(base_dir, "sim_results", sprintf("scenario_%d_runs", scenario_id))
dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)

my_model <- list(
  y1 = list(distribution = "normal", mu = emis$mu1, sigma = emis$sigma1),
  y2 = list(distribution = "normal", mu = emis$mu2, sigma = emis$sigma2)
)

# true dwell-time probabilities
dist_d <- t(sapply(sojourn, function(d) d/sum(d)))

cond_tpm <- matrix(1/(J - 1), J, J)
diag(cond_tpm) <- 0

d_r <- lapply(1:J, function(i) dist_d[i, 1:R_vec[i]])

stopifnot(all(sapply(d_r, function(d) min(1 - cumsum(d)[-length(d)])) > 0))
stopifnot(all(sapply(d_r, sum) < 1))

check_tpm_hmm <- construct_tpm(J, cond_tpm, d_r, R_vec)
K <- sum(R_vec)
delta <- solve(t(diag(K) - check_tpm_hmm + 1), rep(1, K))
stopifnot(all(is.finite(delta)))

indexing_states <- cumsum(c(0, R_vec))
hsmm_delta <- sapply(1:J, function(i) sum(delta[indexing_states[i] + (1:R_vec[i])]))

S <- crossprod(makeD(R_vec[1], order_diff)) 

################################################################################
# global init

# sticky
tpm_0 <- matrix(0.1/(J - 1), J, J)
diag(tpm_0) <- 0.9

cond_tpm_0 <- cond_tpm

# this way all initial expectations of mean state duration match
nb_size_0 <- rep(5, J)
nb_mu_0 <- rep(10, J)

tail_floor <- 1e-4

# generic fallback dwell init (used only if a cascade stage fails)
p_list_0 <- lapply(seq_along(R_vec), function(i) {
  pmf <- dnbinom(0:(R_vec[i] - 1), mu = nb_mu_0[i], size = nb_size_0[i])
  pmf <- pmf * (1 - tail_floor)/max(sum(pmf), tail_floor)
  c(pmf, 1 - sum(pmf))
})

################################################################################
################################################################################

run_one <- function(i) {

  out_file <- file.path(run_dir, sprintf("run_%03d.rds", i))
  if (file.exists(out_file)) {
    cat(sprintf("run %03d already present, skipping\n", i))
    return(invisible(NULL))
  }
  assign(".Random.seed", run_seeds[[i]], envir = .GlobalEnv)

  # simulate data and format
  all_data <- simulate_hsmm(cond_tpm, T_len, N, hsmm_delta, my_model, dist_d)
  test_data <- simulate_hsmm(cond_tpm, T_len, N_test, hsmm_delta, my_model, dist_d)
  data_list <- lapply(1:N, function(n) {
    data.frame(
      var1 = as.vector(all_data$observations$y1[n, ]),
      var2 = as.vector(all_data$observations$y2[n, ])
    )
  })
  test_data_list <- lapply(1:N_test, function(n) {
    data.frame(
      var1 = as.vector(test_data$observations$y1[n, ]),
      var2 = as.vector(test_data$observations$y2[n, ])
    )
  })
  # N x T - true hidden states of test subjects (not used for training, used to obtain accuracy)
  hidden_states <- test_data$hidden_states

  em <- init_emissions(data_list, J) # data-driven means and shared variance

  # multistart for the HMM since most basic, keep best log likelihood
  t_hmm <- system.time(
    f_hmm <- safe(fit_hmm_multistart(em, tpm_0, J, data_list, n_start = 4, pert_sd = 0.5))
  )
  hmm_pars <- extract_relevant_pars(f_hmm, model = "hmm", J = J)

  # warm start FALSE!
  t_nb <- system.time({
    if (FALSE) {
      cond_hmm <- hmm_pars$hmm_tpm
      diag(cond_hmm) <- 0
      cond_hmm <- cond_hmm / rowSums(cond_hmm)
      a_ii <- pmin(diag(hmm_pars$hmm_tpm), 0.999)
      nb_mu_init <- a_ii / (1 - a_ii)
      nb_size_init <- rep(1, J) # geometric
      nb_mu1 <- hmm_pars$mu1
      nb_mu2 <- hmm_pars$mu2
      nb_s1 <- hmm_pars$sigma1
      nb_s2 <- hmm_pars$sigma2
      nb_cond <- cond_hmm
    } else { # fallback generic init
      nb_mu_init <- nb_mu_0
      nb_size_init <- nb_size_0
      nb_mu1 <- em$mu1_0
      nb_mu2 <- em$mu2_0
      nb_s1 <- em$sigma1_0
      nb_s2 <- em$sigma2_0
      nb_cond <- cond_tpm_0
    }
    f_nb <- safe(mle_nbhsmm_multiple(nb_mu1, nb_mu2, nb_s1, nb_s2,
                                     nb_size_init, nb_mu_init, nb_cond, J, data_list, R_vec))
  })
  nb_pars <- extract_relevant_pars(f_nb, model = "nb", R_vec = R_vec, J = J)

  # warm start for nonparametric based on NB
  t_nonp <- system.time({
    if (!is.null(nb_pars)){
      p_list_init <- lapply(seq_along(R_vec), function(k){
        pmf <- dnbinom(0:(R_vec[k] - 1), mu = nb_pars$nb_mu[k], size = nb_pars$nb_size[k])
        # rescale so the reference category keeps at least tail_floor of the mass
        pmf <- pmf * (1 - tail_floor)/max(sum(pmf), tail_floor)
        c(pmf, 1 - sum(pmf)) # length R_k + 1
      })
      np_mu1 <- nb_pars$mu1
      np_mu2 <- nb_pars$mu2
      np_s1 <- nb_pars$sigma1
      np_s2 <- nb_pars$sigma2
      np_cond <- nb_pars$cond_tpm
    } else {
      p_list_init <- p_list_0
      np_mu1 <- em$mu1_0
      np_mu2 <- em$mu2_0
      np_s1 <- em$sigma1_0
      np_s2 <- em$sigma2_0
      np_cond <- cond_tpm_0
    }
    par_0 <- natural_to_unc_nonp_hsmm(np_mu1, np_mu2, np_s1, np_s2, p_list_init, np_cond, J)

    dat <- list(
      data_list = data_list,
      J = J,
      R_vec = R_vec,
      S = S,
      lambda = rep(10, J))
    f_nonp <- safe(qreml_capped(pen_nll, par_0, dat, random = "unc_p"))

    safeguard_ran <- FALSE   
    safeguard_swapped <- FALSE
    LAMBDA_CAP <- 1e6
    LAMBDA_LO <- 1e-3
    nonp_pars <- extract_relevant_pars(f_nonp, model = "nonp", R_vec = R_vec, J = J)
    if (!is.null(nonp_pars)) {
      lam <- nonp_pars$lambda
      hit_hi <- lam >= LAMBDA_CAP * 0.999
      hit_lo <- lam <= LAMBDA_LO
      if (any(hit_hi) || any(hit_lo)) {
        safeguard_ran <- TRUE
        init_lam <- lam 
        init_lam[hit_hi] <- 0.1   
        init_lam[hit_lo] <- 1e4   
        dat_rs <- modifyList(dat, list(lambda = init_lam))
        f_rs <- safe(qreml_capped(pen_nll, par_0, dat_rs, random = "unc_p"))
        if (!is.null(f_rs) && is.finite(rll(f_rs)) &&
            (!is.finite(rll(f_nonp)) || rll(f_rs) > rll(f_nonp))) {
          f_nonp <- f_rs
          nonp_pars <- extract_relevant_pars(f_rs, model = "nonp", R_vec = R_vec, J = J)
          safeguard_swapped <- TRUE
        }  
      }   
    }     
})

  # CV selected
  t_cv_nonparametric <- system.time({
    if (do_cv_aic) {
      cv_lambda <- safe(greedy_lambda_cv(data_list, J, R_vec, S, par_0,
                                         exp_grid = 0:6, start_exp = rep(1, J), K = 5))
      cv_model <- if (is.null(cv_lambda)) NULL else {
        init_nat <- unc_to_natural_nonp_hsmm(cv_lambda$seed_par, J, R_vec)
        safe(mle_nonp_multiple(init_nat$mu1, init_nat$mu2, init_nat$sigma1, init_nat$sigma2,
                               init_nat$p_list, init_nat$cond_tpm, J, data_list, R_vec,
                               cv_lambda$lambda, S))
      }
    } else {
      cv_lambda <- NULL
      cv_model <- NULL
    }
  })

  # AIC 
  t_aic_nonparametric <- system.time({
    if (do_cv_aic) {
      aic_lambda <- safe(greedy_lambda_aic(data_list, J, R_vec, S, par_0,
                                           exp_grid = 0:6, start_exp = rep(1, J)))
      aic_model <- if (is.null(aic_lambda)) NULL else {
        init_nat <- unc_to_natural_nonp_hsmm(aic_lambda$seed_par, J, R_vec)
        safe(mle_nonp_multiple(init_nat$mu1, init_nat$mu2, init_nat$sigma1, init_nat$sigma2,
                               init_nat$p_list, init_nat$cond_tpm, J, data_list, R_vec,
                               aic_lambda$lambda, S))
      }
    } else {
      aic_lambda <- NULL
      aic_model <- NULL
    }
  })

  cv_pars <- extract_relevant_pars(cv_model, model = "nonp_factory", R_vec = R_vec, J = J)
  aic_pars <- extract_relevant_pars(aic_model, model = "nonp_factory", R_vec = R_vec, J = J)
  #nonp_pars <- extract_relevant_pars(f_nonp, model = "nonp", R_vec = R_vec, J = J)

  # emission based permutation once (to handle label switching)
  perm_hmm <- if (is.null(hmm_pars)) NA_real_ else align_states(hmm_pars$mu1, hmm_pars$mu2, emis$mu1, emis$mu2)
  perm_nb <- if (is.null(nb_pars)) NA_real_ else align_states(nb_pars$mu1, nb_pars$mu2, emis$mu1, emis$mu2)
  perm_nonp <- if (is.null(nonp_pars)) NA_real_ else align_states(nonp_pars$mu1, nonp_pars$mu2, emis$mu1, emis$mu2)
  perm_cv <- if (is.null(cv_pars)) NA_real_ else align_states(cv_pars$mu1, cv_pars$mu2, emis$mu1, emis$mu2)
  perm_aic <- if (is.null(aic_pars)) NA_real_ else align_states(aic_pars$mu1, aic_pars$mu2, emis$mu1, emis$mu2)

  # global accuracy
  acc_hmm <- if (is.null(hmm_pars)) NA_real_ else my_viterbi_acc(hmm_pars$delta, hmm_pars$hmm_tpm, test_data_list, J,
                                                                 hmm_pars$mu1, hmm_pars$mu2, hmm_pars$sigma1, hmm_pars$sigma2,
                                                                 hidden_states, model = "HMM", perm = perm_hmm)
  acc_nb <- if (is.null(nb_pars)) NA_real_ else my_viterbi_acc(nb_pars$delta, nb_pars$hmm_tpm, test_data_list, J,
                                                               nb_pars$mu1, nb_pars$mu2, nb_pars$sigma1, nb_pars$sigma2,
                                                               hidden_states, model = "HSMM", R_vec = R_vec, perm = perm_nb)
  acc_nonp <- if (is.null(nonp_pars)) NA_real_ else my_viterbi_acc(nonp_pars$delta, nonp_pars$hmm_tpm, test_data_list, J,
                                                                   nonp_pars$mu1, nonp_pars$mu2, nonp_pars$sigma1, nonp_pars$sigma2,
                                                                   hidden_states, model = "HSMM", R_vec = R_vec, perm = perm_nonp)
  acc_cv <- if (is.null(cv_pars)) NA_real_ else my_viterbi_acc(cv_pars$delta, cv_pars$hmm_tpm, test_data_list, J,
                                                               cv_pars$mu1, cv_pars$mu2, cv_pars$sigma1, cv_pars$sigma2,
                                                               hidden_states, model = "HSMM", R_vec = R_vec, perm = perm_cv)
  acc_aic <- if (is.null(aic_pars)) NA_real_ else my_viterbi_acc(aic_pars$delta, aic_pars$hmm_tpm, test_data_list, J,
                                                                 aic_pars$mu1, aic_pars$mu2, aic_pars$sigma1, aic_pars$sigma2,
                                                                 hidden_states, model = "HSMM", R_vec = R_vec, perm = perm_aic)

  # local accuracy
  lacc_hmm <- if (is.null(hmm_pars)) NA_real_ else my_local_acc(hmm_pars$delta, hmm_pars$hmm_tpm, test_data_list, J,
                                                                hmm_pars$mu1, hmm_pars$mu2, hmm_pars$sigma1, hmm_pars$sigma2,
                                                                hidden_states, model = "HMM", perm = perm_hmm)
  lacc_nb <- if (is.null(nb_pars)) NA_real_ else my_local_acc(nb_pars$delta, nb_pars$hmm_tpm, test_data_list, J,
                                                              nb_pars$mu1, nb_pars$mu2, nb_pars$sigma1, nb_pars$sigma2,
                                                              hidden_states, model = "HSMM", R_vec = R_vec, perm = perm_nb)
  lacc_nonp <- if (is.null(nonp_pars)) NA_real_ else my_local_acc(nonp_pars$delta, nonp_pars$hmm_tpm, test_data_list, J,
                                                                  nonp_pars$mu1, nonp_pars$mu2, nonp_pars$sigma1, nonp_pars$sigma2,
                                                                  hidden_states, model = "HSMM", R_vec = R_vec, perm = perm_nonp)
  lacc_cv <- if (is.null(cv_pars)) NA_real_ else my_local_acc(cv_pars$delta, cv_pars$hmm_tpm, test_data_list, J,
                                                              cv_pars$mu1, cv_pars$mu2, cv_pars$sigma1, cv_pars$sigma2,
                                                              hidden_states, model = "HSMM", R_vec = R_vec, perm = perm_cv)
  lacc_aic <- if (is.null(aic_pars)) NA_real_ else my_local_acc(aic_pars$delta, aic_pars$hmm_tpm, test_data_list, J,
                                                                aic_pars$mu1, aic_pars$mu2, aic_pars$sigma1, aic_pars$sigma2,
                                                                hidden_states, model = "HSMM", R_vec = R_vec, perm = perm_aic)

  # ISE
  ise_hmm <- if (is.null(hmm_pars)) NA_real_ else dwell_mise(dist_d, hmm_pars, model = "hmm", J, R_vec, perm_hmm)
  ise_nb <- if (is.null(nb_pars)) NA_real_ else dwell_mise(dist_d, nb_pars, model = "nb", J, R_vec, perm_nb)
  ise_nonp <- if (is.null(nonp_pars)) NA_real_ else dwell_mise(dist_d, nonp_pars, model = "nonp", J, R_vec, perm_nonp)
  ise_cv <- if (is.null(cv_pars)) NA_real_ else dwell_mise(dist_d, cv_pars, model = "nonp", J, R_vec, perm_cv)
  ise_aic <- if (is.null(aic_pars)) NA_real_ else dwell_mise(dist_d, aic_pars, model = "nonp", J, R_vec, perm_aic)

  result <- list(
    run = i,
    scenario_id = scenario_id,
    scenario = sc,
    J = J,
    R_vec = R_vec,
    cv_aic_selected = do_cv_aic,
    em_0 = em,
    hmm = list(pars = hmm_pars, elapsed = unname(t_hmm["elapsed"]), accuracy = acc_hmm, local_accuracy = lacc_hmm, ise = ise_hmm),
    nb = list(pars = nb_pars, elapsed = unname(t_nb["elapsed"]), accuracy = acc_nb, local_accuracy = lacc_nb, ise = ise_nb),
    nonp = list(pars = nonp_pars, elapsed = unname(t_nonp["elapsed"]), accuracy = acc_nonp, local_accuracy = lacc_nonp, ise = ise_nonp, safeguard_ran = safeguard_ran, safeguard_swapped = safeguard_swapped),
    cv = list(pars = cv_pars, elapsed = unname(t_cv_nonparametric["elapsed"]), accuracy = acc_cv, local_accuracy = lacc_cv,
              ise = ise_cv, lambda = cv_lambda$lambda),
    aic = list(pars = aic_pars, elapsed = unname(t_aic_nonparametric["elapsed"]), accuracy = acc_aic, local_accuracy = lacc_aic,
               ise = ise_aic, lambda = aic_lambda$lambda)
  )

  rm(all_data, test_data, data_list, test_data_list, dat, f_hmm, f_nb, f_nonp, aic_model, cv_model)
  gc()

  # atomic write, so a killed task never leaves a half-written rds
  tmp <- paste0(out_file, ".tmp")
  saveRDS(result, tmp, compress = FALSE)
  file.rename(tmp, out_file)

  cat(sprintf("%s  scenario %d run %03d done (cv_aic = %s)\n",
              format(Sys.time()), scenario_id, i, do_cv_aic))
  invisible(NULL)
}

################################################################################

cat(sprintf("scenario %d | run %d of %d | J = %d | R = %d | cv_aic = %s\n",
            scenario_id, run_id, sim_runs, J, R_vec[1], do_cv_aic))

run_one(run_id)



