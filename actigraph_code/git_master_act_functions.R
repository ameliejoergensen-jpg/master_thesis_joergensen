
# My functions, adjusted for the analysis of actigraph data
# basically the same as in master functions, 
# but fit structure of DEPRESJON dataset better (for example the emissions)
## A lot of exact overlap with master functions, but ensures git_master_actigraph script is self-contained!


# Emissions are zero-inflated Gammas

# mu (K x J) > 0 - move to log scale for working scale
# sigma (K x J) > 0 - also move to log scale for working scale
# z (K x J) between 0-1 (probability of observing exact zero count)


################################################################################
################################################################################

# Get emission probabilities
# use custom zero-inflated Gamma function defined below

get_allprobs <- function(mu, sigma, z, X, T_len, J){
  "[<-" <- ADoverload("[<-")
  # X is T_len x K matrix (features in columns)
  # mu, sigma and z are K x J
  K <- ncol(X)
  allprobs <- AD(matrix(1, nrow = T_len, ncol = J))
  for (j in 1:J){
    pj <- AD(rep(1, T_len))
    for (k in 1:K){
      pj <- pj * dzigamma(X[, k], mu[k, j], sigma[k, j], z[k, j])
    }
    allprobs[, j] <- pj # independence across features
  }
  allprobs
}

################################################################################
################################################################################

# convert to unconstrained (working) scale for HMM

natural_to_unc_hmm <- function(mu, sigma, z, tpm, J){
  # mu, sigma and z = K x J
  l_div <- log(tpm / diag(tpm))
  unc_tpm <- as.vector(l_div[!diag(J)])
  c(as.vector(log(mu)), as.vector(log(sigma)), as.vector(qlogis(z)), unc_tpm)
}


################################################################################
################################################################################

# convert to natural scale for HMM

unc_to_natural_hmm <- function(parameter_vector, J, K){
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  parameter_vector <- unname(parameter_vector)
  
  n_mu  <- K * J
  unc_mu <- parameter_vector[1:n_mu]
  unc_sigma <- parameter_vector[(n_mu + 1):(2 * n_mu)]
  unc_z <- parameter_vector[(2 * n_mu + 1):(3 * n_mu)]
  unc_tpm <- parameter_vector[(3 * n_mu + 1):length(parameter_vector)]
  
  mu <- matrix(exp(unc_mu),    nrow = K, ncol = J)
  sigma <- matrix(exp(unc_sigma), nrow = K, ncol = J)
  z <- matrix(plogis(unc_z),  nrow = K, ncol = J)
  
  tpm <- diag(J)
  tpm[!tpm] <- exp(unc_tpm)
  tpm <- tpm / rowSums(tpm)
  delta <- solve(t(diag(J) - tpm + 1), rep(1, J))
  
  list(mu = mu, sigma = sigma, z = z, tpm = tpm, delta = delta)
}


################################################################################
################################################################################

# negative log likelihood for HMM

nll_hmm <- function(dat){
  function(par){
    getAll(par, dat)
    nat <- unc_to_natural_hmm(parameter_vector, J, K)
    nl_total <- 0
    for (X in data_list){
      allprobs <- get_allprobs(nat$mu, nat$sigma, nat$z, X, nrow(X), J)
      nl_total <- nl_total + my_forward_hsmm(nat$delta, nat$tpm, allprobs, nrow(allprobs))
    }
    mu <- nat$mu
    sigma <- nat$sigma
    z <- nat$z
    tpm <- nat$tpm
    delta <- nat$delta
    REPORT(mu)
    REPORT(sigma)
    REPORT(z)
    REPORT(tpm)
    REPORT(delta)
    nl_total
  }
}


################################################################################
################################################################################

# fit HMM for multiple time series

mle_hmm_multiple <- function(mu, sigma, z, tpm, J, data_list){
  K <- nrow(mu)
  pv <- natural_to_unc_hmm(mu, sigma, z, tpm, J)
  obj <- MakeADFun(nll_hmm(list(data_list = data_list, J = J, K = K)),
                   list(parameter_vector = pv), silent = TRUE)
  model <- nlminb(obj$par, obj$fn, obj$gr,
                  control = list(iter.max = 5000, eval.max = 5000))
  est <- unc_to_natural_hmm(model$par, J, K)
  list(mu = est$mu, sigma = est$sigma, z = est$z, tpm = est$tpm, delta = est$delta,
       mll = model$objective, conv = model$convergence, nit = model$iterations,
       gradient = obj$gr(model$par), estimate = model$par)
}


################################################################################
################################################################################

# with random perturbation of starting values

fit_hmm_multistart <- function(em, tpm_0, J, data_list, n_start = 4, pert_sd = 0.5){
  best_f <- NULL
  best_ll <- -Inf
  K <- nrow(em$mu0)
  for (s in seq_len(n_start)){
    mu <- if (s == 1L){
      em$mu0
    }else{
      exp(log(em$mu0) + matrix(rnorm(K * J, 0, pert_sd), nrow = K, ncol = J))}
    
    f <- mle_hmm_multiple(mu, em$sigma0, em$z0, tpm_0, J, data_list)
    if (is.null(f)) next
    ll <- -f$mll
    if (!is.finite(ll)) next
    if (ll > best_ll){ best_ll <- ll
    best_f <- f }
  }
  best_f
}


################################################################################
################################################################################


my_viterbi_decoder <- function(delta, tpm, data_list, J, mu, sigma, z,
                               model = c("HMM", "HSMM"), R_vec = NULL){
  model <- match.arg(model)
  out <- vector("list", length(data_list))
  pars_ok <- all(is.finite(delta)) && all(is.finite(tpm)) &&
    all(is.finite(mu)) && all(is.finite(sigma)) && all(sigma > 0) && all(mu > 0) &&
    all(is.finite(z)) && all(z >= 0 & z < 1)
  if (!pars_ok) return(out)
  state_map <- if (model == "HSMM") rep(seq_len(J), times = R_vec) else seq_len(J)
  
  for (i in seq_along(data_list)){
    X <- data_list[[i]]
    T_len <- nrow(X)
    allprobs <- if (model == "HSMM"){
      get_allprobs_hsmm(mu, sigma, z, X, T_len, J, R_vec)
    }else{
      get_allprobs(mu, sigma, z, X, T_len, J)}
    
    n_states <- ncol(allprobs)
    stopifnot(n_states == nrow(tpm), n_states == length(delta))
    
    xi <- matrix(0, nrow = n_states, ncol = T_len)
    backtrace <- matrix(0L, nrow = n_states, ncol = T_len)
    v <- delta * allprobs[1, ]
    s <- sum(v)
    if (s == 0){s <- 1}
    
    xi[, 1] <- v/s
    
    for (t in 2:T_len){
      M <- xi[, t-1] * tpm
      backtrace[, t] <- apply(M, 2, which.max)
      v <- allprobs[t, ] * apply(M, 2, max)
      s <- sum(v)
      if (s == 0) s <- 1
      xi[, t] <- v / s
    }
    
    best_path <- integer(T_len)
    best_path[T_len] <- which.max(xi[, T_len])
    for (t in (T_len - 1):1) best_path[t] <- backtrace[best_path[t + 1], t + 1]
    out[[i]] <- state_map[best_path]
  }
  out
}



################################################################################
################################################################################



extract_relevant_pars <- function(model_fit, model = c("hmm", "nb", "nonp", "nonp_factory"), 
                                  R_vec = NULL, J){
  
  if (is.null(model_fit)){return(NULL)}  
  
  model <- match.arg(model)
  
  if (model == "nonp"){
    
    # initial distribution
    # hsmm_delta is the compressed version for HSMM, delta the one for the expanded representation
    delta <- model_fit$delta
    indexing_states <- cumsum(c(0, R_vec))
    hsmm_delta <- sapply(1:J, function(i) sum(delta[indexing_states[i] + (1:R_vec[i]) ]))
    
    # sojourn time distribution
    d_r <- model_fit$d_r
    p_list <- model_fit$p_list
    
    # the tpms
    cond_tpm <-model_fit$cond_tpm
    hmm_tpm <- construct_tpm(J, cond_tpm, d_r, R_vec)
    
    # emission density
    mu <- model_fit$mu
    sigma <- model_fit$sigma
    z <- model_fit$z
    
    
    # final penalty strength
    lambda <- model_fit$lambda
    
    # log likelihood
    ll <- model_fit$llk
    # edf and df
    edf <- model_fit$edf
    df <- model_fit$df
    n_fixpar <- model_fit$n_fixpar # unpenalized parameters
    
    return(list(
      delta = delta,
      hsmm_delta = hsmm_delta,
      d_r = d_r,
      p_list = p_list,
      cond_tpm = cond_tpm,
      hmm_tpm = hmm_tpm,
      mu = mu, sigma = sigma, z = z,
      lambda = lambda,
      ll = ll, edf = edf, df = df, n_fixpar = n_fixpar
    ))
    
  }
  
  if (model == "nonp_factory"){
    
    delta <- model_fit$delta
    indexing_states <- cumsum(c(0, R_vec))
    hsmm_delta <- sapply(1:J, function(i) sum(delta[indexing_states[i] + (1:R_vec[i]) ]))
    
    # sojourn time distribution
    d_r <- model_fit$d_r
    p_list <- model_fit$p_list
    
    # the tpms
    cond_tpm <-model_fit$cond_tpm
    hmm_tpm <- construct_tpm(J, cond_tpm, d_r, R_vec)
    
    # emission density
    mu <- model_fit$mu
    sigma <- model_fit$sigma
    z <- model_fit$z
    
    
    # log likelihood
    ll <- -(model_fit$mll)
    
    return(list(
      delta = delta,
      hsmm_delta = hsmm_delta,
      d_r = d_r,
      p_list = p_list,
      cond_tpm = cond_tpm,
      hmm_tpm = hmm_tpm,
      mu = mu, sigma = sigma, z = z,
      ll = ll
    ))
    
  }
  
  
  if (model == "nb"){
    
    delta <-model_fit$delta
    indexing_states <- cumsum(c(0, R_vec))
    hsmm_delta <- sapply(1:J, function(i) sum(delta[indexing_states[i] + (1:R_vec[i]) ]))
    
    nb_size <- model_fit$nb_size
    nb_mu <- model_fit$nb_mu
    
    d_r <- vector("list", J)
    for (i in 1:J){
      d_r[[i]] <- dnbinom(0:(R_vec[i]-1), size = nb_size[i], mu = nb_mu[i])
    }
    
    cond_tpm <-model_fit$cond_tpm
    hmm_tpm <- construct_tpm(J, cond_tpm, d_r, R_vec)
    
    mu <- model_fit$mu
    sigma <- model_fit$sigma
    z <- model_fit$z
    
    
    ll <- -(model_fit$mll)
    
    
    return(list(
      delta = delta,
      hsmm_delta = hsmm_delta,
      d_r = d_r,
      nb_size = nb_size,
      nb_mu = nb_mu,
      cond_tpm = cond_tpm,
      hmm_tpm = hmm_tpm,
      mu = mu, sigma = sigma, z = z,
      ll = ll
    ))
    
  }
  
  if (model == "hmm"){
    return(list(delta = model_fit$delta, hmm_tpm = model_fit$tpm,
                mu = model_fit$mu, sigma = model_fit$sigma, z = model_fit$z,
                ll = -(model_fit$mll)))
  }
}



################################################################################
################################################################################

# returns negative log likelihood
# see theoretical foundations chapter of thesis (based on Zucchini, 2016)

my_forward_hsmm <- function(delta, tpm, allprobs, T_len){
  # Initialize
  phi <- delta * allprobs[1,]  # build the first forward vector 
  
  s <- sum(phi) # per step likelihood contribution and scaling factor
  # start summing the log likelihood
  ll <- log(s)
  # scale to avoid numerical underflow
  phi <- phi/s
  
  for (t in 2:T_len){
    phi <- colSums(phi * tpm) * allprobs[t,]  # = (phi %*% tpm), no %*% for AD, slow
    s <- sum(phi)
    ll <- ll + log(s)
    phi <- phi/s
  }
  
  # return negative log likelihood
  -ll
}



################################################################################
################################################################################

# Penalized negative log likelihood for the nonparametric HSMM

pen_nll <- function(par){
  
  getAll(par, dat) # unc_p, unc_mu, unc_sigma, unc_z, (unc_tpm), S, lambda, data_list, J, R_vec
  nat <- unc_to_natural_nonp_hsmm(par, J, R_vec)
  
  unc_p_list <- lapply(seq_len(J), function(i) par[[paste0("unc_p", i)]])
  pen <- penalty(unc_p_list, S, lambda)   # S = list of J matrices, lambda length J
  
  nl_total <- 0
  for (data in data_list){
    allprobs <- get_allprobs_hsmm(nat$mu, nat$sigma, nat$z, data, nrow(data), J, R_vec)
    nl_total <- nl_total + my_forward_hsmm(nat$delta, nat$tpm_hmm, allprobs, nrow(allprobs))
  }
  
  p_list <- nat$p_list
  d_r <- nat$d_r
  cond_tpm <- nat$cond_tpm
  delta <- nat$delta
  mu <- nat$mu
  sigma <- nat$sigma
  z <- nat$z
  REPORT(p_list)
  REPORT(d_r)
  REPORT(cond_tpm)
  REPORT(delta)
  REPORT(mu)
  REPORT(sigma)
  REPORT(z)
  
  nl_total + pen
}


################################################################################
################################################################################

# Construct the tpm of the state-expanded HMM approximating the HSMM

# see Pohle et al. (2022)

construct_tpm <- function(J, cond_tpm, d_r, R_vec){
  
  # AD aware operators
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  K <- sum(R_vec) # number of states in expanded state space
  
  # build vector to know where columns of each state begin
  idx <- cumsum(c(0, R_vec))
  G <- AD(matrix(0, K, K))
  
  for (i in 1:J){
    
    Ri <- R_vec[i] # max duration of unstructured start for state i
    di <- d_r[[i]] # extract dwell time probabilities
    
    cdf <- c(0, cumsum(di)[-Ri]) # holds P(dwell time <= r-1) at position r in other terms F(r-1)
    ci <- di / (1 - cdf) # the hazard function - given that the state has lasted this long, what is the probability of leaving exactly now
    cim <- 1 - ci # continuation probability
    rows_i <- idx[i] + (1:Ri) # the indices for G for sub states of i
    
    # continuation
    for (r in seq_len(Ri - 1)) {
      from <- rows_i[r]
      to <- rows_i[r + 1]
      G[from, to] <- cim[r]
    }
    
    # geometric tail
    last <- rows_i[Ri]
    G[last, last] <- cim[Ri]
    
    # exit, enter first sub state of next state aggregate
    for (j in seq_len(J)) {
      if (j == i) next
      entry_of_j <- idx[j] + 1
      G[rows_i, entry_of_j] <- cond_tpm[i, j] * ci
    }
  }
  G
}


################################################################################
################################################################################

# m-th order difference penalty matrix

# m is the order
# use R_i to denote length of vector containing the dwell time probabilities
# we fix the geometric tail at R_i + 1 as the reference category, so there are R_i free parameters

makeD <- function(R_i, m){
  D <- diag(R_i)
  for (i in seq_len(m)) {
    D <- diff(D)}
  return(D)
}


################################################################################
################################################################################

# Natural to unconstrained (working) conversion for nonparametric HSMM

# convert the natural parameters to an unconstrained scale for optimization
# skip initial state probabilities delta (inferred from tpm)
# instead, assume delta is the stationary distribution
# reduces number of parameters which need to be estimated and common practice

# use mu of var1 to identify states according to order
# mu2 remains free

natural_to_unc_nonp_hsmm <- function(mu, sigma, z, p_list, cond_tpm, J){
  
  par <- list(
    unc_mu    = log(mu),       # K x J, mu > 0
    unc_sigma = log(sigma),    # K x J
    unc_z     = qlogis(z)      # K x J, z in (0,1)
  )
  # one ragged unconstrained dwell vector per state; tail mass = reference
  for (i in seq_len(J)){
    p_i <- p_list[[i]]
    Ri1 <- length(p_i)
    par[[paste0("unc_p", i)]] <- log(p_i[-Ri1] / p_i[Ri1])   # length R_i
  }
  if (J > 2){
    tpm_r <- matrix(t(cond_tpm)[!diag(J)], J, J - 1, byrow = TRUE)
    par$unc_tpm <- as.vector(log(tpm_r / tpm_r[, 1])[, -1])
  }
  par
}


################################################################################
################################################################################

# Unconstrained to natural conversion for nonparametric HSMM

# RTMB-safe
# Nonparametric HSMM

unc_to_natural_nonp_hsmm <- function(par, J, R_vec){
  
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  
  K <- sum(R_vec)
  par <- lapply(par, unname)
  
  p_list <- lapply(1:J, function(i){
    upi <- par[[paste0("unc_p", i)]]      # length R_i
    z <- numeric(R_vec[i] + 1)
    z[1:R_vec[i]] <- upi
    exp(z) / sum(exp(z))
  })
  
  d_r <- lapply(p_list, function(p_i) p_i[-length(p_i)])
  mu <- exp(par$unc_mu)
  sigma <- exp(par$unc_sigma)
  zi <- plogis(par$unc_z)
  
  if (J > 2){
    cond_tpm <- matrix(0, J, J)
    cond_tpm[!diag(J)] <- t(matrix(c(rep(1, J), exp(par$unc_tpm)), J, J - 1))
    cond_tpm <- t(cond_tpm) / colSums(cond_tpm)
  } else {
    cond_tpm <- matrix(c(0, 1, 1, 0), 2, 2, byrow = TRUE)
  }
  
  tpm_hmm <- construct_tpm(J, cond_tpm, d_r, R_vec)
  delta <- solve(t(diag(K) - tpm_hmm + 1), rep(1, K))
  list(mu = mu, sigma = sigma, z = zi,
       p_list = p_list, d_r = d_r, cond_tpm = cond_tpm,
       tpm_hmm = tpm_hmm, delta = delta)
}


################################################################################
################################################################################

# Handle divergence issues in qreml resulting from huge lambda


make_qreml_capped <- function(cap = 1e6, exclude = 1e5) {
  g <- LaMa::qreml
  src <- deparse(body(g))
  
  src <- gsub("1e\\+?0?8", "__CAP__", src) # working-infinity cap (lambda_mapped > 1e8)
  src <- gsub("1e\\+?0?6", "__EXC__", src) # convergence-check exclusion threshold
  
  src <- sub(
    "lambda_mapped[which(lambda_mapped > __CAP__)] <- __CAP__",
    "lambda_mapped[which(lambda_mapped > __CAP__)] <- __CAP__\n    lambda_mapped[!is.finite(lambda_mapped)] <- __CAP__",
    src, fixed = TRUE)
  
  src <- gsub("__CAP__", format(cap, scientific = TRUE), src, fixed = TRUE)
  src <- gsub("__EXC__", format(exclude, scientific = TRUE), src, fixed = TRUE)
  
  # a transient NA in one lambda must not abort the convergence check
  src <- gsub("all(rel_change[convInd_unmapped] < tol)",
              "isTRUE(all(rel_change[convInd_unmapped] < tol))", src, fixed = TRUE)
  
  # guarded version for gdeterminant
  gdsafe <- c(
    "gdeterminant <- function(x, eps = NULL, log = TRUE) {",
    "    if (is.null(x)) return(NULL)",
    "    ev <- tryCatch(eigen(x)$values, error = function(e) NA_real_)",
    "    if (anyNA(ev)) return(NA_real_)",
    "    if (is.null(eps)) eps <- max(ev) * .Machine$double.eps",
    "    ld <- sum(log(ev[ev > eps]))",
    "    if (log) ld else exp(ld)",
    "}"
  )
  src <- c(src[1], gdsafe, src[-1]) # inject just after the opening "{"
  
  body(g) <- parse(text = paste(src, collapse = "\n"))[[1]]
  g # environment stays LaMa namespace, so all other internals still resolve
}

qreml_capped <- make_qreml_capped(cap = 1e6, exclude = 1e5)


################################################################################
################################################################################


get_allprobs_hsmm <- function(mu, sigma, z, X, T_len, J, R_vec){
  
  # replace assignment operator with AD-aware version
  "[<-" <- ADoverload("[<-")
  # X is T_len x K matrix (features in columns); mu, sigma, z are K x J
  K <- ncol(X)
  # build vector to know where columns of each state begin
  idx <- cumsum(c(0, R_vec))
  # pre-allocate emission probabilities
  allprobs <- AD(matrix(1, T_len, sum(R_vec)))
  
  for (i in 1:J){
    pj <- AD(rep(1, T_len))
    for (k in 1:K){
      pj <- pj * dzigamma(X[, k], mu[k, i], sigma[k, i], z[k, i])
    }
    # one emission distribution per state aggregate, recycled across its R_i sub states
    allprobs[, idx[i] + (1:R_vec[i])] <- pj
  }
  allprobs
}


################################################################################
################################################################################

# Rewrite dnbinom in an AD-safe way

d_nb <- function(ind, R_vec, nb_size, nb_mu){
  x <- 0:(R_vec[ind] - 1)
  r <- nb_size[ind]
  m <- nb_mu[ind]
  log_p <- lgamma(x + r) - lgamma(r) - lgamma(x + 1) +
    r * (log(r) - log(r + m)) + x * (log(m) - log(r + m))
  exp(log_p)
}


################################################################################
################################################################################

# Natural to unconstrained for nb-HSMM

natural_to_unc_nbhsmm <- function(mu, sigma, z, nb_size, nb_mu, cond_tpm, J){
  emis <- c(as.vector(log(mu)), as.vector(log(sigma)), as.vector(qlogis(z)))
  unc_nb_size <- log(nb_size)
  unc_nb_mu <- log(nb_mu)
  
  if (J > 2){
    tpm_r <- matrix(t(cond_tpm)[!diag(J)], J, J - 1, byrow = TRUE)
    unc_tpm <- as.vector(log(tpm_r / tpm_r[, 1])[, -1])
    parameter_vector <- c(emis, unc_nb_size, unc_nb_mu, unc_tpm)
  } else {
    parameter_vector <- c(emis, unc_nb_size, unc_nb_mu)
  }
  parameter_vector
}


################################################################################
################################################################################

# Unconstrained to natural for nb-HSMM

unc_to_natural_nbhsmm <- function(parameter_vector, J, K, R_vec){
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  parameter_vector <- unname(parameter_vector)
  
  n_mu <- K * J
  unc_mu <- parameter_vector[1:n_mu]
  unc_sigma <- parameter_vector[(n_mu + 1):(2 * n_mu)]
  unc_z <- parameter_vector[(2 * n_mu + 1):(3 * n_mu)]
  
  unc_nb_size <- parameter_vector[(3 * n_mu + 1):(3 * n_mu + J)]
  unc_nb_mu <- parameter_vector[(3 * n_mu + J + 1):(3 * n_mu + 2 * J)]
  
  mu <- matrix(exp(unc_mu),    nrow = K, ncol = J)
  sigma <- matrix(exp(unc_sigma), nrow = K, ncol = J)
  z <- matrix(plogis(unc_z),  nrow = K, ncol = J)
  
  nb_size <- exp(unc_nb_size)
  nb_mu <- exp(unc_nb_mu)
  
  if (J > 2){
    unc_tpm <- parameter_vector[(3 * n_mu + 2 * J + 1):length(parameter_vector)]
    cond_tpm <- matrix(0, J, J)
    cond_tpm[!diag(J)] <- as.vector(t(matrix(c(rep(1, J), exp(unc_tpm)), J, J - 1)))
    cond_tpm <- t(cond_tpm) / colSums(cond_tpm)
  } else {
    cond_tpm <- matrix(c(0, 1, 1, 0), 2, 2, byrow = TRUE)
  }
  
  d_r <- lapply(1:J, d_nb, R_vec = R_vec, nb_size = nb_size, nb_mu = nb_mu)
  tpm_hmm <- construct_tpm(J, cond_tpm, d_r, R_vec)
  delta <- solve(t(diag(sum(R_vec)) - tpm_hmm + 1), rep(1, sum(R_vec)))
  
  list(mu = mu, sigma = sigma, z = z,
       nb_size = nb_size, nb_mu = nb_mu,
       cond_tpm = cond_tpm, tpm_hmm = tpm_hmm, delta = delta)
}


################################################################################
################################################################################

# AD-safe negative log likelihood for nb-HSMM

nll_nbhsmm <- function(dat){
  function(par){
    getAll(par, dat)
    nat <- unc_to_natural_nbhsmm(parameter_vector, J, K, R_vec)
    nl_total <- 0
    for (X in data_list){
      allprobs <- get_allprobs_hsmm(nat$mu, nat$sigma, nat$z, X, nrow(X), J, R_vec)
      nl_total <- nl_total + my_forward_hsmm(nat$delta, nat$tpm_hmm, allprobs, nrow(allprobs))
    }
    mu <- nat$mu
    sigma <- nat$sigma
    z <- nat$z
    nb_size <- nat$nb_size
    nb_mu <- nat$nb_mu
    cond_tpm <- nat$cond_tpm
    tpm_hmm <- nat$tpm_hmm
    delta <- nat$delta
    REPORT(mu)
    REPORT(sigma)
    REPORT(z)
    REPORT(nb_size)
    REPORT(nb_mu)
    REPORT(cond_tpm)
    REPORT(tpm_hmm)
    REPORT(delta)
    nl_total
  }
}


################################################################################
################################################################################

# MLE estimate for nb-HSMM

mle_nbhsmm_multiple <- function(mu, sigma, z, nb_size, nb_mu, cond_tpm, J, data_list, R_vec){
  K <- nrow(mu)
  pv <- natural_to_unc_nbhsmm(mu, sigma, z, nb_size, nb_mu, cond_tpm, J)
  dat <- list(data_list = data_list, J = J, K = K, R_vec = R_vec)
  obj <- MakeADFun(nll_nbhsmm(dat), list(parameter_vector = pv), silent = TRUE)
  
  model <- nlminb(obj$par, obj$fn, obj$gr,
                  control = list(iter.max = 5000, eval.max = 5000))
  est <- unc_to_natural_nbhsmm(model$par, J, K, R_vec)
  
  indexing_states <- cumsum(c(0, R_vec))
  hsmm_delta <- sapply(1:J, function(i) sum(est$delta[indexing_states[i] + (1:R_vec[i])]))
  
  list(mu = est$mu, sigma = est$sigma, z = est$z,
       cond_tpm = est$cond_tpm, tpm_hmm = est$tpm_hmm,
       delta = est$delta, hsmm_delta = hsmm_delta,
       nb_size = est$nb_size, nb_mu = est$nb_mu,
       mll = model$objective, conv = model$convergence, nit = model$iterations,
       gradient = obj$gr(model$par), estimate = model$par)
}


################################################################################
################################################################################

# zero-inflated Gamma density, AD-safe
# x may contain exact zeros and NA

dzigamma <- function(x, mu, sigma, z){
  "[<-" <- ADoverload("[<-") # make AD safe
  ind <- as.numeric(!is.na(x)) # 1 = observed, 0 = missing
  x[ind == 0] <- 1 # missing contributes factor 1 later
  is0 <- as.numeric(x == 0) # 1 = exact zero
  x[is0 == 1] <- 1  # to keep dgamma away from 0
  shape <- (mu / sigma)^2
  scale <- sigma^2 / mu
  dens <- is0 * z + (1 - is0) * (1 - z) * dgamma(x, shape = shape, scale = scale)
  ind * dens + (1 - ind)
}


################################################################################
################################################################################

# Build the dwell time distribution as inferred by the nonparametric HSMM

# this accounts for the unstructured start (free parameters estimated up to R_vec i)
# followed by the geometric tail
# as explained in thesis and construction of tpm (by Pohle et al., 2022)

dwell_nonp <- function(d_r_i, R, R_eval){
  out <- numeric(R_eval)
  out[1:R] <- d_r_i # explicit part (unstructured start)
  tail_mass <- max(1 - sum(d_r_i), 0) # reference category mass beyond Rvec_i
  
  if (R_eval > R && tail_mass > 0){
    # continuation prob of the last sub-state self-loop (see construction of tpm)
    cdf_Rm1 <- sum(d_r_i[1:(R-1)]) # F(R-1)
    haz_R <- d_r_i[R] / (1 - cdf_Rm1) # ci[R]
    cont <- 1 - haz_R # cim[R], geometric continuation rate
    rr <- seq_len(R_eval - R)
    out[(R+1):R_eval] <- tail_mass * (1 - cont) * cont^(rr - 1)
  }
  return(out)
}


################################################################################
################################################################################

# Build dwell-time matrix to compare to truth (DGP fed to simulator)

# explicitly rebuild the whole matrix fed to the simulator, also the tail mass
# if R_vec was chosen too small, MISE will reflect this
# but there should be little mass beyond R_vec anyway

build_dwell_matrix <- function(pars, model, J, R_vec, R_eval){
  D <- matrix(0, J, R_eval)
  for (i in 1:J){
    
    if (model == "nonp"){
      D[i, ] <- dwell_nonp(pars$d_r[[i]], R_vec[i], R_eval)
      
    } else if (model == "nb"){
      D[i, ] <- dnbinom(0:(R_eval - 1), mu = pars$nb_mu[i], size = pars$nb_size[i])
      
    } else if (model == "hmm"){
      a_ii   <- diag(pars$hmm_tpm)[i]
      D[i, ] <- dgeom(0:(R_eval - 1), p = (1 - a_ii)) 
    }
    D[i, ] <- D[i, ]/sum(D[i, ]) # renormalize
  }
  D
}


################################################################################
################################################################################



