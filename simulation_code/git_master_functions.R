# Amelie Jörgensen
# Master thesis
# Automatic smoothness selection for nonparametric dwell-time distributions in hidden semi-Markov models

################################################################################
################################################################################

# Notations:

# N: number of subjects (number of independent sequences)
# T_len: length of the time series
# J: number of hidden states

# delta: the stationary distribution of the Markov chain, here also assumed to be the initial distribution
# cond_tpm: the conditional transition probability matrix of the HSMM
# tpm: transition probability matrix of the HMM
# var_model: specifies the emission distribution
# dist_d: contains the dwell time probabilities (true DGP)
# d_r: contains the dwell time probabilities but stored in a list to handle different R_i (estimated)
# p_list: basically the same as d_r, but R_i + 1, where R_i + 1 contains the probability mass of the geometric tail
# R_vec: a vector of length J, specifying the state aggregate size for the state expansion trick


################################################################################
################################################################################


# 1) Function to simulate data

# var_model needs to be a list of the observed variables, 
# specifying their distribution and parameters

# dist_d contains the dwell time probabilities
# for a HMM, set dist_d to a geometric distribution

simulate_hsmm <- function(cond_tpm, T_len, N, delta, var_model, dist_d){
  # hidden number of states
  J <- nrow(cond_tpm)
  # max duration
  D_max <- ncol(dist_d)
  # state matrix for all subjects
  states <- matrix(NA, nrow = N, ncol = T_len)
  durations <- vector("list", N)
  # sample initial state according to init prob given by pi
  # we allow right-censoring but assume the first state starts at t = 1 (no left-censoring, standard) for simulation
  states[, 1] <- sample.int(J, size = N, replace = TRUE, prob = delta)
  # sample duration of first state
  # note: the simulator assumes the sequence starts at t = 1, even though this is not strictly necessary for the models
  for (n in 1:N){
    durations[[n]] <- sample.int(D_max, size = 1, prob = dist_d[states[n, 1], ])
    states[n, 1:durations[[n]]] <- states[n,1]
  }
  
  for (n in 1:N){
    while (length(states[n, !is.na(states[n, ])]) < T_len){
      
      length_existing_states <- length(states[n, !is.na(states[n, ])])
      next_state <- sample(J, size = 1, prob = cond_tpm[tail(states[n, !is.na(states[n, ])], 1), ])
      next_duration <- sample.int(D_max, size = 1, prob = dist_d[next_state, ])
      durations[[n]] <- c(durations[[n]], next_duration)
      
      if ((length_existing_states + next_duration) <= T_len ){
        
        states[n, (length_existing_states + 1):(length_existing_states + next_duration)] <- next_state
        
      }else{
        
        states[n, (length_existing_states + 1):T_len] <- next_state
      }
    }
  }
  
  # generate observations
  flat_states <- as.vector(states)
  obs <- lapply(var_model, function(var){
    
    y <- switch(var$distribution,
                "normal" = rnorm(n = length(flat_states), mean = var$mu[flat_states], sd = var$sigma[flat_states]),
                "poisson" = rpois(n = length(flat_states), lambda = var$lambda[flat_states]),
                "bernoulli" = rbinom(n = length(flat_states), size = 1, prob = var$p[flat_states])
    )
    
    matrix(y, nrow = N, ncol = T_len)
  })
  
  simulated_hsmm <- list("hidden_states" = states, observations = obs, duration = durations)
  
  return(simulated_hsmm)
  
}


################################################################################
################################################################################

# 2) Reshape observations from simulator

# used before for previous HMM implementation (EM), but maybe still useful for log emissions or something else
# reshape simulated data into list format
# all_vars should be a list of matrices, one for each variable of shape N x T obtained from the current simulation generator
# output is a nested list [subject][var]
# simulation generates same T for now, but if we want variable T just change all_vars to list of lists

reshape_observations <- function(all_vars){
  
  N <- nrow(all_vars[[1]]) # all subjects observed across all variables
  
  all_observations <- lapply(1:N, function(i) {
    lapply(all_vars, function(var) var[i, ])
    
  })
  
  return(all_observations)
  
}


################################################################################
################################################################################

# 3) Compute log emissions

# compute log emissions
# user input:
# observations stored in a list for one subject, so from output above: reshaped observations [subject]
# var_model specifying distribution and parameters for observed vars

# for now: normal, poisson, bernoulli

# Output: J x T vector for joint log emission probabilities of all vars per state per time point

get_log_emissions <- function(observations, var_model, J){
  
  number_of_variables <- length(var_model)
  
  T <- length(observations[[1]]) # all variables observed across all time points of that subject
  
  # we need this per variable, then add them together
  log_emissions_one_var <- list()
  
  for (var in names(var_model)){
    
    log_emissions_one_var[[var]] <- matrix(NA_real_, J, T)
    
    dist <- var_model[[var]]$distribution
    y <- observations[[var]]
    
    if(!dist %in% c("normal", "poisson", "bernoulli")) stop("Distribution must be normal, poisson, or bernoulli")
    
    if (dist == "normal"){
      
      for (i in 1:J){
        # should both be of length J
        mu <- var_model[[var]]$mu
        sigma <- var_model[[var]]$sigma 
        
        log_emissions_one_var[[var]][i, ] <- dnorm(x = y, mean = mu[i], sd = sigma[i], log = TRUE)
      }
    }
    
    if (dist == "poisson"){
      
      for (i in 1:J){
        
        lambda <- var_model[[var]]$lambda
        
        log_emissions_one_var[[var]][i, ] <- dpois(x = y, lambda = lambda[i], log = TRUE)
      }
    }
    
    if (dist == "bernoulli"){
      
      for (i in 1:J){
        
        p <- var_model[[var]]$p
        
        log_emissions_one_var[[var]][i, ] <- dbinom(x = y, size = 1, prob = p[i], log = TRUE)
      }
    }
    
  }
  
  # for all vars - assumes independence
  log_emissions <- Reduce(`+`, log_emissions_one_var) # should be J x T
  
  return(log_emissions)
  
}


################################################################################
################################################################################

# 4) Compute emissions

# simpler, efficient version of get log emissions
# on probability scale
# assume two observed variables, both normally distributed
# can add different distributions later
# furthermore, we assume that mu1 < mu2 < mu3 for states 1, ..., 3 on one variable to avoid label switching
# but technically not necessary anymore since now implemented emission based permutation after model fit

get_allprobs <- function(mu1, mu2, sigma1, sigma2, data, T_len, J){
  # if an observation is missing, it stays at 1 and essentially contributes nothing
  allprobs <- AD(matrix(1, nrow = T_len, ncol = J))
  var1_index <- as.numeric(!is.na(data$var1)) # non-missing observations
  var2_index <- as.numeric(!is.na(data$var2))
  
  var1 <- data$var1
  var1[var1_index == 0] <- 0 # to handle missing values later
  var2 <- data$var2
  var2[var2_index == 0] <- 0
  
  for (i in 1:J){
    # missing values get assigned 1 and contribute nothing
    var1_probs <- var1_index * dnorm(x = var1, mean = mu1[i], sd = sigma1[i]) + (1 - var1_index)
    var2_probs <- var2_index * dnorm(x = var2, mean = mu2[i], sd = sigma2[i]) + (1 - var2_index)
    # assumes independence
    allprobs[, i] <- var1_probs * var2_probs
  }
  
  return(allprobs)
  
}


################################################################################
################################################################################

# 5) Compute emissions - expanded HMM

# same as the above, but accounting for the state aggregates
# RTMB safe

get_allprobs_hsmm <- function(mu1, mu2, sigma1, sigma2, data, T_len, J, R_vec){
  
  # replace assignment operator with AD-aware version
  "[<-" <- ADoverload("[<-")
  # build vector to know where columns of each state begin
  idx <- cumsum(c(0, R_vec))
  # pre-allocate emission probabilities
  allprobs <- AD(matrix(1, T_len, sum(R_vec)))
  
  var1_index <- as.numeric(!is.na(data$var1)) # non-missing observations
  var2_index <- as.numeric(!is.na(data$var2))
  
  var1 <- data$var1
  var1[var1_index == 0] <- 0
  var2 <- data$var2
  var2[var2_index == 0] <- 0
  
  for (i in 1:J){
    # missing values get assigned 1 and contribute nothing
    var1_probs <- var1_index * dnorm(x = var1, mean = mu1[i], sd = sigma1[i]) + (1 - var1_index)
    var2_probs <- var2_index * dnorm(x = var2, mean = mu2[i], sd = sigma2[i]) + (1 - var2_index)
    # one emission distribution per state, recycled across its R_i sub states
    # indexing ensures it is matched to state aggregate, recycled column-wise
    allprobs[, idx[i] + (1:R_vec[i])] <- var1_probs * var2_probs
  }
  allprobs
}

################################################################################
################################################################################

# 6) The Viterbi algorithm for hidden Markov models

# Global state decoding

my_viterbi <- function(delta, tpm, observations, var_model){
  
  T_len <- length(observations[[1]])
  J <- nrow(tpm)
  
  log_v <- matrix(NA, nrow = J, ncol = T_len)
  backtrace <- matrix(NA, nrow = J, ncol = T_len)
  
  log_emit <- get_log_emissions(observations = observations, var_model = var_model, J = J)
  log_tpm <- log(tpm)
  
  # Initialize
  log_v[, 1] <- log(delta) + log_emit[, 1]
  
  # Recursion
  for (t in 2:T_len){
    
    # J x J
    vals <- log_v[, t-1] + log_tpm
    
    # max over j for each i
    log_v[, t] <- log_emit[, t] + apply(vals, 2, max)
    # keep track of previous best state
    backtrace[, t] <- apply(vals, 2, which.max)
    
  }
  
  # back-trace
  best_path <- numeric(T_len)
  best_path[T_len] <- which.max(log_v[, T_len])
  
  for (t in (T_len-1):1){
    best_path[t] <- backtrace[best_path[t+1], t+1]
  }
  
  return(list(best_path = best_path))
  
}


################################################################################
################################################################################

# 7) Local decoder for hidden Markov models

my_local_acc <- function(delta, tpm, data_list, J,
                         mu1, mu2, sigma1, sigma2,
                         hidden_states, 
                         model = c("HMM", "HSMM"), R_vec = NULL, perm){
  
  model <- match.arg(model)
  stopifnot(length(data_list) == nrow(hidden_states))
  
  # to match expanded state space
  if (model == "HSMM"){
    state_map <- rep(seq_len(J), times = R_vec)
  }else{
    state_map <- seq_len(J)
  }
  
  accuracy_per_subject <- numeric(length(data_list))
  
  for (i in seq_along(data_list)){
    
    data <- data_list[[i]]
    true_states <- hidden_states[i, ]
    T_len <- nrow(data)
    
    if (model == "HSMM"){
      allprobs <- get_allprobs_hsmm(mu1, mu2, sigma1, sigma2, data, T_len, J, R_vec)
    }else{
      allprobs <- get_allprobs(mu1, mu2, sigma1, sigma2, data, T_len, J)
    }
    
    n_states <- ncol(allprobs)
    stopifnot(n_states == nrow(tpm), n_states == length(delta))
    
    # scaled forward pass
    alpha <- matrix(0, nrow = T_len, ncol = n_states)
    cscale <- numeric(T_len)
    
    phi <- delta * allprobs[1,]  
    cscale[1] <- sum(phi) 
    
    alpha[1, ] <- phi/cscale[1]
    
    for (t in 2:T_len){
      phi <- colSums(alpha[t-1, ] * tpm) * allprobs[t,]
      cscale[t] <- sum(phi)
      alpha[t, ] <- phi/cscale[t]
    }
    
    # scaled backward pass
    beta <- matrix(0, nrow = T_len, ncol = n_states)
    beta[T_len, ] <- 1
    
    for (t in (T_len - 1):1){
      b <- colSums(t(tpm) * (beta[t+1, ] * allprobs[t+1, ]))
      beta[t, ] <- b/cscale[t+1]
    }
    
    gamma <- alpha * beta
    gamma <- gamma/rowSums(gamma)
    
    # for HSMM, take sum within each state block before choosing max
    if (model == "HSMM"){
      # T x J
      gamma_state <- t(rowsum(t(gamma), group = state_map))   
    }else{
      gamma_state <- gamma
    }
    
    decoded_fitted <- max.col(gamma_state, ties.method = "first")
    decoded_true <- perm[decoded_fitted]
    accuracy_per_subject[i] <- mean(decoded_true == true_states)
  }
  
  return(accuracy_per_subject)
  
}  




################################################################################
################################################################################

# 8) Forward algorithm

# in case of HSMM, delta refers to the expanded state space and has length of the sum of all R_i
# the same applies for tpm, it is the tpm of the expanded HMM

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

# 9) Construct the tpm of the state-expanded HMM approximating the HSMM

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

# 9) m-th order difference penalty matrix

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

# 9) Natural to unconstrained (working) conversion for nonparametric HSMM

# convert the natural parameters to an unconstrained scale for optimization
# skip initial state probabilities delta (inferred from tpm)
# instead, assume delta is the stationary distribution
# reduces number of parameters which need to be estimated and common practice

# use mu of var1 to identify states according to order
# mu2 remains free

# for now needs equal R_i across states

natural_to_unc_nonp_hsmm <- function(mu1, mu2, sigma1, sigma2, p_list, cond_tpm, J){
  
  unc_p <- t(sapply(p_list, function(p_i){
    Ri1 <- length(p_i) # Ri1 = R_i + 1
    log(p_i[-Ri1] / p_i[Ri1]) # tail mass is the reference category   
    # so p_i is exactly length R_i
  }))
  
  if (!is.matrix(unc_p)){
    unc_p <- matrix(unc_p, nrow = J)}
  
  par <- list(
    unc_p = unc_p, # J x R matrix
    # take log of positive differences, since mu1 < mu2 < mu3 for variable 1 for identification
    unc_mu1 = c(mu1[1], log(diff(mu1))),
    unc_mu2 = mu2,
    unc_sigma1 = log(sigma1),
    unc_sigma2 = log(sigma2)
  )
  if (J > 2){
    tpm_r <- matrix(t(cond_tpm)[!diag(J)], J, J - 1, byrow = TRUE) # get rid of empty diagonal
    par$unc_tpm <- as.vector(log(tpm_r / tpm_r[, 1])[, -1]) # fix first off-diagonal entry as reference category
  }
  par
}


################################################################################
################################################################################

# 10) Unconstrained to natural conversion for nonparametric HSMM

# RTMB-safe
# Nonparametric HSMM

unc_to_natural_nonp_hsmm <- function(par, J, R_vec){
  
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  
  K <- sum(R_vec)
  par <- lapply(par, unname)
  
  p_list <- lapply(1:J, function(i){
    z <- numeric(R_vec[i] + 1)
    z[1:R_vec[i]] <- par$unc_p[i, ]
    exp(z) / sum(exp(z))
  })
  
  d_r <- lapply(p_list, function(p_i) p_i[-length(p_i)])
  mu1 <- par$unc_mu1[1] + c(0, cumsum(exp(par$unc_mu1[-1])))
  mu2 <- par$unc_mu2
  sigma1 <- exp(par$unc_sigma1)
  sigma2 <- exp(par$unc_sigma2)
  
  if (J > 2){
    cond_tpm <- matrix(0, J, J)
    cond_tpm[!diag(J)] <- t(matrix(c(rep(1, J), exp(par$unc_tpm)), J, J - 1))
    cond_tpm <- t(cond_tpm) / colSums(cond_tpm)
  } else {
    cond_tpm <- matrix(c(0, 1, 1, 0), 2, 2, byrow = TRUE)
  }
  
  tpm_hmm <- construct_tpm(J, cond_tpm, d_r, R_vec)
  # solve for stationary distribution
  delta <- solve(t(diag(K) - tpm_hmm + 1), rep(1, K))
  list(mu1 = mu1, mu2 = mu2, sigma1 = sigma1, sigma2 = sigma2,
       p_list = p_list, d_r = d_r, cond_tpm = cond_tpm, tpm_hmm = tpm_hmm, delta = delta)
}




################################################################################
################################################################################

# 11) Penalized negative log likelihood for the nonparametric HSMM

# RTMB-safe
# handles multiple observation sequences stored in data_list
# idx[i] is always end of last state aggregate block

pen_nll <- function(par){
  
  getAll(par, dat) # exposes unc_p (each vector length R_i), S, lambda, data_list, J, R_vec
  nat <- unc_to_natural_nonp_hsmm(par, J, R_vec)
  
  pen <- penalty(unc_p, S, lambda) # penalty term (from LaMa, works with qreml function)
  
  nl_total <- 0
  
  for (data in data_list){
    allprobs <- get_allprobs_hsmm(nat$mu1, nat$mu2, nat$sigma1, nat$sigma2, data, nrow(data), J, R_vec)
    nl_total <- nl_total + my_forward_hsmm(nat$delta, nat$tpm_hmm, allprobs, nrow(allprobs))
  }
  
  # assign names
  p_list <- nat$p_list
  d_r <- nat$d_r
  cond_tpm <- nat$cond_tpm
  delta <- nat$delta
  mu1 <- nat$mu1
  mu2 <- nat$mu2
  sigma1 <- nat$sigma1
  sigma2 <- nat$sigma2
  REPORT(p_list) 
  REPORT(d_r)
  REPORT(cond_tpm)
  REPORT(delta)
  REPORT(mu1)
  REPORT(mu2)
  REPORT(sigma1)
  REPORT(sigma2)
  
  nl_total + pen
}


################################################################################
################################################################################

# 12) Natural to unconstrained for HMM

natural_to_unc_hmm <- function(mu1, mu2, sigma1, sigma2, tpm, J){
  # transform transition probability matrix
  l_div <- log(tpm/diag(tpm))
  unc_tpm <- as.vector(l_div[!diag(J)]) # keep the diagonal fixed
  
  # transform emission parameters
  # take log of positive differences, since mu1 < mu2 < mu3 for variable 1 for identification
  unc_mu1 <- c(mu1[1], log(diff(mu1)))
  unc_mu2 <- mu2
  
  unc_sigma1 <- log(sigma1)
  unc_sigma2 <- log(sigma2)
  
  parameter_vector <- c(unc_mu1, unc_mu2, unc_sigma1, unc_sigma2, unc_tpm)
  return(parameter_vector)
}



################################################################################
################################################################################

# 13) Unconstrained to natural HMM

unc_to_natural_hmm <- function(parameter_vector, J){
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  
  parameter_vector <- unname(parameter_vector)
  
  unc_mu1 <- parameter_vector[1:J]
  unc_mu2 <- parameter_vector[(J+1):(2*J)]
  
  unc_sigma1 <- parameter_vector[(2*J + 1): (3*J)]
  unc_sigma2 <- parameter_vector[(3*J + 1): (4*J)]
  
  unc_tpm <- parameter_vector[(4*J + 1): length(parameter_vector)]
  
  mu1 <- unc_mu1[1] + c(0, cumsum(exp(unc_mu1[-1])))
  mu2 <- unc_mu2
  
  sigma1 <- exp(unc_sigma1)
  sigma2 <- exp(unc_sigma2)
  
  tpm <- diag(J)
  tpm[!tpm] <- exp(unc_tpm)
  tpm <- tpm/rowSums(tpm)
  
  delta <- solve(t(diag(J)-tpm+1),rep(1,J))
  
  return(list(mu1 = mu1, mu2 = mu2, sigma1 = sigma1, sigma2 = sigma2, tpm = tpm, delta = delta))
}


################################################################################
################################################################################

# 14) Natural to unconstrained for nb-HSMM

natural_to_unc_nbhsmm <- function(mu1, mu2, sigma1, sigma2, nb_size, nb_mu, cond_tpm, J){
  
  # transform emission parameters
  # same as before
  unc_mu1 <- c(mu1[1], log(diff(mu1)))
  unc_mu2 <- mu2
  
  unc_sigma1 <- log(sigma1)
  unc_sigma2 <- log(sigma2)
  
  # for nb sojourn time
  unc_nb_size <- log(nb_size)
  unc_nb_mu <- log(nb_mu)
  
  if (J > 2){
    tpm_r <- matrix(t(cond_tpm)[!diag(J)], J, J-1, byrow=TRUE) # discard the 0 diagonal
    unc_tpm <- as.vector(log(tpm_r/tpm_r[,1])[,-1])  # use the first off-diagonal entry as reference category
    parameter_vector <- c(unc_mu1, unc_mu2, unc_sigma1, unc_sigma2, unc_nb_size, unc_nb_mu, unc_tpm)
  }else{
    parameter_vector <- c(unc_mu1, unc_mu2, unc_sigma1, unc_sigma2, unc_nb_size, unc_nb_mu)
  }
  
  return(parameter_vector)
}


################################################################################
################################################################################

# 15) Unconstrained to natural for nb-HSMM

unc_to_natural_nbhsmm <- function(parameter_vector, J, R_vec){
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  
  parameter_vector <- unname(parameter_vector)
  
  unc_mu1 <- parameter_vector[1:J]
  unc_mu2 <- parameter_vector[(J+1):(2*J)]
  
  unc_sigma1 <- parameter_vector[(2*J + 1): (3*J)]
  unc_sigma2 <- parameter_vector[(3*J + 1): (4*J)]
  
  unc_nb_size <- parameter_vector[(4*J +1) : (5*J)]
  unc_nb_mu <- parameter_vector[(5*J +1) : (6*J)]
  
  mu1 <- unc_mu1[1] + c(0, cumsum(exp(unc_mu1[-1])))
  mu2 <- unc_mu2
  
  sigma1 <- exp(unc_sigma1)
  sigma2 <- exp(unc_sigma2)
  
  nb_size <- exp(unc_nb_size)
  nb_mu <- exp(unc_nb_mu)
  
  if (J > 2){
    unc_tpm <- parameter_vector[(6*J +1): length(parameter_vector)]
    cond_tpm <- matrix(0, J, J)
    cond_tpm[!diag(J)] <- as.vector(t(matrix(c(rep(1,J), exp(unc_tpm)), J, J-1)))
    cond_tpm <- t(cond_tpm)/colSums(cond_tpm)
  }else{
    cond_tpm <- matrix(c(0, 1, 1, 0), 2, 2, byrow = TRUE)
  }
  
  d_r <- lapply(1:J, d_nb, R_vec = R_vec, nb_size = nb_size, nb_mu = nb_mu)
  tpm_hmm <- construct_tpm(J, cond_tpm, d_r, R_vec)
  
  delta <- solve(t(diag(sum(R_vec)) - tpm_hmm + 1), rep(1, sum(R_vec)))
  
  return(list(mu1 = mu1, mu2 = mu2, sigma1 = sigma1, sigma2 = sigma2,
              nb_size = nb_size, nb_mu = nb_mu,
              cond_tpm = cond_tpm, tpm_hmm = tpm_hmm,
              delta = delta))
}



################################################################################
################################################################################

# 16) AD-safe negative log likelihood for HMM

nll_hmm <- function(dat){
  function(par){
    getAll(par, dat)
    nat <- unc_to_natural_hmm(parameter_vector, J)
    nl_total <- 0
    for (data in data_list){
      allprobs <- get_allprobs(nat$mu1, nat$mu2, nat$sigma1, nat$sigma2, data, nrow(data), J)
      nl_total <- nl_total + my_forward_hsmm(nat$delta, nat$tpm, allprobs, nrow(allprobs))
    }
    mu1 <- nat$mu1
    mu2 <- nat$mu2
    sigma1 <- nat$sigma1
    sigma2 <- nat$sigma2
    tpm <- nat$tpm
    delta <- nat$delta
    REPORT(mu1)
    REPORT(mu2)
    REPORT(sigma1)
    REPORT(sigma2)
    REPORT(tpm)
    REPORT(delta)
    nl_total
  }
}


################################################################################
################################################################################

# 17) AD-safe negative log likelihood for nb-HSMM

nll_nbhsmm <- function(dat){
  function(par){
    getAll(par, dat)
    nat <- unc_to_natural_nbhsmm(parameter_vector, J, R_vec)
    nl_total <- 0
    for (data in data_list){
      allprobs <- get_allprobs_hsmm(nat$mu1, nat$mu2, nat$sigma1, nat$sigma2, data, nrow(data), J, R_vec)
      nl_total <- nl_total + my_forward_hsmm(nat$delta, nat$tpm_hmm, allprobs, nrow(allprobs))
    }
    mu1 <- nat$mu1
    mu2 <- nat$mu2
    sigma1 <- nat$sigma1
    sigma2 <- nat$sigma2
    nb_size <- nat$nb_size
    nb_mu <- nat$nb_mu
    cond_tpm <- nat$cond_tpm
    tpm_hmm <- nat$tpm_hmm
    delta <- nat$delta
    REPORT(mu1)
    REPORT(mu2)
    REPORT(sigma1)
    REPORT(sigma2)
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

# 18) MLE estimate for HMM

mle_hmm_multiple <- function(mu1, mu2, sigma1, sigma2, tpm, J, data_list){
  parameter_vector <- natural_to_unc_hmm(mu1, mu2, sigma1, sigma2, tpm, J)
  
  dat <- list(data_list = data_list, J = J)
  par <- list(parameter_vector = parameter_vector)
  obj <- MakeADFun(nll_hmm(dat), par, silent = TRUE)
  
  model <- nlminb(obj$par, obj$fn, obj$gr,
                  control = list(iter.max = 5000, eval.max = 5000))
  mle_estimate <- unc_to_natural_hmm(model$par, J)
  
  return(list(mu1 = mle_estimate$mu1, mu2 = mle_estimate$mu2,
              sigma1 = mle_estimate$sigma1, sigma2 = mle_estimate$sigma2,
              tpm = mle_estimate$tpm, delta = mle_estimate$delta,
              mll = model$objective, # was model$minimum (in nlm)
              conv = model$convergence, # was model$code
              nit = model$iterations,
              gradient = obj$gr(model$par), # was model$gradient
              estimate = model$par # was model$estimate
  ))
  
  
}



################################################################################
################################################################################

# 19) MLE estimate for nb-HSMM

mle_nbhsmm_multiple <- function(mu1, mu2, sigma1, sigma2, nb_size, nb_mu, cond_tpm, J, data_list, R_vec){
  parameter_vector <- natural_to_unc_nbhsmm(mu1, mu2, sigma1, sigma2, nb_size, nb_mu, cond_tpm, J)
  
  dat <- list(data_list = data_list, J = J, R_vec = R_vec)
  par <- list(parameter_vector = parameter_vector)
  obj <- MakeADFun(nll_nbhsmm(dat), par, silent = TRUE)
  
  model <- nlminb(obj$par, obj$fn, obj$gr,
                  control = list(iter.max = 5000, eval.max = 5000))
  mle_estimate <- unc_to_natural_nbhsmm(model$par, J, R_vec)
  
  indexing_states <- cumsum(c(0, R_vec))
  hsmm_delta <- sapply(1:J, function(i) sum(mle_estimate$delta[indexing_states[i] + (1:R_vec[i]) ]))
  
  return(list(mu1 = mle_estimate$mu1, mu2 = mle_estimate$mu2, sigma1 = mle_estimate$sigma1,
              sigma2 = mle_estimate$sigma2,
              cond_tpm = mle_estimate$cond_tpm, tpm_hmm = mle_estimate$tpm_hmm,
              delta = mle_estimate$delta, hsmm_delta = hsmm_delta,
              nb_size = mle_estimate$nb_size, nb_mu = mle_estimate$nb_mu,
              mll = model$objective,
              conv = model$convergence,
              nit = model$iterations,
              gradient = obj$gr(model$par),
              estimate = model$par
  ))
}



################################################################################
################################################################################

# 20) Rewrite dnbinom in an AD-safe way

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

# 21) Obtain the stationary distribution


# the stationary distribution is by definition the eigen-vector of the transposed tpm
# corresponding to the largest eigenvalue (1)
# otherwise solve(A, b) solves Ax = b, AD-safe

get_stationary <- function(A) {
  eig <- eigen(t(A))
  pi <- Re(eig$vectors[, 1])
  pi / sum(pi)
}


################################################################################
################################################################################

# 22) Get data driven initialization

init_emissions <- function(dat, J){
  
  X <- do.call(rbind, lapply(dat, function(d) as.matrix(d[, c("var1","var2")])))
  X <- X[stats::complete.cases(X), , drop = FALSE]
  
  km <- stats::kmeans(X, centers = J, nstart = 10)
  ord <- order(km$centers[, 1]) # order by var1 mean
  cen <- km$centers[ord, , drop = FALSE]
  
  mu1 <- cen[, 1]
  mu2 <- cen[, 2]
  
  # pooled sd
  sigma1 <- rep(stats::sd(X[, 1]), J)
  sigma2 <- rep(stats::sd(X[, 2]), J)
  
  list(mu1_0 = mu1, mu2_0 = mu2, sigma1_0 = sigma1, sigma2_0 = sigma2)
}


################################################################################
################################################################################

# 23) Handle divergence issues in qreml resulting from huge lambda

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

# 24) A simpler version of the Viterbi algorithm 

# directly computes accuracy
# handles state aggregates from HSMM

my_viterbi_acc <- function(delta, tpm, data_list, J,
                           mu1, mu2, sigma1, sigma2,
                           hidden_states, 
                           model = c("HMM", "HSMM"), R_vec = NULL, perm){
  
  model <- match.arg(model)
  stopifnot(length(data_list) == nrow(hidden_states))
  out <- rep(NA_real_, length(data_list)) 
  
  # add check 
  pars_ok <- all(is.finite(delta)) && all(is.finite(tpm)) &&
    all(is.finite(c(mu1, mu2, sigma1, sigma2))) &&
    all(sigma1 > 0) && all(sigma2 > 0)
  if (!pars_ok){return(out)}
  
  # to match expanded state space
  if (model == "HSMM"){
    state_map <- rep(seq_len(J), times = R_vec)
  }else{
    state_map <- seq_len(J)
  }
  
  accuracy_per_subject <- numeric(length(data_list))
  
  for (i in seq_along(data_list)){
    
    data <- data_list[[i]]
    true_states <- hidden_states[i, ]
    T_len <- nrow(data)
    
    if (model == "HSMM"){
      allprobs <- get_allprobs_hsmm(mu1, mu2, sigma1, sigma2, data, T_len, J, R_vec)
    }else{
      allprobs <- get_allprobs(mu1, mu2, sigma1, sigma2, data, T_len, J)
    }
    
    n_states <- ncol(allprobs)
    stopifnot(n_states == nrow(tpm), n_states == length(delta))
    
    xi <- matrix(0,  nrow = n_states, ncol = T_len)   # scaled Viterbi values
    backtrace <- matrix(0L, nrow = n_states, ncol = T_len)
    
    # initialise + scale
    v <- delta * allprobs[1, ]
    s <- sum(v)
    if (s == 0){s <- 1}                  
    xi[, 1] <- v/s
    
    for (t in 2:T_len){
      M <- xi[, t-1] * tpm
      backtrace[, t] <- apply(M, 2, which.max)
      v <- allprobs[t, ] * apply(M, 2, max)
      s <- sum(v)
      if (s == 0){s <- 1}
      xi[, t] <- v/s
    }
    
    # back-trace
    best_path <- integer(T_len)
    best_path[T_len] <- which.max(xi[, T_len])
    for (t in (T_len - 1):1){
      best_path[t] <- backtrace[best_path[t + 1], t + 1]
    }
    
    decoded_fitted <- state_map[best_path]
    # use permutation
    decoded_true <- perm[decoded_fitted]
    accuracy_per_subject[i] <- mean(decoded_true == true_states)
  }
  
  accuracy_per_subject
  
}


################################################################################
################################################################################

# 25) Extract parameters of interest

# only extract the relevant parameters from the final model fit to save memory on HPC

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
    mu1 <- model_fit$mu1
    mu2 <- model_fit$mu2
    sigma1 <- model_fit$sigma1
    sigma2 <- model_fit$sigma2
    
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
      mu1 = mu1, mu2 = mu2, sigma1 = sigma1, sigma2 = sigma2,
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
    mu1 <- model_fit$mu1
    mu2 <- model_fit$mu2
    sigma1 <- model_fit$sigma1
    sigma2 <- model_fit$sigma2
    
    # log likelihood
    ll <- -(model_fit$mll)
    
    return(list(
      delta = delta,
      hsmm_delta = hsmm_delta,
      d_r = d_r,
      p_list = p_list,
      cond_tpm = cond_tpm,
      hmm_tpm = hmm_tpm,
      mu1 = mu1, mu2 = mu2, sigma1 = sigma1, sigma2 = sigma2,
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
    
    mu1 <- model_fit$mu1
    mu2 <- model_fit$mu2
    sigma1 <- model_fit$sigma1
    sigma2 <- model_fit$sigma2
    
    ll <- -(model_fit$mll)
    
    
    return(list(
      delta = delta,
      hsmm_delta = hsmm_delta,
      d_r = d_r,
      nb_size = nb_size,
      nb_mu = nb_mu,
      cond_tpm = cond_tpm,
      hmm_tpm = hmm_tpm,
      mu1 = mu1, mu2 = mu2, sigma1 = sigma1, sigma2 = sigma2,
      ll = ll
    ))
    
  }
  
  if (model == "hmm"){
    
    delta <- model_fit$delta
    
    hmm_tpm <- model_fit$tpm
    
    mu1 <- model_fit$mu1
    mu2 <- model_fit$mu2
    
    sigma1 <- model_fit$sigma1
    sigma2 <- model_fit$sigma2
    
    ll <- -(model_fit$mll)
    
    return(list(
      delta = delta,
      hmm_tpm = hmm_tpm,
      mu1 = mu1, mu2 = mu2, sigma1 = sigma1, sigma2 = sigma2,
      ll = ll
    ))
    
  }
}




################################################################################
################################################################################

# 26) Safe function call for HPC

safe <- function(expr) tryCatch(expr, error = function(e) NULL)

################################################################################
################################################################################

# 27) Multistart the HMM and keep the best one

fit_hmm_multistart <- function(em, tpm_0, J, data_list, n_start = 4, pert_sd = 0.5) {
  best_f <- NULL
  best_ll <- -Inf
  
  for (s in seq_len(n_start)) {
    if (s == 1L) { # just use data based emissions without permutation if multistart set to 1 (no multistart)
      mu1 <- em$mu1_0
      mu2 <- em$mu2_0
    } else {
      mu1p <- em$mu1_0 + rnorm(J, 0, pert_sd)
      mu2p <- em$mu2_0 + rnorm(J, 0, pert_sd)
      ord  <- order(mu1p) # keep ascending order
      mu1  <- mu1p[ord]
      mu2 <- mu2p[ord]
    }
    f <- safe(mle_hmm_multiple(mu1, mu2, em$sigma1_0, em$sigma2_0, tpm_0, J, data_list))
    if (is.null(f)) next
    p <- extract_relevant_pars(f, model = "hmm", J = J)
    
    ll <- if (is.null(p)){
      NULL}else{p$ll}
    
    if (is.null(ll) || !is.finite(ll)) next
    if (ll > best_ll) { 
      best_ll <- ll # keep HMM with highest log likelihood
      best_f <- f }
  }
  best_f
}


################################################################################
################################################################################

# 28) Factory negative log likelihood for penalized HSMM

# use this when lambda is given (CV, AIC, etc. for each candidate lambda)

make_pen_nll <- function(dat) {
  force(dat) # evaluate dat now
  function(par) {
    getAll(par, dat) # exposes unc_p, S, lambda, data_list, J, R_vec
    nat <- unc_to_natural_nonp_hsmm(par, J, R_vec)
    pen <- penalty(unc_p, S, lambda)
    
    nl_total <- 0
    for (data in data_list){
      allprobs <- get_allprobs_hsmm(nat$mu1, nat$mu2, nat$sigma1, nat$sigma2,
                                    data, nrow(data), J, R_vec)
      nl_total <- nl_total + my_forward_hsmm(nat$delta, nat$tpm_hmm, allprobs, nrow(allprobs))
    }
    
    p_list <- nat$p_list
    d_r <- nat$d_r
    cond_tpm <- nat$cond_tpm
    delta <- nat$delta
    mu1 <- nat$mu1
    mu2 <- nat$mu2
    sigma1 <- nat$sigma1
    sigma2 <- nat$sigma2
    REPORT(p_list)
    REPORT(d_r)
    REPORT(cond_tpm)
    REPORT(delta)
    REPORT(mu1)
    REPORT(mu2)
    REPORT(sigma1)
    REPORT(sigma2)
    
    nl_total + pen
  }
}



################################################################################
################################################################################


# 29) MLE for the nonparametric HSMM, given lambda (needs factory pen nll)

# bound for numerical stability issues
mle_nonp_multiple <- function(mu1, mu2, sigma1, sigma2, p_list, cond_tpm,
                              J, data_list, R_vec, lambda, S, bound = 30){
  par <- natural_to_unc_nonp_hsmm(mu1, mu2, sigma1, sigma2, p_list, cond_tpm, J)
  dat <- list(data_list = data_list, J = J, R_vec = R_vec, S = S, lambda = lambda)
  
  obj <- MakeADFun(make_pen_nll(dat), par, silent = TRUE)   # fixed lambda
  model <- nlminb(obj$par, obj$fn, obj$gr,
                  lower = -bound, upper = bound,
                  control = list(iter.max = 5000, eval.max = 5000))
  mle_estimate <- unc_to_natural_nonp_hsmm(obj$env$parList(model$par), J, R_vec)
  
  indexing_states <- cumsum(c(0, R_vec))
  hsmm_delta <- sapply(1:J, function(i) sum(mle_estimate$delta[indexing_states[i] + (1:R_vec[i])]))
  
  list(mu1 = mle_estimate$mu1, mu2 = mle_estimate$mu2,
       sigma1 = mle_estimate$sigma1, sigma2 = mle_estimate$sigma2,
       cond_tpm = mle_estimate$cond_tpm, tpm_hmm = mle_estimate$tpm_hmm,
       delta = mle_estimate$delta, hsmm_delta = hsmm_delta,
       p_list = mle_estimate$p_list, d_r = mle_estimate$d_r,
       mll = model$objective, conv = model$convergence, nit = model$iterations,
       gradient = obj$gr(model$par), estimate = model$par)
}


################################################################################
################################################################################

# 30) Functions for CV validated lambda

# create the time series with blanked out points as described by Langrock (2015)
# outputs N x T_len matrix indicating what fold each point belongs to
# for now splits whole dataset, maybe build subset version later

make_cv_folds <- function(N, T_len, K = 5, seed = 1L) {
  
  # leave global seed untouched for reproducibility
  old <- if (exists(".Random.seed", .GlobalEnv)) get(".Random.seed", .GlobalEnv) else NULL
  on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = .GlobalEnv), add = TRUE)
  set.seed(seed)
  
  fold <- sample(rep_len(1:K, N * T_len)) # sample random assignment to fold
  matrix(fold, nrow = N, ncol = T_len) # put into matrix form
}

# 30.1)

# set held out points to NA
# mask is logical, 1 = blanked
blank_points <- function(data_list, mask) {
  lapply(seq_along(data_list), function(n) {
    d <- data_list[[n]]
    na_t <- which(mask[n, ])
    if (length(na_t)) d[na_t, ] <- NA # blank both var1 and var2
    d # local copy
  })
}

# 30.2)

# fit nonp with warm start, lambda fixed
# start with unconstrained parameter list (no conversions, save time + not necessary)
# bound for numerical stability issues
fit_nonp_at_lambda <- function(par_init, data_list, J, R_vec, S, lambda, bound = 30) {
  dat <- list(data_list = data_list, J = J, R_vec = R_vec, S = S, lambda = lambda)
  obj <- MakeADFun(make_pen_nll(dat), par_init, silent = TRUE)
  m <- nlminb(obj$par, obj$fn, obj$gr,
              lower = -bound, upper = bound,
              control = list(iter.max = 5000, eval.max = 5000))
  list(par_list = obj$env$parList(m$par), mll = m$objective, conv = m$convergence)
}

# 30.3)

# evaluate on held-out points
# data list score has calibration points blanked
# evaluate on validation points
# return unpenalized log likelihood of fitted model with lambda
score_validation <- function(par_list, data_list_score, J, R_vec, S) {
  dat <- list(data_list = data_list_score, J = J, R_vec = R_vec,
              S = S, lambda = rep(0, J))
  obj <- MakeADFun(make_pen_nll(dat), par_list, silent = TRUE)
  -as.numeric(obj$fn(obj$par))
}

# 30.4)

# look up warm start
nearest_cached <- function(fold_fits, target_key, seed_par) {
  keys <- names(fold_fits)
  if (length(keys) == 0L) return(seed_par)
  # parse target key back into numeric vector
  tgt <- as.numeric(strsplit(target_key, "_", fixed = TRUE)[[1]])
  # get distance to neighboring keys, choose closest one for warm starting new fit
  d <- vapply(keys, function(kk)
    sum(abs(as.numeric(strsplit(kk, "_", fixed = TRUE)[[1]]) - tgt)), numeric(1))
  fold_fits[[ keys[ which.min(d) ] ]]
}


################################################################################
################################################################################

# 31) CV score

# get full CV score for one lambda across all folds
cv_score_lambda <- function(exp_vec, exp_key, fold_id, K, data_list,
                            J, R_vec, S, cache, seed_par) {
  lambda <- 10^exp_vec
  total <- 0
  for (k in 1:K) {
    val_mask <- (fold_id == k) # mark fold as validation set
    data_fit <- blank_points(data_list, val_mask) # set validation NA
    # blank calibration points (score on this, the validation points)
    data_score <- blank_points(data_list, !val_mask) 
    # choose warm start from this folds cache
    init <- nearest_cached(cache$fit[[k]], exp_key, seed_par)
    
    fit <- tryCatch(fit_nonp_at_lambda(init, data_fit, J, R_vec, S, lambda),
                    error = function(e) NULL)
    if (is.null(fit) || fit$conv != 0) return(-Inf)
    
    # needed for fits with future lambdas to warm start this fold (save time!)
    cache$fit[[k]][[exp_key]] <- fit$par_list # persists (cache is env)
    
    # score on validation points
    s <- tryCatch(score_validation(fit$par_list, data_score, J, R_vec, S),
                  error = function(e) NA_real_)
    if (!is.finite(s)) return(-Inf)
    total <- total + s
  }
  mean_ll <- total/K # average log likelihood across all folds
  return(mean_ll)
}


################################################################################
################################################################################

# 32) Greedy neighboring algorithm for CV

greedy_lambda_cv <- function(data_list, J, R_vec, S, par_init,
                             exp_grid = 0:6,
                             start_exp = rep(2, J),
                             K = 5, fold_seed = 1L, max_iter = 30) {
  
  N <- length(data_list)
  T_len <- nrow(data_list[[1]]) # assumes same T_len
  # important - fold partition fixed across lambda otherwise not comparable
  fold_id <- make_cv_folds(N, T_len, K, seed = fold_seed) # build the fold partition once
  
  # fit once on full data to get initial warm start + fall-back warm start
  # par_init must be already a warm start from nb-HSMM!
  seed_fit <- fit_nonp_at_lambda(par_init, data_list, J, R_vec, S, 10^start_exp)
  seed_par <- seed_fit$par_list
  
  cache <- new.env(parent = emptyenv())
  cache$score <- list() # exp_key = CV score if computed already
  cache$fit <- replicate(K, list(), simplify = FALSE) # per fold, to use key to get par list
  
  key_of <- function(e) paste(e, collapse = "_")
  
  # only recompute if not done yet, otherwise return CV score by key
  eval_exp <- function(e) {
    kk <- key_of(e)
    if (!is.null(cache$score[[kk]])){
      return(cache$score[[kk]])
      }else{
    sc <- cv_score_lambda(e, kk, fold_id, K, data_list, J, R_vec, S, cache, seed_par)
    cache$score[[kk]] <- sc
    return(sc)}
  }
  
  # move one exp up or down by one, keeping the other two fixed
  # that is 6 neighbors (2 per state, only direct neighbors)
  neighbours <- function(e) {
    out <- list()
    for (j in seq_len(J)) for (step in c(-1L, 1L)) {
      e2 <- e
      e2[j] <- e[j] + step
      if (e2[j] >= min(exp_grid) && e2[j] <= max(exp_grid))
        out[[length(out) + 1L]] <- e2
    }
    out
  }
  
  cur <- start_exp
  cur_sc <- eval_exp(cur)
  trace <- list(list(exp = cur, score = cur_sc))
  
  for (it in seq_len(max_iter)) {
    nb <- neighbours(cur)
    nb_sc <- vapply(nb, eval_exp, numeric(1))
    best <- which.max(nb_sc)
    if (length(best) == 0L || nb_sc[best] <= cur_sc) break   # no strict improvement
    cur <- nb[[best]]
    cur_sc <- nb_sc[best]
    trace[[length(trace) + 1L]] <- list(exp = cur, score = cur_sc)
  }
  
  list(lambda = 10^cur,
       exp = cur,
       cv_score = cur_sc,
       trace = trace, # inspect for oscillation, raise K if noisy
       n_evaluated = length(cache$score),
       fold_id = fold_id,
       seed_par = seed_par) # return initial warm start for full model fit once we have lambda
}



################################################################################
################################################################################

# 33) Compute edf for AIC

compute_edf_aic <- function(par_list, data_list, J, R_vec, S, lambda) {
  
  # penalized information H + lambda*S, at theta_hat
  dat_pen <- list(data_list = data_list, J = J, R_vec = R_vec, S = S, lambda = lambda)
  obj_pen <- MakeADFun(make_pen_nll(dat_pen), par_list, silent = TRUE)
  theta <- obj_pen$par
  H_pen <- obj_pen$he(theta)  # = H + lambda*S
  
  # unpenalized ll and Hessian at the same theta_hat (lambda = 0)
  dat_unp <- list(data_list = data_list, J = J, R_vec = R_vec, S = S, lambda = rep(0, J))
  obj_unp <- MakeADFun(make_pen_nll(dat_unp), par_list, silent = TRUE)
  ll <- -as.numeric(obj_unp$fn(theta))
  H_exact <- obj_unp$he(theta)
  
  # flat-index of each state's unc_p block, for per-state edf
  idx_tmpl <- obj_pen$env$parList(seq_along(theta))
  unc_p_index <- lapply(1:J, function(i) as.integer(idx_tmpl$unc_p[i, ]))
  
  # edf defined as trace of product of inverse penalized Hessian times unpenalized Hessian at theta hat
  ridge <- 1e-6 * mean(abs(diag(H_pen))) # stability
  M <- solve(H_pen + diag(ridge, nrow(H_pen)), H_exact)
  dM <- diag(M)
  
  edf <- sum(dM)
  edf_per_state <- vapply(unc_p_index, function(ix) sum(dM[ix]), numeric(1))
  
  list(
    ll = ll,
    npar = length(theta),
    lambda = lambda,
    edf = edf,
    edf_per_state = edf_per_state,
    aic = -2 * ll + 2 * edf
  )
}


################################################################################
################################################################################

# 34) Greedy neighboring algorithm for AIC

greedy_lambda_aic <- function(data_list, J, R_vec, S, par_init,
                              exp_grid = 0:6, start_exp = rep(2, J),
                              max_iter = 30) {
  
  # warm start, use NB warm start for par_init
  seed_fit <- fit_nonp_at_lambda(par_init, data_list, J, R_vec, S, 10^start_exp)
  seed_par <- seed_fit$par_list
  
  cache <- new.env(parent = emptyenv())
  cache$aic <- list()                     
  cache$fit <- list()                          
  cache$info <- list()                                
  
  key_of <- function(e) paste(e, collapse = "_")
  
  nearest <- function(target_key) { # nearest already fitted lambda
    keys <- names(cache$fit)
    
    if (length(keys) == 0L){
      return(seed_par)}
    
    tgt <- as.numeric(strsplit(target_key, "_", fixed = TRUE)[[1]])
    d <- vapply(keys, function(kk)
      sum(abs(as.numeric(strsplit(kk, "_", fixed = TRUE)[[1]]) - tgt)), numeric(1))
    cache$fit[[keys[which.min(d)]]]
  }
  
  eval_exp <- function(e) {
    kk <- key_of(e)
    if (!is.null(cache$aic[[kk]])) return(cache$aic[[kk]])
    lambda <- 10^e
    fit <- tryCatch(fit_nonp_at_lambda(nearest(kk), data_list, J, R_vec, S, lambda),
                    error = function(err) NULL)
    
    # check if all parameters are finite but bound worked, then conv = 1 is ok
    # otherwise diverged fit, shall discard
    ok <- !is.null(fit) && fit$conv %in% c(0L, 1L) && all(is.finite(unlist(fit$par_list)))
    if (!ok) {
      cache$aic[[kk]] <- Inf
      return(Inf) }
    
    cache$fit[[kk]] <- fit$par_list
    
    info <- tryCatch(
      compute_edf_aic(fit$par_list, data_list, J, R_vec, S, lambda),
      error = function(err) NULL)
    if (is.null(info) || !is.finite(info$aic)) { 
      cache$aic[[kk]] <- Inf
      return(Inf) }
    
    cache$info[[kk]] <- info
    cache$aic[[kk]] <- info$aic
    info$aic
  }
  
  # same as for greedy CV
  neighbours <- function(e) {
    out <- list()
    for (j in seq_len(J)) for (step in c(-1L, 1L)) {
      e2 <- e
      e2[j] <- e[j] + step
      if (e2[j] >= min(exp_grid) && e2[j] <= max(exp_grid))
        out[[length(out) + 1L]] <- e2
    }
    out
  }
  
  cur <- start_exp
  cur_aic <- eval_exp(cur)
  trace <- list(list(exp = cur, aic = cur_aic))
  
  for (it in seq_len(max_iter)) {
    nb <- neighbours(cur)
    nb_aic <- vapply(nb, eval_exp, numeric(1))
    best <- which.min(nb_aic)       # smaller AIC is better
    if (length(best) == 0L || nb_aic[best] >= cur_aic) break
    cur <- nb[[best]]
    cur_aic <- nb_aic[best]
    trace[[length(trace) + 1L]] <- list(exp = cur, aic = cur_aic)
  }
  
  list(lambda = 10^cur,
       exp = cur,
       aic = cur_aic,
       info = cache$info[[key_of(cur)]],      
       trace = trace,
       n_evaluated = length(cache$aic),
       seed_par = seed_par)
}


################################################################################
################################################################################

# 36) Emission-based state alignment

# compute state permutation once based on closest means, then use for both Viterbi and MISE
# emission-based permutation standard and not biased towards any single outcome measure

align_states <- function(fit_mu1, fit_mu2, true_mu1, true_mu2){
  J <- length(true_mu1)
  # cost[i, j] = mismatch between fitted state i and true state j
  cost <- matrix(0, J, J)
  for (i in 1:J){
    for (j in 1:J){
      cost[i, j] <- (fit_mu1[i] - true_mu1[j])^2 + (fit_mu2[i] - true_mu2[j])^2
    }
  }
  # Hungarian to minimize total cost
  perm <- as.integer(clue::solve_LSAP(cost))
  return(perm)
}


################################################################################
################################################################################

# 37) Build the dwell time distribution as inferred by the nonparametric HSMM

# this accounts for the unstructured start (free parameters estimated up to R_vec i)
# followed by the geometric tail
# as explained in thesis and construction of tpm (by Pohle et al.)

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

# 38) Build dwell-time matrix to compare to truth (DGP fed to simulator)

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

# 39) Calculate MSSE 
# (mean summed squared error, prev labelled mean integrated squared error but adjusted to fit discrete distribution,
# hence the name mismatch)

dwell_mise <- function(dist_d, pars, model, J, R_vec, perm, R_eval = ncol(dist_d)){
  
  fit_D <- build_dwell_matrix(pars, model, J, R_vec, R_eval)
  
  # true dwell time 
  # should already be normalized if using same dist d but safeguard
  true_D <- dist_d[, 1:R_eval, drop = FALSE]
  true_D <- true_D/rowSums(true_D)
  
  ise <- numeric(J)
  for (i in 1:J){
    j <- perm[i] # use same permutation
    
    # calculate integrated squared error
    ise[j] <- sum((fit_D[i, ] - true_D[j, ])^2)
  }
  
  names(ise) <- paste0("true_state_", 1:J)
  list(per_state = ise, total = sum(ise))
}


################################################################################
################################################################################


# 40) small helper for convergence safeguard

rll <- function(f) {
  v <- f$llk_restricted
  if (is.null(v) || length(v) == 0) NA_real_ else v[length(v)]   # last = converged
}



