# Simulation scenarios
# This script stores the specification for the dwell-time distributions and emission distributions


###############################################################################
# custom nonparametric peaked distribution

D <- 200 # same as support used for simulation
J <- 3
R <- 20


softspike <- function(m, s, eps = 0.005) {
  d <- dnorm(1:D, m, s)
  d <- d/sum(d)
  fl <- c(rep(1/(2*R), 2*R), numeric(D - 2*R))
  d <- (1 - eps)*d + eps*fl
  d/sum(d)
}

###############################################################################

sojourn_distributions <- list(
  all_geometric = list(dgeom(0:199, p = 0.05), dgeom(0:199, p = 0.1), dgeom(0:199, p = 0.2)),
  
  # mu of nb will be 15, as this gets shifted
  one_nb = list(dnbinom(0:199, size = 60, mu = 14), dgeom(0:199, p = 0.1), dgeom(0:199, p = 0.2)),
  
  # Poisson gets shifted
  one_bimodal = list(0.4 * dgeom(0:199, p = 0.5) + 0.6 * dpois(0:199, lambda = 14), 
                     dgeom(0:199, p = 0.1), dgeom(0:199, p = 0.2)),
  
  # nb gets shifted again
  one_nb_one_spike = list(dnbinom(0:199, mu = 7, size = 60), softspike(16, 1), dgeom(0:199, p = 0.1)),
  
  spikes = list(softspike(5, 2), softspike(10, 2), softspike(15, 2)))



emission_distributions <- list(
  medium_overlap = list(mu1 = c(0, 1, 2), mu2 = c(2, 3, 4), sigma1 = c(1, 1, 1), sigma2 = c(1, 1, 1)),
  high_overlap = list(mu1 = c(0, 0.7, 1.4), mu2 = c(2, 2.7, 3.4), sigma1 = c(1, 1, 1), sigma2 = c(1, 1, 1)))







