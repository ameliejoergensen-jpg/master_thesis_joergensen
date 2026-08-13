
# stores the emission and dwell specification
source("git_master_store_scenarios_complex.R")

scenarios <- expand.grid(
  T_len = 1000, # length of time series
  N = c(30), # subjects
  order_diff = 2,
  sojourn = names(sojourn_distributions),
  emission = names(emission_distributions),
  stringsAsFactors = FALSE
)

scenarios$scenario_id <- seq_len(nrow(scenarios))
scenarios$seed <- 5000 + scenarios$scenario_id # assign seed to each scenario for reproducibility

# write the csv specifying the simulation scenarios
write.csv(scenarios, "master_scenarios_complex.csv", row.names = FALSE)




