source("2_run_sims.R")

nsims <- 5
results_complex_covars <- run_full_set(sim_type = "complex_covars", nsims = nsims, seed = NA)
saveRDS(results_complex_covars, "results/results_complex_covars.rds")




# Small edits
# ACI now.