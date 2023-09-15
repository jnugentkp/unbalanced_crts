source("2_run_sims.R")

nsims <- 500
results_complex_covars <- run_full_set(sim_type = "complex_covars",
                                       cores = 30,
                                       nsims = nsims, seed = NA)
saveRDS(results_complex_covars, "results/results_complex_covars.rds")

results_complex_covars <- run_full_set(sim_type = "main_terms",
                                       cores = 30,
                                       nsims = nsims, seed = NA)
saveRDS(results_complex_covars, "results/results_main_terms.rds")

results_complex_covars <- run_full_set(sim_type = "overfit",
                                       cores = 30,
                                       nsims = nsims, seed = NA)
saveRDS(results_complex_covars, "results/results_overfit.rds")

results_complex_covars <- run_full_set(sim_type = "HTE",
                                       cores = 30,
                                       nsims = nsims, seed = NA)
saveRDS(results_complex_covars, "results/results_HTE.rds")

# Small edits
# ACI now.