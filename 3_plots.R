library(tidyverse)
source("3_plot_functions.R")

results <- readRDS("results/temp_results.rds")
results <- readRDS("results/results_complex_covars.rds")

bin <- summarize_binary(results)
cont <- summarize_continuous(results)

generate_TIE_plot(cont)
generate_TIE_plot(bin)

generate_coverage_plot(cont)
generate_coverage_plot(bin)

generate_bias_plot(cont)
generate_bias_plot(bin)

generate_power_plot(cont)
generate_power_plot(bin)
















