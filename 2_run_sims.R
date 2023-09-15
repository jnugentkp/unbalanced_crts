library(tidyverse)
library(ltmle)
library(doParallel)
library(foreach)
library(sandwich)
library(lme4)
library(gee)
source("0_DGP.R")
source("1_fit_models.R")

run_full_set <- function(sim_type = "complex_covars", nsims = 5, cores = 8, seed = NA){
  if(!is.na(seed)){
    set.seed(seed)
  }
  
  trt_clusts <- c(25, 50)
  trt_clust_size <- c(5, 10)
  control_clusts <- c(25, 50)
  control_clust_size <- c(1, 5, 10)
  txt_eff <- c(.25, 0, -.25)
  v <- 20
  binary <- c(T, F)
  sim_type <- sim_type
  params <- expand_grid(trt_clusts, trt_clust_size,
                        control_clusts, control_clust_size,
                        txt_eff, binary, v, sim_type)
  params
  
  one_set <- function(DGP_specification, nsims = 5, cores = cores){
    
    # Get true PATE values by averaging...
    pop <- bind_rows(lapply(as.list(1:1000), function(x) gen_obs_dat(LPS_only = T,
                                                                     trt_clusts = DGP_specification[[1,"trt_clusts"]],
                                                                     trt_clust_size = DGP_specification[[1,"trt_clust_size"]],
                                                                     control_clusts = DGP_specification[[1,"control_clusts"]],
                                                                     control_clust_size = DGP_specification[[1,"control_clust_size"]],
                                                                     txt_eff = DGP_specification[[1,"txt_eff"]],
                                                                     binary = DGP_specification[[1,"binary"]],
                                                                     sim_type = DGP_specification[[1,"sim_type"]])))
    PATE <- generate_summary_measures(LP0 = pop$LP0, LP1 = pop$LP1, binary = DGP_specification[[1,"binary"]])
    
    registerDoParallel(cores = cores)
    results <- data.frame(foreach(j = 1:nsims, .combine = rbind, .errorhandling = "remove") %dopar% {
      
      fit_models(dat = gen_obs_dat(trt_clusts = DGP_specification[[1,"trt_clusts"]],
                                   trt_clust_size = DGP_specification[[1,"trt_clust_size"]],
                                   control_clusts = DGP_specification[[1,"control_clusts"]],
                                   control_clust_size = DGP_specification[[1,"control_clust_size"]],
                                   txt_eff = DGP_specification[[1,"txt_eff"]],
                                   binary = DGP_specification[[1,"binary"]],
                                   sim_type = DGP_specification[[1,"sim_type"]]),
                 binary = DGP_specification[[1,"binary"]],
                 control_clust_size = DGP_specification[[1,"control_clust_size"]],
                 v = DGP_specification[[1,"v"]])
    })
    stopImplicitCluster()
    return(cbind.data.frame(results, PATE, DGP_specification, cores = cores, nsims = nsims))
  }

  return(data.frame(foreach(j = 1:nrow(params), .combine = rbind) %do%
                          {one_set(DGP_specification = params[j,], nsims = nsims)}))
}





