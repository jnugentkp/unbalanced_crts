require(tidyverse)
`%nin%` <- negate(`%in%`)


generate_summary_measures <- function(LP0, LP1, binary = F){
  if(binary){
    P0 <- plogis(LP0)
    P1 <- plogis(LP1)
    glmm_logistic_regression_beta <- mean(LP1 - LP0)
    glmm_true_cond_or <- exp(glmm_logistic_regression_beta)
    gee_true_marg_or <- glmm_true_cond_or # or same as tmle_true_marg_or? seems not
    tmle_true_marg_or <- (mean(P1) / (1-mean(P1))) / (mean(P0) / (1-mean(P0)))
    tmle_true_marg_rd <- true_marg_mean_diff <- mean(P1 - P0)
    tmle_true_marg_rr <- mean(P1) / mean(P0)
  } else {
    true_marg_mean_diff <- mean(LP1) - mean(LP0)
    glmm_logistic_regression_beta <- NA
    glmm_true_cond_or <- NA
    gee_true_marg_or <- NA
    tmle_true_marg_or <- NA
    tmle_true_marg_rd <- NA
    tmle_true_marg_rr <- NA
  }
  return(cbind.data.frame(true_marg_mean_diff = true_marg_mean_diff,
                          glmm_logistic_regression_beta = glmm_logistic_regression_beta,
                          glmm_true_cond_or = glmm_true_cond_or,
                          gee_true_marg_or = gee_true_marg_or,
                          tmle_true_marg_or = tmle_true_marg_or,
                          tmle_true_marg_rd = tmle_true_marg_rd,
                          tmle_true_marg_rr = tmle_true_marg_rr))
}

log_prob_from_linear_predictor <- function(LP1, LP0, epsilon = .0001){
  lo_goal_lp <- exp(epsilon)
  hi_goal_lp <- exp(1 - epsilon)
  
  minLP <- min(c(LP0, LP1))
  maxLP <- max(c(LP0, LP1))
  
  LP0 <- LP0 - minLP + epsilon
  LP1 <- LP1 - minLP + epsilon
  
  maxLP <- max(c(LP0, LP1))
  scale_factor <- hi_goal_lp / maxLP
  
  LP0 <- LP0 * scale_factor + 1
  LP1 <- LP1 * scale_factor + 1
  
  return(list(LP1 = LP1, LP0 = LP0))
}

gen_obs_dat <- function(binary = F,
                        sim_type = "main_terms",
                        trt_clust_size = 5, trt_clusts = 50,
                        control_clust_size = 1, control_clusts = 50,
                        txt_eff = .3, ranef_sd = .25, covar1_coef = .25, seed = NA,
                        LPS_only = F){
  if(!is.na(seed)){ set.seed(seed) }
  total_units_txt <- trt_clust_size*trt_clusts
  total_units_ctr <- control_clust_size*control_clusts
  total_units <- total_units_txt + total_units_ctr
  clust_id_txt <- rep(1:trt_clusts, each = trt_clust_size)
  clust_id_ctr <- rep((max(clust_id_txt)+1):(max(clust_id_txt) + control_clusts), each = control_clust_size)
  clust_id <- c(clust_id_txt, clust_id_ctr)
  total_clusts <- max(clust_id)
  clust_lengths <- c(rep(trt_clust_size, times = trt_clusts),
                     rep(control_clust_size, times = control_clusts))
  
  
  A <- c(rep(1, times = total_units_txt), rep(0, times = total_units_ctr))
  if(ranef_sd == 0){
    UE <- 0
  } else {
    UEs <- rnorm(n = total_clusts, mean = 0, sd = ranef_sd)
    UE <- rep(UEs, times = clust_lengths)
  }
  covar0 <- rnorm(n = total_units, mean = UE, sd = 1)
  covar1 <- rbinom(n = total_units, size = 1, prob = .3)
  
  if(binary){
    txt_eff <- 3*txt_eff
    intercept <- -5
    
    if(sim_type == "main_terms"){
      LP0 <- 4.5 + intercept + covar1_coef*covar1 + covar0 + UE
      LP1 <- 4.5 + intercept + covar1_coef*covar1 + covar0 + UE + txt_eff
    } else if(sim_type == "complex_covars"){
      UE <- exp(UE)
      LP0 <- intercept + covar1_coef*covar1 + covar0^2 + covar0*covar1 + UE
      LP1 <- intercept + covar1_coef*covar1 + covar0^2 + covar0*covar1 + UE + txt_eff
    } else if(sim_type == "HTE"){
      LP0 <- 4 + intercept + covar1_coef*covar1 + covar0 + UE
      LP1 <- 4 + intercept + covar1_coef*covar1 + covar0 + UE + txt_eff + 2*txt_eff*covar1
    } else if(sim_type == "overfit"){
      LP0 <- 5 + intercept + UE
      LP1 <- 5 + intercept + UE + .5*txt_eff
    } else {
      print("Simulation structure not entered correctly. Must be main_terms, complex_covars, HTE, or overfit")
    }
    
    
    # IF LOG SCALE... (rather than logistic)
    # trans <- log_prob_from_linear_predictor(LP1 = LP1, LP0 = LP0)
    # LP0 <- trans$LP0
    # LP1 <- trans$LP1
    # P0 <- log(LP0)
    # P1 <- log(LP1)
    
    
    P0 <- plogis(LP0)
    P1 <- plogis(LP1)
    #print(summary(P0))
    #print(summary(P1))
    
    UY <- runif(total_units)
    Y0 <- as.numeric(P0 > UY)
    Y1 <- as.numeric(P1 > UY)
    Y <- ifelse(A == 1, Y1, Y0)
    
  } else {
    # Continuous outcomes
    UY <- rnorm(n = total_units)
    
    if(sim_type == "main_terms"){
      LP0 <- covar1_coef*covar1 + covar0 + UY + UE
      LP1 <- covar1_coef*covar1 + covar0 + UY + UE + txt_eff
    } else if(sim_type == "complex_covars"){
      UE <- exp(UE)
      LP0 <- covar1_coef*covar1 + covar0^2 + covar0*covar1 + UY + UE
      LP1 <- covar1_coef*covar1 + covar0^2 + covar0*covar1 + UY + UE + txt_eff
    } else if(sim_type == "HTE"){
      # print("HTE, continuous")
      LP0 <- covar1_coef*covar1 + covar0 + UY + UE
      LP1 <- covar1_coef*covar1 + covar0 + UY + UE + txt_eff + 2*txt_eff*covar1*covar0
    } else if(sim_type == "overfit"){
      LP0 <- UY + UE
      LP1 <- UY + UE + txt_eff
    } else {
      print("Simulation structure not entered correctly. Must be main_terms, complex_covars, HTE, or overfit")
    }
    
    #print(summary(LP0))
    #print(summary(LP1))
    Y <- ifelse(A == 1, LP1, LP0)  
  }
  
  if(LPS_only){
    return(cbind.data.frame(LP0 = LP0, LP1 = LP1))
  }
  return(cbind.data.frame(clust_id, A, covar0, covar1, Y))
}





#head(gen_all(binary = T, sim_type = "overfit"))
#head(gen_all(binary = T, sim_type = "overfit", LPS_only = T))#p <- gen_all(binary = F)

