library(tidyverse)
library(ltmle)
library(doParallel)
library(foreach)
library(sandwich)
library(lme4)
library(gee)
library(lmerTest)
library(lmtest)

fit_models <- function(dat,
                       binary = F,
                       v = 20,
                       control_clust_size,
                       SLL = c("SL.mean", "SL.glm", "SL.earth", "SL.glm.interaction")){
  
  nclust <- length(unique(dat$clust_id))
  
  n <- nrow(dat)
  
  if(binary){
    ###############################################################################################
    # TMLE
    ltmle_mod <- suppressWarnings(suppressMessages(
      ltmle(data = dat %>% select(all_of(c("covar0", "covar1", "A", "Y"))),
            Ynodes = "Y", Anodes = "A",
            id = dat$clust_id,
            variance.method = "tmle", ####################################################"ic",
            SL.library = SLL,
            SL.cvControl = list(V = v),
            abar = list(1,0))))
    summary_ltmle_mod <- summary(ltmle_mod)
    ltmle_or_pe <- summary_ltmle_mod$effect.measures$OR$estimate
    ltmle_or_se <- summary_ltmle_mod$effect.measures$OR$std.dev
    #ltmle_lo <- summary_ltmle_mod$effect.measures$OR$CI[1]
    #ltmle_hi <- summary_ltmle_mod$effect.measures$OR$CI[2]
    # Use nclust - 2 instead...
    ltmle_or_lo <- exp(log(ltmle_or_pe) - ltmle_or_se*qt(.975, df = nclust - 2))
    ltmle_or_hi <- exp(log(ltmle_or_pe) + ltmle_or_se*qt(.975, df = nclust - 2))
    
    ltmle_rr_pe <- summary_ltmle_mod$effect.measures$RR$estimate
    ltmle_rr_se <- summary_ltmle_mod$effect.measures$RR$std.dev
    ltmle_rr_lo <- exp(log(ltmle_rr_pe) - ltmle_rr_se*qt(.975, df = nclust - 2))
    ltmle_rr_hi <- exp(log(ltmle_rr_pe) + ltmle_rr_se*qt(.975, df = nclust - 2))
    
    ltmle_pe <- summary_ltmle_mod$effect.measures$ATE$estimate
    ltmle_se <- summary_ltmle_mod$effect.measures$ATE$std.dev
    ltmle_lo <- ltmle_pe - ltmle_se*qt(.975, df = nclust - 2)
    ltmle_hi <- ltmle_pe + ltmle_se*qt(.975, df = nclust - 2)
    
    
    ###############################################################################################
    # GEE 
    gee <- tryCatch(
      # 'tryCatch()' will return the last evaluated expression 
      # in case the "try" part was completed successfully
      glm(Y ~ A + covar0 + covar1,# id = dat$clust_id, corstr = "exchangeable",
          #maxiter = 100,#
          data = dat, family = binomial(link = "logit")) %>% coeftest(vcov = sandwich),
      # gee(Y ~ A + covar0 + covar1, id = dat$clust_id, corstr = "exchangeable",
      #     maxiter = 100,
      #     data = dat, family = binomial(link = "logit")),
      error = function(e){
        print(e$message)
        return(NULL)
      }
    )    
    if(is.null(gee)){
      gee_pe <- NA
      gee_se <- NA
      gee_lo <- NA
      gee_hi <- NA
    } else {
      # gee_pe <- exp(summary(gee)$coefficients["A","Estimate"])
      # gee_se <- summary(gee)$coefficients["A","Robust S.E."]
      # gee_lo <- exp(log(gee_pe) - qnorm(0.975) * gee_se)
      # gee_hi <- exp(log(gee_pe) + qnorm(0.975) * gee_se)
      gee_pe <- exp(gee["A","Estimate"])
      gee_se <- gee["A","Std. Error"]
      gee_lo <- exp(log(gee_pe) - qnorm(0.975) * gee_se)
      gee_hi <- exp(log(gee_pe) + qnorm(0.975) * gee_se)
      
    }
    
    ###############################################################################################
    # GLMM
    glmm <- tryCatch(
      {
        if(control_clust_size == 1){
          glmer(Y ~ A + covar0 + covar1 +
                  (0 + A | clust_id),
                control = glmerControl(optCtrl = list(maxit = 1e8, maxfun = 1e8),
                                       boundary.tol = 1e-6),
                family = binomial(link = "logit"),
                data = dat)
        } else {
          glmer(Y ~ A + covar0 + covar1 +
                  (1 | clust_id),
                control = glmerControl(optCtrl = list(maxit = 1e8, maxfun = 1e8),
                                       boundary.tol = 1e-6),
                family = binomial(link = "logit"),
                data = dat)
        }
      },
      error = function(e){
        print(e$message)
        return(NULL)
      }
    )
    if(is.null(glmm)){
      glmm_pe <- NA
      glmm_se <- NA
      glmm_lo <- NA
      glmm_hi <- NA
    } else{
      glmm_pe <- exp(summary(glmm)$coefficients["A","Estimate"])
      glmm_se <- summary(glmm)$coefficients["A","Std. Error"]
      glmm_lo <- exp(log(glmm_pe) - qnorm(0.975) * glmm_se)
      glmm_hi <- exp(log(glmm_pe) + qnorm(0.975) * glmm_se)
    }
    
    
  } else { # CONTINUOUS OUTCOME
    # Individual level effect, aggregating the IC accounting for clust size - this is \Psi^I
    ltmle_mod <- suppressWarnings(suppressMessages(
      ltmle(data = dat %>% select(all_of(c("covar0", "covar1", "A", "Y"))),
            Ynodes = "Y", Anodes = "A",
            id = dat$clust_id,
            variance.method = "ic",
            SL.library = SLL,
            SL.cvControl = list(V = v),
            abar = list(1,0))))
    summary_ltmle_mod <- summary(ltmle_mod)
    ltmle_pe <- summary_ltmle_mod$effect.measures$ATE$estimate
    ltmle_se <- summary_ltmle_mod$effect.measures$ATE$std.dev
    # Use nclust - 2 instead...
    #ltmle_lo <- summary_ltmle_mod$effect.measures$ATE$CI[1]
    #ltmle_hi <- summary_ltmle_mod$effect.measures$ATE$CI[2]
    ltmle_lo <- ltmle_pe - ltmle_se*qt(.975, df = nclust - 2)
    ltmle_hi <- ltmle_pe + ltmle_se*qt(.975, df = nclust - 2)
    
    ltmle_or_pe <- NA
    ltmle_or_se <- NA
    ltmle_or_lo <- NA
    ltmle_or_hi <- NA
    
    ltmle_rr_pe <- NA
    ltmle_rr_se <- NA
    ltmle_rr_lo <- NA
    ltmle_rr_hi <- NA
    
    
    # GEE 
    gee <- gee(Y ~ A + covar0 + covar1, id = dat$clust_id, corstr = "exchangeable",
               data = dat)
    gee_pe <- summary(gee)$coefficients["A","Estimate"]
    gee_se <- summary(gee)$coefficients["A","Robust S.E."]
    gee_lo <- gee_pe - qnorm(0.975) * gee_se
    gee_hi <- gee_pe + qnorm(0.975) * gee_se
    
    
    # GLMM
    # Note using Nugent paper to find optimal DF for LMM in this situation (Satterthwaite)
    if(control_clust_size == 1){
      glmm <- lmer(Y ~ A + covar0 + covar1 +
                     (0 + A | clust_id),
                   data = dat, REML = T)
    } else {
      glmm <- lmer(Y ~ A + covar0 + covar1 +
                     (1 | clust_id), data = dat, REML = T)
    }
    glmm_pe <- summary(glmm)$coefficients["A","Estimate"]
    glmm_se <- summary(glmm)$coefficients["A","Std. Error"]
    satterthwaite_df <- summary(glmm)$coefficients["A","df"]
    glmm_lo <- glmm_pe - qt(0.975, df = satterthwaite_df) * glmm_se
    glmm_hi <- glmm_pe + qt(0.975, df = satterthwaite_df) * glmm_se
  }

  return(cbind.data.frame(
    #v = v,
    ltmle_pe = ltmle_pe,
    ltmle_se = ltmle_se,
    ltmle_lo = ltmle_lo,
    ltmle_hi = ltmle_hi,
    
    ltmle_or_pe = ltmle_or_pe,
    ltmle_or_se = ltmle_or_se,
    ltmle_or_lo = ltmle_or_lo,
    ltmle_or_hi = ltmle_or_hi,
    
    ltmle_rr_pe = ltmle_rr_pe,
    ltmle_rr_se = ltmle_rr_se,
    ltmle_rr_lo = ltmle_rr_lo,
    ltmle_rr_hi = ltmle_rr_hi,
    
    gee_pe = gee_pe,
    gee_se = gee_se,
    gee_lo = gee_lo,
    gee_hi = gee_hi,
    
    glmm_pe = glmm_pe,
    glmm_se = glmm_se,
    glmm_lo = glmm_lo,
    glmm_hi = glmm_hi))
}
#run_sim(binary = F, sim_type = "HTE")
#run_sim(binary = T)





