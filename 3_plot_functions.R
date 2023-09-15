library(tidyverse)

summarize_binary <- function(results){
  results %>% filter(binary) %>% 
    mutate(across(contains("clust"), as.factor)) %>% 
    mutate(across(contains("txt_eff"), as.factor)) %>% 
    group_by(binary,
             trt_clusts,
             trt_clust_size,
             control_clusts,
             control_clust_size,
             txt_eff,
             v, true_marg_mean_diff,
             sim_type, nsims) %>%
    summarize(tmle_powerTIE_or = mean(ltmle_or_lo > 1 | ltmle_or_hi < 1, na.rm = F),
              tmle_powerTIE_rr = mean(ltmle_rr_lo > 1 | ltmle_rr_hi < 1, na.rm = F),
              tmle_powerTIE_rd = mean(ltmle_lo > 0 | ltmle_hi < 0, na.rm = F),
              glmm_powerTIE = mean(glmm_lo > 1 | glmm_hi < 1, na.rm = T),
              gee_powerTIE = round(mean(gee_lo > 1 | gee_hi < 1, na.rm = T), digits = 4),
              tmle_cov_or = mean(ltmle_or_lo < tmle_true_marg_or & ltmle_or_hi > tmle_true_marg_or, na.rm = F),
              tmle_cov_rr = mean(ltmle_rr_lo < tmle_true_marg_rr & ltmle_rr_hi > tmle_true_marg_rr, na.rm = F),
              tmle_cov_rd = mean(ltmle_lo < tmle_true_marg_rd & ltmle_hi > tmle_true_marg_rd, na.rm = F),
              glmm_cov = mean(glmm_lo < glmm_true_cond_or & glmm_hi > glmm_true_cond_or, na.rm = T),
              gee_cov = round(mean(gee_lo < gee_true_marg_or & gee_hi > gee_true_marg_or, na.rm = T), digits = 4),
              tmle_bias_or = mean(log(ltmle_or_pe) - log(tmle_true_marg_or), na.rm = F),
              tmle_bias_rr = mean(log(ltmle_rr_pe) - log(tmle_true_marg_rr), na.rm = F),
              tmle_bias_rd = mean(ltmle_pe - tmle_true_marg_rd, na.rm = F),
              glmm_bias = mean(plogis(glmm_pe) - plogis(glmm_true_cond_or), na.rm = T),
              gee_bias = round(mean(plogis(gee_pe) - plogis(gee_true_marg_or), na.rm = T), digits = 4)) %>% 
    ungroup()
}

summarize_continuous <- function(results){
  results %>% filter(!binary) %>%
    mutate(across(contains("clust"), as.factor)) %>% 
    mutate(across(contains("txt_eff"), as.factor)) %>% 
    group_by(binary,
             trt_clusts,
             trt_clust_size,
             control_clusts,
             control_clust_size,
             txt_eff,
             v, true_marg_mean_diff,
             sim_type, nsims) %>%
    summarize(tmle_powerTIE = mean(ltmle_lo > 0 | ltmle_hi < 0),
              glmm_powerTIE = mean(glmm_lo > 0 | glmm_hi < 0),
              gee_powerTIE = mean(gee_lo > 0 | gee_hi < 0),
              tmle_cov = mean(ltmle_lo < true_marg_mean_diff & ltmle_hi > true_marg_mean_diff),
              glmm_cov = mean(glmm_lo < true_marg_mean_diff & glmm_hi > true_marg_mean_diff),
              gee_cov = mean(gee_lo < true_marg_mean_diff & gee_hi > true_marg_mean_diff),
              tmle_bias = mean(ltmle_pe - true_marg_mean_diff),
              glmm_bias = mean(glmm_pe - true_marg_mean_diff),
              gee_bias = mean(gee_pe - true_marg_mean_diff)) %>% 
    ungroup() 
}

##########################################
# TIE
##########################################
generate_TIE_plot <- function(dat, dodge_width = .4){
  
  if(dat$binary[1] == T){
    dodge_width <- .6
  } else {
    dodge_width <- .4
  } 
  
  plotdat_TIE <- dat %>% filter(txt_eff == 0) %>% 
    pivot_longer(cols = contains("TIE"), names_to = "method", values_to = "TIE") %>%
    mutate(id = row_number()) %>% group_by(id) %>% 
    mutate(TIE_lo = prop.test(x = TIE*nsims, n = nsims)$conf.int[1],
           TIE_hi = prop.test(x = TIE*nsims, n = nsims)$conf.int[2]) %>%
    ungroup() %>% select(-id)
  
  ggplot(data = plotdat_TIE,
         aes(x = interaction(control_clusts, trt_clusts),
             shape = txt_eff,
             color = method, group = interaction(method, txt_eff))) +
    geom_hline(aes(yintercept = .05), linetype = 2, alpha = .5) +
    geom_point(aes(y = TIE), position = position_dodge(width = dodge_width)) +
    geom_errorbar(aes(ymin = TIE_lo, ymax = TIE_hi, linetype = txt_eff), width = .2, position = position_dodge(width = dodge_width)) +
    facet_grid(control_clust_size ~ trt_clust_size, scales = "free_y", labeller = "label_both")
}
############################################ CONT



##########################################
# coverage
##########################################
generate_coverage_plot <- function(dat, dodge_width = .4){
  
  if(dat$binary[1] == T){
    dodge_width <- .6
  } else {
    dodge_width <- .4
  } 
  
  plotdat_cov <- dat %>% filter(txt_eff != 0) %>% 
    pivot_longer(cols = contains("cov"), names_to = "method", values_to = "coverage") %>%
    mutate(id = row_number()) %>% group_by(id) %>% 
    mutate(cov_lo = prop.test(x = coverage*nsims, n = nsims)$conf.int[1],
           cov_hi = prop.test(x = coverage*nsims, n = nsims)$conf.int[2]) %>%
    ungroup() %>% select(-id)
  
  ggplot(data = plotdat_cov,
         aes(x = interaction(control_clusts, trt_clusts),
             shape = txt_eff,
             color = method, group = interaction(method, txt_eff))) +
    geom_hline(aes(yintercept = .95), linetype = 2, alpha = .5) +
    geom_point(aes(y = coverage), position = position_dodge(width = dodge_width)) +
    geom_errorbar(aes(ymin = cov_lo, ymax = cov_hi, linetype = txt_eff), width = .2, position = position_dodge(width = dodge_width)) +
    facet_grid(control_clust_size ~ trt_clust_size, scales = "free_y", labeller = "label_both")
}


##########################################
# bias
##########################################
generate_bias_plot <- function(dat, dodge_width = .4){
  
  if(dat$binary[1] == T){
    dodge_width <- .6
  } else {
    dodge_width <- .4
  } 
  
  plotdat_bias <- dat %>% filter(txt_eff != 0) %>% 
    pivot_longer(cols = contains("bias"), names_to = "method", values_to = "bias")
  
  ggplot(data = plotdat_bias,
         aes(x = interaction(control_clusts, trt_clusts),
             shape = txt_eff,
             color = method, group = interaction(method, txt_eff))) +
    geom_hline(aes(yintercept = 0), linetype = 2, alpha = .5) +
    geom_point(aes(y = bias), position = position_dodge(width = dodge_width)) +
    facet_grid(control_clust_size ~ trt_clust_size, scales = "free_y", labeller = "label_both")
}


##########################################
# power
##########################################
generate_power_plot <- function(dat, dodge_width = .4){
  
  if(dat$binary[1] == T){
    dodge_width <- .6
  } else {
    dodge_width <- .4
  } 
  
  plotdat_power <- dat %>% filter(txt_eff != 0) %>% 
    pivot_longer(cols = contains("power"), names_to = "method", values_to = "power") %>%
    mutate(id = row_number()) %>% group_by(id) %>% 
    mutate(power_lo = prop.test(x = power*nsims, n = nsims)$conf.int[1],
           power_hi = prop.test(x = power*nsims, n = nsims)$conf.int[2]) %>%
    ungroup() %>% select(-id)
  
  ggplot(data = plotdat_power,
         aes(x = interaction(control_clusts, trt_clusts),
             shape = txt_eff,
             color = method, group = interaction(method, txt_eff))) +
    geom_point(aes(y = power), position = position_dodge(width = dodge_width)) +
    geom_errorbar(aes(ymin = power_lo, ymax = power_hi, linetype = txt_eff),
                  width = .3, position = position_dodge(width = dodge_width)) +
    facet_grid(control_clust_size ~ trt_clust_size, scales = "free_y", labeller = "label_both")
}


























