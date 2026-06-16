source("non-lin-ln-k-iterative-3-pt-approx.R")
library(dplyr)
library(delaydiscount)
library(nlme)

# Input parameters, can be changed

npg <- 140  # Number of subjects per group
nstudies <- 1  # Number of studies to simulate

# Mean parameters for the two groups
EFT_mean = seq(from = -5.5, to = -7, by = -0.1)
NCC_mean = seq(from = -5.5, to = -7, by = -0.1)

result <- data.frame(matrix(nrow = 0, ncol = 6))
names(result) <- c("subj", "group", "true_ln_k", "nl_ln_k", "lin_ln_k", "nli_ln_k")

for(i in 1:length(EFT_mean)){
  # Set parameters
  conditions = c("EFT", "NCC")
  num_subj = c(npg, npg)
  time_points = c(30, 90, 180, 365, 1095, 1825, 3650)
  # mean_ln_k = c(-6, -7)
  mean_ln_k = c(EFT_mean[i], NCC_mean[i])  # vary EFT mean from -5.5 to -7, spaced by 0.1
  sigma_sq = 2
  g = 10.4
  
  # Set up output
  
  simulated_data <- simulate_dataset(conditions, num_subj*nstudies, time_points, 
                                     mean_ln_k, sigma_sq, g)
  simulated_data_p <- prepare_data_frame(simulated_data)
  simulated_data_p$study = 1
  
  # Do nonlinear analysis
  simulated_nl_ln_k = iter_ln_k_3_pt_approx_w_init(simulated_data_p,  
                                                   init_ln_k=simulated_data_p$lin_ln_k[seq(from = 7, to = length(simulated_data_p$lin_ln_k), by = 7)],
                                                   default_delta = 1)
  
  nli_ln_k <- simulated_nl_ln_k$nonlin_ln_k
  
  simulated_data_p$subj = as.factor(simulated_data_p$subj)
  simulated_data_p$group = as.factor(simulated_data_p$group)
  
  # Get true ln_k df
  true_ln_k_df <- simulated_data %>%
    filter(delay == 30) %>%
    select(subj, group, true_ln_k)
  
  # Fit nonlinear hierarchical model
  rFullFit<-nlme(indiff~1/(1+exp(logk)*delay), fixed=logk~group, random=logk~1|subj, data=simulated_data_p, start=c(-4, -4))
  
  # Get the est ln_ks for subjects
  nl_group_est <- cumsum(rFullFit$coefficients$fixed)
  names(nl_group_est) <- c("EFT", "NCC")
  nlm_ln_k <- rFullFit$coefficients$random[[1]] + nl_group_est[true_ln_k_df$group]
  
  # Fit linearized model
  # Get est ln_ks for subjects
  lin_est_ln_k <- get_subj_est_ln_k(simulated_data_p)
  lin_ln_k <- lin_est_ln_k$ln_k
  
  true_ln_k_df$nlm_ln_k <- as.numeric(nlm_ln_k)
  true_ln_k_df$lin_ln_k <- lin_ln_k
  true_ln_k_df$nli_ln_k <- nli_ln_k
  result <- rbind(result, true_ln_k_df)
}

write.csv(result, "sim_subj_ln_k_ests.csv", row.names = FALSE)