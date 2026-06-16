source("non-lin-ln-k-iterative-3-pt-approx.R")
library(dplyr)
library(delaydiscount)
library(nlme)


out_path <- "Simulation Results/"  # Folder for simulation results

# Input parameters, can be changed

npg <- 70  # Number of subjects per group
nstudies <- 1000  # Number of studies to simulate

# Mean parameters for the two groups
EFT_mean = seq(from = -5.5, to = -7, by = -0.1)
NCC_mean = seq(from = -5.5, to = -7, by = -0.1)

# For a power simulation, use different mean parameters
# EFT_mean = seq(from = -5.5, to = -7, by = -0.1)
# NCC_mean = rep(-6.5, 16)

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
  sim_res <- data.frame(matrix(nrow = nstudies, ncol = 33))
  names(sim_res) = c("true_eft_mean",	"true_ncc_mean", "sigma_sq", "g",
                     "n_eft",	"n_ncc",
                     "nl_est_eft",	"nl_est_ncc",	"p_val_nl",	"nl_mse_eft",	"nl_mse_ncc",
                     "nl_est_sd", "nl_est_re_sd",
                     "lin_est_eft",	"lin_est_ncc",	"p_val_lin", "lin_mse_eft",	"lin_mse_ncc",
                     "lin_est_sigma_sq", "lin_est_g",
                     "nli_est_eft", "nli_est_ncc", "p_val_nli", "nli_mse_eft", "nli_mse_ncc",
                     "lin_rc_est_eft", "lin_rc_est_ncc", "p_val_lin_rc",
                     "nli_rc_est_eft", "nli_rc_est_ncc", "p_val_nli_rc",
                     "n_eft_rc", "n_ncc_rc"
                      )
  sim_res$true_eft_mean <- mean_ln_k[1]
  sim_res$true_ncc_mean <- mean_ln_k[2]
  sim_res$sigma_sq <- sigma_sq
  sim_res$g <- g
  sim_res$n_eft <- num_subj[1]
  sim_res$n_ncc <- num_subj[2]
  
  simulated_data <- simulate_dataset(conditions, num_subj*nstudies, time_points, 
                                     mean_ln_k, sigma_sq, g)
  simulated_data_p <- prepare_data_frame(simulated_data)
  
  # Separate subjects by study
  simulated_data_p$study <- ceiling((simulated_data_p$subj %% (npg*nstudies))/npg)

  # Do nonlinear analysis
  simulated_nl_ln_k = iter_ln_k_3_pt_approx_w_init(simulated_data_p,  
                                                   init_ln_k=simulated_data_p$lin_ln_k[seq(from = 7, to = length(simulated_data_p$lin_ln_k), by = 7)],
                                                   default_delta = 1)
  
  simulated_data_p$nli_ln_k <- rep(simulated_nl_ln_k$nonlin_ln_k, each = 7)
  
  for(sim in 1:nstudies){
    single_study <- simulated_data_p %>%
      filter(study == sim)
    
    single_study$subj = as.factor(single_study$subj)
    single_study$group = as.factor(single_study$group)
    
    # Get true ln_k df
    true_ln_k_df <- simulated_data %>%
      filter(study == sim, delay == 30) %>%
      select(subj, group, true_ln_k)
    
    # Fit nonlinear hierarchical model
    rFullFit<-nlme(indiff~1/(1+exp(logk)*delay), fixed=logk~group, random=logk~1|subj, data=single_study, start=c(-4, -4))
    # Save parameter estimates
    sim_res$nl_est_eft[sim] <- rFullFit$coefficients$fixed[1]
    sim_res$nl_est_ncc[sim] <- sum(rFullFit$coefficients$fixed)
    sim_res$nl_est_sd[sim] <- rFullFit$sigma
    sim_res$nl_est_re_sd[sim] <- as.numeric(VarCorr(rFullFit)[1,2])
    
    # Save p-value
    sim_res$p_val_nl[sim] <- summary(rFullFit)$tTable[2,5]
    # Calculate MSEs
    #  Get the est ln_ks for subjects
    nl_group_est <- cumsum(rFullFit$coefficients$fixed)
    names(nl_group_est) <- c("EFT", "NCC")
    nl_com_est_ln_k <- rFullFit$coefficients$random[[1]] + nl_group_est[true_ln_k_df$group]
    nl_sq_err <- (nl_com_est_ln_k - true_ln_k_df$true_ln_k)^2
    sim_res$nl_mse_eft[sim] <- mean(nl_sq_err[true_ln_k_df$group == "EFT"])
    sim_res$nl_mse_ncc[sim] <- mean(nl_sq_err[true_ln_k_df$group == "NCC"])
    
    # Fit linearized model
    lin_model <- dd_hyperbolic_model(single_study)
    # Save parameter estimates
    sim_res$lin_est_eft[sim] <- lin_model$ln_k_mean$ln_k_mean[1]
    sim_res$lin_est_ncc[sim] <- lin_model$ln_k_mean$ln_k_mean[2]
    sim_res$lin_est_sigma_sq[sim] <- lin_model$var[1]
    sim_res$lin_est_g[sim] <- lin_model$var[2]
    # Save p-value
    sim_res$p_val_lin[sim] <- lin_model$pairwise_f_tests$p_value
    # Get est ln_ks for subjects
    lin_est_ln_k <- get_subj_est_ln_k(single_study)
    lin_sq_err <- (lin_est_ln_k$ln_k - true_ln_k_df$true_ln_k)^2
    sim_res$lin_mse_eft[sim] <- mean(lin_sq_err[true_ln_k_df$group == "EFT"])
    sim_res$lin_mse_ncc[sim] <- mean(lin_sq_err[true_ln_k_df$group == "NCC"])
    
    # Analyze nonlinear individual model
    nli_ln_k_df <- single_study %>%
      filter(delay == 30)
    sim_res$nli_est_eft[sim] <- mean(nli_ln_k_df$nli_ln_k[nli_ln_k_df$group == "EFT"])
    sim_res$nli_est_ncc[sim] <- mean(nli_ln_k_df$nli_ln_k[nli_ln_k_df$group == "NCC"])
    sim_res$p_val_nli[sim] <- t.test(nli_ln_k_df$nli_ln_k[nli_ln_k_df$group == "EFT"],
                                     nli_ln_k_df$nli_ln_k[nli_ln_k_df$group == "NCC"],
                                     var.equal = TRUE)$p.value
    lin_sq_err <- (lin_est_ln_k$ln_k - true_ln_k_df$true_ln_k)^2
    sim_res$nli_mse_eft[sim] <- mean((true_ln_k_df$true_ln_k[true_ln_k_df$group == "EFT"] -
                                       nli_ln_k_df$nli_ln_k[nli_ln_k_df$group == "EFT"])^2)
    sim_res$nli_mse_ncc[sim] <- mean((true_ln_k_df$true_ln_k[true_ln_k_df$group == "NCC"] -
                                        nli_ln_k_df$nli_ln_k[nli_ln_k_df$group == "NCC"])^2)
    
    # Apply rule check
    rc_res <- jb_rule_check(single_study)
    single_study_rc <- merge(single_study, rc_res) %>%
      filter(C1, C2)
    
    # Fit linearized model
    lin_rc_model <- dd_hyperbolic_model(single_study_rc)
    
    # Save parameter estimates
    sim_res$lin_rc_est_eft[sim] <- lin_rc_model$ln_k_mean$ln_k_mean[1]
    sim_res$lin_rc_est_ncc[sim] <- lin_rc_model$ln_k_mean$ln_k_mean[2]
    # Save p-value
    sim_res$p_val_lin_rc[sim] <- lin_rc_model$pairwise_f_tests$p_value
    
    # Save non-linear individual parameter estimates
    nli_rc_ln_k_df <- single_study_rc %>%
      filter(delay == 30)
    rc_eft <- nli_rc_ln_k_df$nli_ln_k[nli_rc_ln_k_df$group == "EFT"]
    rc_ncc <- nli_rc_ln_k_df$nli_ln_k[nli_rc_ln_k_df$group == "NCC"]
    rc_t <- t.test(rc_eft, rc_ncc, var.equal = TRUE)
    sim_res$nli_rc_est_eft[sim] <- rc_t$estimate[1]
    sim_res$nli_rc_est_ncc[sim] <- rc_t$estimate[2]
    sim_res$p_val_nli_rc[sim] <- rc_t$p.value
      
    # Save rule check sample size
    sim_res$n_eft_rc[sim] <- length(rc_eft)
    sim_res$n_ncc_rc[sim] <- length(rc_ncc)
  }
  
  
  # Save simulation
  out_name <- paste0("eft", EFT_mean[i]*10, "ncc", NCC_mean[i]*10,
                     "n_eft", num_subj[1], "n_ncc", npg)
  out_files <- list.files(path = out_path, pattern = out_name)
  if(length(out_files) > 0){
    out_name <- paste0(out_name, "_", length(out_files))
  }
  write.csv(sim_res, paste0(out_path, out_name, ".csv"))
}

