library(dplyr)

path <- "Simulation Results/"

sim_files <- list.files(path)

summary_df <- data.frame(matrix(nrow = length(sim_files), ncol = 32))
names(summary_df) <- c("sims", "n_eft", "n_ncc",
                       "true_eft_mean", "true_ncc_mean", "sigma_sq", "g",
                       "nl_power", "nl_mean_est_eft_param", "nl_mean_est_ncc_param", "nl_mse_eft_param", "nl_mse_ncc_param", "nl_mse_eft_subj", "nl_mse_ncc_subj",
                       "lin_power", "lin_mean_est_eft_param", "lin_mean_est_ncc_param", "lin_mse_eft_param", "lin_mse_ncc_param", "lin_mse_eft_subj", "lin_mse_ncc_subj",
                       "nli_power", "nli_mean_est_eft_param", "nli_mean_est_ncc_param", "nli_mse_eft_param", "nli_mse_ncc_param", "nli_mse_eft_subj", "nli_mse_ncc_subj",
                       "lin_rc_power", "nli_rc_power", "prop_rc_lost_eft", "prop_rc_lost_ncc")

for(i in 1:length(sim_files)){
  sim_res <- read.csv(paste0(path, sim_files[i]))
  
  # Save number of sims
  summary_df$sims[i] <- length(sim_res$true_eft_mean)
  
  summary_df$n_eft[i] <- sim_res$n_eft[1]
  summary_df$n_ncc[i] <- sim_res$n_ncc[1]
  
  # Save the parameters
  summary_df$true_eft_mean[i] <- sim_res$true_eft_mean[1]
  summary_df$true_ncc_mean[i] <- sim_res$true_ncc_mean[1]
  summary_df$sigma_sq[i] <- sim_res$sigma_sq[1]
  summary_df$g[i] <- sim_res$g[1]
  
  # Get the power
  summary_df$nl_power[i] <- mean(sim_res$p_val_nl <= 0.05)
  summary_df$lin_power[i] <- mean(sim_res$p_val_lin <= 0.05)
  summary_df$nli_power[i] <- mean(sim_res$p_val_nli <= 0.05)
  summary_df$lin_rc_power[i] <- mean(sim_res$p_val_lin_rc <= 0.05)
  summary_df$nli_rc_power[i] <- mean(sim_res$p_val_nli_rc <= 0.05)
  
  # Save the mean of the parameter estimates
  summary_df$nl_mean_est_eft_param[i] <- mean(sim_res$nl_est_eft)
  summary_df$nl_mean_est_ncc_param[i] <- mean(sim_res$nl_est_ncc)
  summary_df$lin_mean_est_eft_param[i] <- mean(sim_res$lin_est_eft)
  summary_df$lin_mean_est_ncc_param[i] <- mean(sim_res$lin_est_ncc)
  summary_df$nli_mean_est_eft_param[i] <- mean(sim_res$nli_est_eft)
  summary_df$nli_mean_est_ncc_param[i] <- mean(sim_res$nli_est_ncc)
  
  # Get the parameter MSE
  summary_df$nl_mse_eft_param[i] <- mean((sim_res$nl_est_eft - sim_res$true_eft_mean)^2)
  summary_df$nl_mse_ncc_param[i] <- mean((sim_res$nl_est_ncc - sim_res$true_ncc_mean)^2)
  summary_df$lin_mse_eft_param[i] <- mean((sim_res$lin_est_eft - sim_res$true_eft_mean)^2)
  summary_df$lin_mse_ncc_param[i] <- mean((sim_res$lin_est_ncc - sim_res$true_ncc_mean)^2)
  summary_df$nli_mse_eft_param[i] <- mean((sim_res$nli_est_eft - sim_res$true_eft_mean)^2)
  summary_df$nli_mse_ncc_param[i] <- mean((sim_res$nli_est_ncc - sim_res$true_ncc_mean)^2)
  
  
  # Get mean subject MSEs
  summary_df$nl_mse_eft_subj[i] <- mean(sim_res$nl_mse_eft)
  summary_df$nl_mse_ncc_subj[i] <- mean(sim_res$nl_mse_ncc)
  summary_df$lin_mse_eft_subj[i] <- mean(sim_res$lin_mse_eft)
  summary_df$lin_mse_ncc_subj[i] <- mean(sim_res$lin_mse_ncc)
  summary_df$nli_mse_eft_subj[i] <- mean(sim_res$nli_mse_eft)
  summary_df$nli_mse_ncc_subj[i] <- mean(sim_res$nli_mse_ncc)
  
  # Get proportion of subjects lost due to rule check
  summary_df$prop_rc_lost_eft[i] <- 1-mean(sim_res$n_eft_rc)/mean(sim_res$n_eft)
  summary_df$prop_rc_lost_ncc[i] <- 1-mean(sim_res$n_eft_rc)/mean(sim_res$n_ncc)
}

# Merge observations from identical cases
summary_df <- summary_df %>%
  group_by(n_eft, n_ncc, true_eft_mean, true_ncc_mean, sigma_sq, g) %>%
  summarise(nl_power = sum(sims*nl_power)/sum(sims),
            nl_mean_est_eft_param = sum(sims*nl_mean_est_eft_param)/sum(sims),
            nl_mean_est_ncc_param = sum(sims*nl_mean_est_ncc_param)/sum(sims),
            nl_mse_eft_param = sum(sims*nl_mse_eft_param)/sum(sims),
            nl_mse_ncc_param = sum(sims*nl_mse_ncc_param)/sum(sims),
            nl_mse_eft_subj = sum(sims*nl_mse_eft_subj)/sum(sims),
            nl_mse_ncc_subj = sum(sims*nl_mse_ncc_subj)/sum(sims),
            lin_power = sum(sims*lin_power)/sum(sims),
            lin_mean_est_eft_param = sum(sims*lin_mean_est_eft_param)/sum(sims),
            lin_mean_est_ncc_param = sum(sims*lin_mean_est_ncc_param)/sum(sims),
            lin_mse_eft_param = sum(sims*lin_mse_eft_param)/sum(sims),
            lin_mse_ncc_param = sum(sims*lin_mse_ncc_param)/sum(sims),
            lin_mse_eft_subj = sum(sims*lin_mse_eft_subj)/sum(sims),
            lin_mse_ncc_subj = sum(sims*lin_mse_ncc_subj)/sum(sims),
            nli_power = sum(sims*nli_power)/sum(sims),
            nli_mean_est_eft_param = sum(sims*nli_mean_est_eft_param)/sum(sims),
            nli_mean_est_ncc_param = sum(sims*nli_mean_est_ncc_param)/sum(sims),
            nli_mse_eft_param = sum(sims*nli_mse_eft_param)/sum(sims),
            nli_mse_ncc_param = sum(sims*nli_mse_ncc_param)/sum(sims),
            nli_mse_eft_subj = sum(sims*nli_mse_eft_subj)/sum(sims),
            nli_mse_ncc_subj = sum(sims*nli_mse_ncc_subj)/sum(sims),
            lin_rc_power = sum(sims*lin_rc_power)/sum(sims),
            nli_rc_power = sum(sims*nli_rc_power)/sum(sims),
            prop_rc_lost_eft = sum(sims*prop_rc_lost_eft)/sum(sims),
            prop_rc_lost_ncc = sum(sims*prop_rc_lost_ncc)/sum(sims),
            sims = sum(sims)
  ) %>%
  mutate(d = sqrt(7)*(true_ncc_mean - true_eft_mean)/sqrt(sigma_sq*(g+1)))

sdf140_t1e <- summary_df %>%
  filter(n_eft == 140, true_eft_mean == true_ncc_mean) %>%
  arrange(true_eft_mean)

sdf140_pow <- summary_df %>%
  filter(n_eft == 140, true_ncc_mean == -6.5) %>%
  arrange(true_eft_mean)

sdf140_pow_55 <- summary_df %>%
  filter(n_eft == 140, true_ncc_mean == -5.5) %>%
  arrange(true_eft_mean)

sdf35_pow <- summary_df %>%
  filter(n_eft == 35, true_ncc_mean == -6.5) %>%
  arrange(true_eft_mean)

sdf35_pow_55 <- summary_df %>%
  filter(n_eft == 35, true_ncc_mean == -5.5) %>%
  arrange(true_eft_mean)

sdf35_t1e <- summary_df %>%
  filter(n_eft == 35, true_eft_mean == true_ncc_mean) %>%
  arrange(true_eft_mean)

sdf70_pow <- summary_df %>%
  filter(n_eft == 70, true_ncc_mean == -6.5) %>%
  arrange(true_eft_mean)

sdf70_pow_55 <- summary_df %>%
  filter(n_eft == 70, true_ncc_mean == -5.5) %>%
  arrange(true_eft_mean)

sdf70_t1e <- summary_df %>%
  filter(n_eft == 70, true_eft_mean == true_ncc_mean) %>%
  arrange(true_eft_mean)


##############################
# Figure 3 Code
##############################

# Plot mean square errors for group hyperparameter and subject parameter estimates

# N=35
# Hyperparameter Estimate
plot(sdf35_t1e$true_eft_mean, (sdf35_t1e$nl_mse_eft_param+sdf35_t1e$nl_mse_ncc_param)/2, type = "l", col = "red", main = "Aggregate Param Mean Square Error", sub = "EFT = NCC", ylim = c(0.07, 0.1))
lines(sdf35_t1e$true_eft_mean, (sdf35_t1e$lin_mse_eft_param+sdf35_t1e$lin_mse_ncc_param)/2, col = "blue")
lines(sdf35_t1e$true_eft_mean, (sdf35_t1e$nli_mse_eft_param+sdf35_t1e$nli_mse_ncc_param)/2, col = "orange")
legend(x = "bottomleft", 
       legend = c("Young", "Lin. Hier.", "Nonlin."),
       col = c("red", "blue", "orange"), lwd = 2)

# Subject Parameter Estimates
plot(sdf35_t1e$true_eft_mean, (sdf35_t1e$nl_mse_eft_subj+sdf35_t1e$nl_mse_ncc_subj)/2, type = "l", col = "red", main = "Aggregate Subj Mean Square Error", sub = "EFT = NCC", ylim = c(0, 0.7))
lines(sdf35_t1e$true_eft_mean, (sdf35_t1e$lin_mse_eft_subj+sdf35_t1e$lin_mse_ncc_subj)/2, col = "blue")
lines(sdf35_t1e$true_eft_mean, (sdf35_t1e$nli_mse_eft_subj+sdf35_t1e$nli_mse_ncc_subj)/2, col = "orange")
legend(x = "bottomleft", 
       legend = c("Young", "Lin. Hier.", "Nonlin."),
       col = c("red", "blue", "orange"), lwd = 2)

# N=70
# Hyperparameter Estimate
plot(sdf70_t1e$true_eft_mean, (sdf70_t1e$nl_mse_eft_param+sdf70_t1e$nl_mse_ncc_param)/2, type = "l", col = "red", main = "Aggregate Param Mean Square Error", sub = "EFT = NCC")
lines(sdf70_t1e$true_eft_mean, (sdf70_t1e$lin_mse_eft_param+sdf70_t1e$lin_mse_ncc_param)/2, col = "blue")
lines(sdf70_t1e$true_eft_mean, (sdf70_t1e$nli_mse_eft_param+sdf70_t1e$nli_mse_ncc_param)/2, col = "orange")

# Subject Parameter Estimates
plot(sdf70_t1e$true_eft_mean, (sdf70_t1e$nl_mse_eft_subj+sdf70_t1e$nl_mse_ncc_subj)/2, type = "l", col = "red", main = "Aggregate Subj Mean Square Error", sub = "EFT = NCC", ylim = c(0, 0.7))
lines(sdf70_t1e$true_eft_mean, (sdf70_t1e$lin_mse_eft_subj+sdf70_t1e$lin_mse_ncc_subj)/2, col = "blue")
lines(sdf70_t1e$true_eft_mean, (sdf70_t1e$nli_mse_eft_subj+sdf70_t1e$nli_mse_ncc_subj)/2, col = "orange")
legend(x = "bottomleft", 
       legend = c("Young", "Lin. Hier.", "Nonlin."),
       col = c("red", "blue", "orange"), lwd = 2)

# N=140
# Hyperparameter Estimate
plot(sdf140_t1e$true_eft_mean, (sdf140_t1e$nl_mse_eft_param+sdf140_t1e$nl_mse_ncc_param)/2, type = "l", col = "red", main = "Aggregate Param Mean Square Error", sub = "EFT = NCC")
lines(sdf140_t1e$true_eft_mean, (sdf140_t1e$lin_mse_eft_param+sdf140_t1e$lin_mse_ncc_param)/2, col = "blue")
lines(sdf140_t1e$true_eft_mean, (sdf140_t1e$nli_mse_eft_param+sdf140_t1e$nli_mse_ncc_param)/2, col = "orange")

# Subject Parameter Estimates
plot(sdf140_t1e$true_eft_mean, (sdf140_t1e$nl_mse_eft_subj+sdf140_t1e$nl_mse_ncc_subj)/2, type = "l", col = "red", main = "Aggregate Subj Mean Square Error", sub = "EFT = NCC", ylim = c(0, 0.7))
lines(sdf140_t1e$true_eft_mean, (sdf140_t1e$lin_mse_eft_subj+sdf140_t1e$lin_mse_ncc_subj)/2, col = "blue")
lines(sdf140_t1e$true_eft_mean, (sdf140_t1e$nli_mse_eft_subj+sdf140_t1e$nli_mse_ncc_subj)/2, col = "orange")
legend(x = "bottomleft", 
       legend = c("Young", "Lin. Hier.", "Nonlin."),
       col = c("red", "blue", "orange"), lwd = 2)



##############################
# Figure 4 Code
##############################

# Produce scatterplots of individual subject true vs estimated ln(k) values,
#  with loess fit curves

big_true_ln_k_df <- read.csv("sim_subj_ln_k_ests.csv", row.names = F)

lin_loess <- loess(lin_ln_k ~ true_ln_k, data = big_true_ln_k_df)
nlm_loess <- loess(nlm_ln_k ~ true_ln_k, data = big_true_ln_k_df)
nli_loess <- loess(nli_ln_k ~ true_ln_k, data = big_true_ln_k_df)

plot(big_true_ln_k_df$true_ln_k, big_true_ln_k_df$lin_ln_k)
abline(a=0, b=1, col = "black")
lines((-120:0)/10, predict(lin_loess, (-120:0)/10), col = "blue")

plot(big_true_ln_k_df$true_ln_k, big_true_ln_k_df$nlm_ln_k)
abline(a=0, b=1, col = "black")
lines((-120:0)/10, predict(nlm_loess, (-120:0)/10), col = "blue")

plot(big_true_ln_k_df$true_ln_k, big_true_ln_k_df$nli_ln_k)
abline(a=0, b=1, col = "black")
lines((-120:0)/10, predict(nli_loess, (-120:0)/10), col = "blue")


##############################
# Figure 5 Code
##############################

# Plot power differences between methods

# Power difference plot (N=35)
plot(sdf35_pow_55$true_eft_mean, sdf35_pow_55$lin_power-sdf35_pow_55$nl_power, main = "Power improvement of linear test w/ no rule check", type = "l", ylim = c(-0.007, 0.15), col = "red", lwd=2, sub = "N=35, NCC_mean=-5.5", ylab = "Power Improvement of lin. test")
lines(sdf35_pow_55$true_eft_mean, sdf35_pow_55$lin_power-sdf35_pow_55$nli_power, col = "orange", lwd=2)
lines(sdf35_pow_55$true_eft_mean, sdf35_pow_55$lin_power-sdf35_pow_55$lin_rc_power, col = "purple", lwd=2)
lines(sdf35_pow_55$true_eft_mean, sdf35_pow_55$lin_power-sdf35_pow_55$nli_rc_power, col = "darkgreen", lwd=2)
abline(h=0)
legend(x = "topright", 
       legend = c("Young", "Nonlin.", "Lin. Hier. w/RC", "Nonlin. w/RC"),
       col = c("red", "orange", "purple", "darkgreen"), lwd = 2)

# Power difference plot (N=70)
plot(sdf70_pow_55$true_eft_mean, sdf70_pow_55$lin_power-sdf70_pow_55$nl_power, main = "Power improvement of linear test w/ no rule check", type = "l", ylim = c(-0.007, 0.15), col = "red", lwd=2, sub = "N=70, NCC_mean=-5.5", ylab = "Power Improvement of lin. test")
lines(sdf70_pow_55$true_eft_mean, sdf70_pow_55$lin_power-sdf70_pow_55$nli_power, col = "orange", lwd=2)
lines(sdf70_pow_55$true_eft_mean, sdf70_pow_55$lin_power-sdf70_pow_55$lin_rc_power, col = "purple", lwd=2)
lines(sdf70_pow_55$true_eft_mean, sdf70_pow_55$lin_power-sdf70_pow_55$nli_rc_power, col = "darkgreen", lwd=2)
abline(h=0)
legend(x = "topright", 
       legend = c("Young", "Nonlin.", "Lin. Hier. w/RC", "Nonlin. w/RC"),
       col = c("red", "orange", "purple", "darkgreen"), lwd = 2)

# Power difference plot (N=140)
plot(sdf140_pow_55$true_eft_mean, sdf140_pow_55$lin_power-sdf140_pow_55$nl_power, main = "Power improvement of linear test w/ no rule check", type = "l", ylim = c(-0.007, 0.15), col = "red", lwd=2, sub = "N=140, NCC_mean=-5.5", ylab = "Power Improvement of lin. test")
lines(sdf140_pow_55$true_eft_mean, sdf140_pow_55$lin_power-sdf140_pow_55$nli_power, col = "orange", lwd=2)
lines(sdf140_pow_55$true_eft_mean, sdf140_pow_55$lin_power-sdf140_pow_55$lin_rc_power, col = "purple", lwd=2)
lines(sdf140_pow_55$true_eft_mean, sdf140_pow_55$lin_power-sdf140_pow_55$nli_rc_power, col = "darkgreen", lwd=2)
abline(h=0)
legend(x = "topleft", 
       legend = c("Young", "Nonlin.", "Lin. Hier. w/RC", "Nonlin. w/RC"),
       col = c("red", "orange", "purple", "darkgreen"), lwd = 2)



