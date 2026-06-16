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

# Type 1 Error Plot (N=140)

# Type 1 Error (N = 35)
plot(sdf35_t1e$true_eft_mean, sdf35_t1e$nl_power, main = "Type 1 Error by ln(k) mean", type = "l", col = "red", ylim = c(0, 0.06), ylab = "Type 1 Error", xlab = "ln(k) mean",
     sub = "N_EFT = 35, N_NCC = 35", lwd = 2)
lines(sdf35_t1e$true_eft_mean, sdf35_t1e$lin_power, col = "blue", lwd = 2)
lines(sdf35_t1e$true_eft_mean, sdf35_t1e$nli_power, col = "orange", lwd=2)
lines(sdf35_t1e$true_eft_mean, sdf35_t1e$lin_rc_power, col = "purple",lwd=2)
lines(sdf35_t1e$true_eft_mean, sdf35_t1e$nli_rc_power, col = "darkgreen", lwd=2)
abline(h=0.05)
legend(x = "bottomleft", 
       legend = c("Young", "Lin. Hier.", "Nonlin.", "Lin. Hier. w/RC", "Nonlin. w/RC"),
       col = c("red", "blue", "orange", "purple", "darkgreen"), lwd = 2)

# Type 1 Error (N = 70)
plot(sdf70_t1e$true_eft_mean, sdf70_t1e$nl_power, main = "Type 1 Error by ln(k) mean", type = "l", col = "red", ylim = c(0, 0.06), ylab = "Type 1 Error", xlab = "ln(k) mean",
     sub = "N_EFT = 70, N_NCC = 70", lwd = 2)
lines(sdf70_t1e$true_eft_mean, sdf70_t1e$lin_power, col = "blue", lwd = 2)
lines(sdf70_t1e$true_eft_mean, sdf70_t1e$nli_power, col = "orange", lwd=2)
lines(sdf70_t1e$true_eft_mean, sdf70_t1e$lin_rc_power, col = "purple",lwd=2)
lines(sdf70_t1e$true_eft_mean, sdf70_t1e$nli_rc_power, col = "darkgreen", lwd=2)
abline(h=0.05)
legend(x = "bottomleft", 
       legend = c("Young", "Lin. Hier.", "Nonlin.", "Lin. Hier. w/RC", "Nonlin. w/RC"),
       col = c("red", "blue", "orange", "purple", "darkgreen"), lwd = 2)

# Type 1 Error (N = 140)
plot(sdf140_t1e$true_eft_mean, sdf140_t1e$nl_power, main = "Type 1 Error by ln(k) mean", type = "l", col = "red", ylim = c(0, 0.06), ylab = "Type 1 Error", xlab = "ln(k) mean",
     sub = "N_EFT = 140, N_NCC = 140", lwd = 2)
lines(sdf140_t1e$true_eft_mean, sdf140_t1e$lin_power, col = "blue", lwd = 2)
lines(sdf140_t1e$true_eft_mean, sdf140_t1e$nli_power, col = "orange", lwd=2)
lines(sdf140_t1e$true_eft_mean, sdf140_t1e$lin_rc_power, col = "purple",lwd=2)
lines(sdf140_t1e$true_eft_mean, sdf140_t1e$nli_rc_power, col = "darkgreen", lwd=2)
abline(h=0.05)
legend(x = "bottomleft", 
       legend = c("Young", "Lin. Hier.", "Nonlin.", "Lin. Hier. w/RC", "Nonlin. w/RC"),
       col = c("red", "blue", "orange", "purple", "darkgreen"), lwd = 2)


# Power Plot (N = 35)
plot(sdf35_pow$true_eft_mean, sdf35_pow$nl_power, main = "Power by ln(k) mean", type = "l", col = "red", ylim = c(0, 1), lwd=2, sub = "N=35, NCC_mean=-6.5")
lines(sdf35_pow$true_eft_mean, sdf35_pow$lin_power, col = "blue", lwd=2)
lines(sdf35_pow$true_eft_mean, sdf35_pow$nli_power, col = "orange", lwd=2)
lines(sdf35_pow$true_eft_mean, sdf35_pow$lin_rc_power, col = "purple", lwd=2)
lines(sdf35_pow$true_eft_mean, sdf35_pow$nli_rc_power, col = "darkgreen", lwd=2)

plot(sdf35_pow_55$true_eft_mean, sdf35_pow_55$nl_power, main = "Power by ln(k) mean", type = "l", col = "red", ylim = c(0, 1), lwd =2, sub = "N=35, NCC_mean=-5.5")
lines(sdf35_pow_55$true_eft_mean, sdf35_pow_55$lin_power, col = "blue", lwd=2)
lines(sdf35_pow_55$true_eft_mean, sdf35_pow_55$nli_power, col = "orange", lwd=2)
lines(sdf35_pow_55$true_eft_mean, sdf35_pow_55$lin_rc_power, col = "purple", lwd=2)
lines(sdf35_pow_55$true_eft_mean, sdf35_pow_55$nli_rc_power, col = "darkgreen", lwd=2)
legend(x = "topright", 
       legend = c("Young", "Lin. Hier.", "Nonlin.", "Lin. Hier. w/RC", "Nonlin. w/RC"),
       col = c("red", "blue", "orange", "purple", "darkgreen"), lwd = 2)

# Power difference plot
plot(sdf35_pow_55$true_eft_mean, sdf35_pow_55$lin_power-sdf35_pow_55$nl_power, main = "Power improvement of linear test w/ no rule check", type = "l", ylim = c(-0.007, 0.15), col = "red", lwd=2, sub = "N=35, NCC_mean=-5.5", ylab = "Power Improvement of lin. test")
lines(sdf35_pow_55$true_eft_mean, sdf35_pow_55$lin_power-sdf35_pow_55$nli_power, col = "orange", lwd=2)
lines(sdf35_pow_55$true_eft_mean, sdf35_pow_55$lin_power-sdf35_pow_55$lin_rc_power, col = "purple", lwd=2)
lines(sdf35_pow_55$true_eft_mean, sdf35_pow_55$lin_power-sdf35_pow_55$nli_rc_power, col = "darkgreen", lwd=2)
abline(h=0)
legend(x = "topright", 
       legend = c("Young", "Nonlin.", "Lin. Hier. w/RC", "Nonlin. w/RC"),
       col = c("red", "orange", "purple", "darkgreen"), lwd = 2)

# Power Plot (N = 70)
plot(sdf70_pow$true_eft_mean, sdf70_pow$nl_power, main = "Power by ln(k) mean", type = "l", col = "red", ylim = c(0, 1), lwd=2, sub = "N=70, NCC_mean=-6.5")
lines(sdf70_pow$true_eft_mean, sdf70_pow$lin_power, col = "blue", lwd=2)
lines(sdf70_pow$true_eft_mean, sdf70_pow$nli_power, col = "orange", lwd=2)
lines(sdf70_pow$true_eft_mean, sdf70_pow$lin_rc_power, col = "purple", lwd=2)
lines(sdf70_pow$true_eft_mean, sdf70_pow$nli_rc_power, col = "darkgreen", lwd=2)
legend(x = "topleft", 
       legend = c("Young", "Lin. Hier.", "Nonlin.", "Lin. Hier. w/RC", "Nonlin. w/RC"),
       col = c("red", "blue", "orange", "purple", "darkgreen"), lwd = 2)


plot(sdf70_pow_55$true_eft_mean, sdf70_pow_55$nl_power, main = "Power by ln(k) mean", type = "l", col = "red", ylim = c(0, 1), lwd=2, sub = "N=70, NCC_mean=-5.5")
lines(sdf70_pow_55$true_eft_mean, sdf70_pow_55$lin_power, col = "blue", lwd=2)
lines(sdf70_pow_55$true_eft_mean, sdf70_pow_55$nli_power, col = "orange", lwd=2)
lines(sdf70_pow_55$true_eft_mean, sdf70_pow_55$lin_rc_power, col = "purple", lwd=2)
lines(sdf70_pow_55$true_eft_mean, sdf70_pow_55$nli_rc_power, col = "darkgreen", lwd=2)
legend(x = "topright", 
       legend = c("Young", "Lin. Hier.", "Nonlin.", "Lin. Hier. w/RC", "Nonlin. w/RC"),
       col = c("red", "blue", "orange", "purple", "darkgreen"), lwd = 2)

# Power difference plot
plot(sdf70_pow_55$true_eft_mean, sdf70_pow_55$lin_power-sdf70_pow_55$nl_power, main = "Power improvement of linear test w/ no rule check", type = "l", ylim = c(-0.007, 0.15), col = "red", lwd=2, sub = "N=70, NCC_mean=-5.5", ylab = "Power Improvement of lin. test")
lines(sdf70_pow_55$true_eft_mean, sdf70_pow_55$lin_power-sdf70_pow_55$nli_power, col = "orange", lwd=2)
lines(sdf70_pow_55$true_eft_mean, sdf70_pow_55$lin_power-sdf70_pow_55$lin_rc_power, col = "purple", lwd=2)
lines(sdf70_pow_55$true_eft_mean, sdf70_pow_55$lin_power-sdf70_pow_55$nli_rc_power, col = "darkgreen", lwd=2)
abline(h=0)
legend(x = "topright", 
       legend = c("Young", "Nonlin.", "Lin. Hier. w/RC", "Nonlin. w/RC"),
       col = c("red", "orange", "purple", "darkgreen"), lwd = 2)

# Power Plot (N = 140)
plot(sdf140_pow$true_eft_mean, sdf140_pow$nl_power, main = "Power by ln(k) mean", type = "l", col = "red", ylim = c(0, 1), lwd=2, sub = "N=140, NCC_mean=-6.5")
lines(sdf140_pow$true_eft_mean, sdf140_pow$lin_power, col = "blue", lwd=2)
lines(sdf140_pow$true_eft_mean, sdf140_pow$nli_power, col = "orange", lwd=2)
lines(sdf140_pow$true_eft_mean, sdf140_pow$lin_rc_power, col = "purple", lwd=2)
lines(sdf140_pow$true_eft_mean, sdf140_pow$nli_rc_power, col = "darkgreen", lwd=2)
legend(x = "topleft", 
       legend = c("Young", "Lin. Hier.", "Nonlin.", "Lin. Hier. w/RC", "Nonlin. w/RC"),
       col = c("red", "blue", "orange", "purple", "darkgreen"), lwd = 2)


plot(sdf140_pow_55$true_eft_mean, sdf140_pow_55$nl_power, main = "Power by ln(k) mean", type = "l", col = "red", ylim = c(0, 1), lwd=2, sub = "N=140, NCC_mean=-5.5")
lines(sdf140_pow_55$true_eft_mean, sdf140_pow_55$lin_power, col = "blue", lwd=2)
lines(sdf140_pow_55$true_eft_mean, sdf140_pow_55$nli_power, col = "orange", lwd=2)
lines(sdf140_pow_55$true_eft_mean, sdf140_pow_55$lin_rc_power, col = "purple", lwd=2)
lines(sdf140_pow_55$true_eft_mean, sdf140_pow_55$nli_rc_power, col = "darkgreen", lwd=2)
legend(x = "bottomleft", 
       legend = c("Young", "Lin. Hier.", "Nonlin.", "Lin. Hier. w/RC", "Nonlin. w/RC"),
       col = c("red", "blue", "orange", "purple", "darkgreen"), lwd = 2)

# Power difference plot
plot(sdf140_pow_55$true_eft_mean, sdf140_pow_55$lin_power-sdf140_pow_55$nl_power, main = "Power improvement of linear test w/ no rule check", type = "l", ylim = c(-0.007, 0.15), col = "red", lwd=2, sub = "N=140, NCC_mean=-5.5", ylab = "Power Improvement of lin. test")
lines(sdf140_pow_55$true_eft_mean, sdf140_pow_55$lin_power-sdf140_pow_55$nli_power, col = "orange", lwd=2)
lines(sdf140_pow_55$true_eft_mean, sdf140_pow_55$lin_power-sdf140_pow_55$lin_rc_power, col = "purple", lwd=2)
lines(sdf140_pow_55$true_eft_mean, sdf140_pow_55$lin_power-sdf140_pow_55$nli_rc_power, col = "darkgreen", lwd=2)
abline(h=0)
legend(x = "topleft", 
       legend = c("Young", "Nonlin.", "Lin. Hier. w/RC", "Nonlin. w/RC"),
       col = c("red", "orange", "purple", "darkgreen"), lwd = 2)

# Save datasets that makes power plots
pow_35 <- sdf35_pow_55 %>%
  select(n_eft, n_ncc, true_eft_mean, true_ncc_mean, sigma_sq, g, 
         nl_power, lin_power, nli_power, lin_rc_power, nli_rc_power) %>%
  mutate(d = sqrt(7)*abs(true_eft_mean - true_ncc_mean)/sqrt(sigma_sq*(g+1)))
# write.csv(pow_35, "Plots/pow_n35.csv", row.names = FALSE)

pow_70 <- sdf70_pow_55 %>%
  select(n_eft, n_ncc, true_eft_mean, true_ncc_mean, sigma_sq, g, 
         nl_power, lin_power, nli_power, lin_rc_power, nli_rc_power) %>%
  mutate(d = sqrt(7)*abs(true_eft_mean - true_ncc_mean)/sqrt(sigma_sq*(g+1)))
# write.csv(pow_70, "Plots/pow_n70.csv", row.names = FALSE)

pow_140 <- sdf140_pow_55 %>%
  select(n_eft, n_ncc, true_eft_mean, true_ncc_mean, sigma_sq, g, 
         nl_power, lin_power, nli_power, lin_rc_power, nli_rc_power) %>%
  mutate(d = sqrt(7)*abs(true_eft_mean - true_ncc_mean)/sqrt(sigma_sq*(g+1)))
# write.csv(pow_140, "Plots/pow_n140.csv", row.names = FALSE)



# MSE Plot

# Type 1 Error scenarios (N=35)
plot(sdf35_t1e$true_eft_mean, (sdf35_t1e$nl_mse_eft_param+sdf35_t1e$nl_mse_ncc_param)/2, type = "l", col = "red", main = "Aggregate Param Mean Square Error", sub = "EFT = NCC", ylim = c(0.07, 0.1))
lines(sdf35_t1e$true_eft_mean, (sdf35_t1e$lin_mse_eft_param+sdf35_t1e$lin_mse_ncc_param)/2, col = "blue")
lines(sdf35_t1e$true_eft_mean, (sdf35_t1e$nli_mse_eft_param+sdf35_t1e$nli_mse_ncc_param)/2, col = "orange")
legend(x = "bottomleft", 
       legend = c("Young", "Lin. Hier.", "Nonlin."),
       col = c("red", "blue", "orange"), lwd = 2)

# Type 1 Error scenarios (N=70)
plot(sdf70_t1e$true_eft_mean, sdf70_t1e$nl_mse_eft_param, type = "l", col = "red", main = "EFT Param Mean Square Error", sub = "EFT = NCC")
lines(sdf70_t1e$true_eft_mean, sdf70_t1e$lin_mse_eft_param, col = "blue")
lines(sdf70_t1e$true_eft_mean, sdf70_t1e$nli_mse_eft_param, col = "orange")

plot(sdf70_t1e$true_eft_mean, sdf70_t1e$nl_mse_ncc_param, type = "l", col = "red", main = "NCC Param Mean Square Error", sub = "EFT = NCC")
lines(sdf70_t1e$true_eft_mean, sdf70_t1e$lin_mse_ncc_param, col = "blue")
lines(sdf70_t1e$true_eft_mean, sdf70_t1e$nli_mse_ncc_param, col = "orange")

plot(sdf70_t1e$true_eft_mean, (sdf70_t1e$nl_mse_eft_param+sdf70_t1e$nl_mse_ncc_param)/2, type = "l", col = "red", main = "Aggregate Param Mean Square Error", sub = "EFT = NCC")
lines(sdf70_t1e$true_eft_mean, (sdf70_t1e$lin_mse_eft_param+sdf70_t1e$lin_mse_ncc_param)/2, col = "blue")
lines(sdf70_t1e$true_eft_mean, (sdf70_t1e$nli_mse_eft_param+sdf70_t1e$nli_mse_ncc_param)/2, col = "orange")


# Type 1 Error scenarios (N=140)
plot(sdf140_t1e$true_eft_mean, sdf140_t1e$nl_mse_eft_param, type = "l", col = "red", main = "EFT Param Mean Square Error", sub = "EFT = NCC")
lines(sdf140_t1e$true_eft_mean, sdf140_t1e$lin_mse_eft_param, col = "blue")
lines(sdf140_t1e$true_eft_mean, sdf140_t1e$nli_mse_eft_param, col = "orange")

plot(sdf140_t1e$true_eft_mean, sdf140_t1e$nl_mse_ncc_param, type = "l", col = "red", main = "NCC Param Mean Square Error", sub = "EFT = NCC")
lines(sdf140_t1e$true_eft_mean, sdf140_t1e$lin_mse_ncc_param, col = "blue")
lines(sdf140_t1e$true_eft_mean, sdf140_t1e$nli_mse_ncc_param, col = "orange")

plot(sdf140_t1e$true_eft_mean, (sdf140_t1e$nl_mse_eft_param+sdf140_t1e$nl_mse_ncc_param)/2, type = "l", col = "red", main = "Aggregate Param Mean Square Error", sub = "EFT = NCC")
lines(sdf140_t1e$true_eft_mean, (sdf140_t1e$lin_mse_eft_param+sdf140_t1e$lin_mse_ncc_param)/2, col = "blue")
lines(sdf140_t1e$true_eft_mean, (sdf140_t1e$nli_mse_eft_param+sdf140_t1e$nli_mse_ncc_param)/2, col = "orange")

# Power scenarios (N=140)
plot(sdf140_pow_55$true_eft_mean, sdf140_pow_55$nl_mse_eft_param, type = "l", col = "red", main = "EFT Param Mean Square Error", sub = "NCC = -5.5")
lines(sdf140_pow_55$true_eft_mean, sdf140_pow_55$lin_mse_eft_param, col = "blue")
lines(sdf140_pow_55$true_eft_mean, sdf140_pow_55$nli_mse_eft_param, col = "orange")

plot(sdf140_pow$true_eft_mean, sdf140_pow$nl_mse_ncc_param, type = "l", col = "red", main = "NCC Param Mean Square Error", sub = "NCC = -5.5")
lines(sdf140_pow$true_eft_mean, sdf140_pow$lin_mse_ncc_param, col = "blue")

# Power scenarios (N=35)
plot(sdf35_pow$true_eft_mean, sdf35_pow$nl_mse_eft_param, type = "l", col = "red", main = "EFT Param Mean Square Error", ylim = c(0.07,0.1))
lines(sdf35_pow$true_eft_mean, sdf35_pow$lin_mse_eft_param, col = "blue")
lines(sdf35_pow$true_eft_mean, sdf35_pow$nli_mse_eft_param, col = "orange")


plot(sdf35_pow$true_eft_mean, sdf35_pow$nl_mse_ncc_param, type = "l", col = "red", main = "NCC Param Mean Square Error", ylim = c(0.07,0.1))
lines(sdf35_pow$true_eft_mean, sdf35_pow$lin_mse_ncc_param, col = "blue")



# Bias plot

# Type 1 Error scenarios
# N=35
plot(sdf35_t1e$true_eft_mean, sdf35_t1e$nl_mean_est_eft_param - sdf35_t1e$true_eft_mean, type = "l", col = "red", main = "EFT Param Bias",
     xlab="True EFT Mean", ylab = "Bias of EFT mean estimator", ylim = c(-0.06, 0.12))
lines(sdf35_t1e$true_eft_mean, sdf35_t1e$lin_mean_est_eft_param - sdf35_t1e$true_eft_mean, type = "l", col = "blue")
lines(sdf35_t1e$true_eft_mean, sdf35_t1e$nli_mean_est_eft_param - sdf35_t1e$true_eft_mean, type = "l", col = "orange")

# N=70
plot(sdf70_t1e$true_eft_mean, sdf70_t1e$nl_mean_est_eft_param - sdf70_t1e$true_eft_mean, type = "l", col = "red", main = "EFT Param Bias",
     xlab="True EFT Mean", ylab = "Bias of EFT mean estimator", ylim = c(-0.06, 0.12))
lines(sdf70_t1e$true_eft_mean, sdf70_t1e$lin_mean_est_eft_param - sdf70_t1e$true_eft_mean, type = "l", col = "blue")
lines(sdf70_t1e$true_eft_mean, sdf70_t1e$nli_mean_est_eft_param - sdf70_t1e$true_eft_mean, type = "l", col = "orange")


# N=140
plot(sdf140_t1e$true_eft_mean, sdf140_t1e$nl_mean_est_eft_param - sdf140_t1e$true_eft_mean, type = "l", col = "red", main = "EFT Param Bias",
     xlab="True EFT Mean", ylab = "Bias of EFT mean estimator", ylim = c(-0.06, 0.12), lwd=2)
lines(sdf140_t1e$true_eft_mean, sdf140_t1e$lin_mean_est_eft_param - sdf140_t1e$true_eft_mean, type = "l", col = "blue", lwd=2)
lines(sdf140_t1e$true_eft_mean, sdf140_t1e$nli_mean_est_eft_param - sdf140_t1e$true_eft_mean, type = "l", col = "orange", lwd=2)






# Power scenarios

plot(sdf140_pow$true_eft_mean, sdf140_pow$nl_mean_est_eft_param - sdf140_pow$true_eft_mean, type = "l", col = "red", main = "EFT Param Bias")
lines(sdf140_pow$true_eft_mean, sdf140_pow$lin_mean_est_eft_param - sdf140_pow$true_eft_mean, type = "l", col = "blue")



plot(sdf140_pow$true_eft_mean, sdf140_pow$nl_mean_est_ncc_param - sdf140_pow$true_ncc_mean, type = "l", col = "red", main = "NCC Param Bias", ylim = c(-0.06, 0.06))
lines(sdf140_pow$true_eft_mean, sdf140_pow$lin_mean_est_ncc_param - sdf140_pow$true_ncc_mean, type = "l", col = "blue")


plot(sdf35_pow$true_eft_mean, sdf35_pow$nl_mean_est_eft_param - sdf35_pow$true_eft_mean, type = "l", col = "red", main = "EFT Param Bias")
lines(sdf35_pow$true_eft_mean, sdf35_pow$lin_mean_est_eft_param - sdf35_pow$true_eft_mean, type = "l", col = "blue")

plot(sdf35_t1e$true_eft_mean, sdf35_t1e$nl_mean_est_eft_param - sdf35_t1e$true_eft_mean, type = "l", col = "red", main = "EFT Param Bias")
lines(sdf35_t1e$true_eft_mean, sdf35_t1e$lin_mean_est_eft_param - sdf35_t1e$true_eft_mean, type = "l", col = "blue")
lines(sdf35_t1e$true_eft_mean, sdf35_t1e$nli_mean_est_eft_param - sdf140_t1e$true_eft_mean, type = "l", col = "orange")


plot(sdf35_pow$true_eft_mean, sdf35_pow$nl_mean_est_ncc_param - sdf35_pow$true_ncc_mean, type = "l", col = "red", main = "NCC Param Bias", ylim = c(-0.06, 0.06))
lines(sdf35_pow$true_eft_mean, sdf35_pow$lin_mean_est_ncc_param - sdf35_pow$true_ncc_mean, type = "l", col = "blue")



# Subject MSEs

# Type 1 Error Scenarios
# N=35
plot(sdf35_t1e$true_eft_mean, sdf35_t1e$nl_mse_eft_subj, type = "l", col = "red", main = "EFT subj Mean Square Error", sub = "EFT = NCC", ylim = c(0, 0.7))
lines(sdf35_t1e$true_eft_mean, sdf35_t1e$lin_mse_eft_subj, col = "blue")
lines(sdf35_t1e$true_eft_mean, sdf35_t1e$nli_mse_eft_subj, col = "orange")
legend(x = "bottomleft", 
       legend = c("Young", "Lin. Hier.", "Nonlin."),
       col = c("red", "blue", "orange"), lwd = 2)


plot(sdf35_t1e$true_eft_mean, sdf35_t1e$nl_mse_ncc_subj, type = "l", col = "red", main = "NCC subj Mean Square Error", sub = "EFT = NCC", ylim = c(0,0.7))
lines(sdf35_t1e$true_ncc_mean, sdf35_t1e$lin_mse_ncc_subj, col = "blue")
lines(sdf35_t1e$true_ncc_mean, sdf35_t1e$nli_mse_ncc_subj, col = "orange")
legend(x = "bottomleft", 
       legend = c("Young", "Lin. Hier.", "Nonlin."),
       col = c("red", "blue", "orange"), lwd = 2)


# N=70
plot(sdf70_t1e$true_eft_mean, sdf70_t1e$nl_mse_eft_subj, type = "l", col = "red", main = "EFT subj Mean Square Error", sub = "EFT = NCC", ylim = c(0, 0.7))
lines(sdf70_t1e$true_eft_mean, sdf70_t1e$lin_mse_eft_subj, col = "blue")
lines(sdf70_t1e$true_eft_mean, sdf70_t1e$nli_mse_eft_subj, col = "orange")


plot(sdf70_t1e$true_eft_mean, sdf70_t1e$nl_mse_ncc_subj, type = "l", col = "red", main = "NCC subj Mean Square Error", sub = "EFT = NCC", ylim = c(0,0.7))
lines(sdf70_t1e$true_ncc_mean, sdf70_t1e$lin_mse_ncc_subj, col = "blue")
lines(sdf70_t1e$true_ncc_mean, sdf70_t1e$nli_mse_ncc_subj, col = "orange")

# N=140
plot(sdf140_t1e$true_eft_mean, sdf140_t1e$nl_mse_eft_subj, type = "l", col = "red", main = "EFT subj Mean Square Error", sub = "EFT = NCC", ylim = c(0, 0.7))
lines(sdf140_t1e$true_eft_mean, sdf140_t1e$lin_mse_eft_subj, col = "blue")
lines(sdf140_t1e$true_eft_mean, sdf140_t1e$nli_mse_eft_subj, col = "orange")


plot(sdf140_t1e$true_eft_mean, sdf140_t1e$nl_mse_ncc_subj, type = "l", col = "red", main = "NCC subj Mean Square Error", sub = "EFT = NCC", ylim = c(0,0.7))
lines(sdf140_t1e$true_ncc_mean, sdf140_t1e$lin_mse_ncc_subj, col = "blue")
lines(sdf140_t1e$true_ncc_mean, sdf140_t1e$nli_mse_ncc_subj, col = "orange")


# Power scenarios
plot(sdf140_pow$true_eft_mean, sdf140_pow$nl_mse_eft_subj, type = "l", col = "red", main = "EFT subj Mean Square Error", sub = "NCC = -5.5", ylim = c(0, 0.5))
lines(sdf140_pow$true_eft_mean, sdf140_pow$lin_mse_eft_subj, col = "blue")

plot(sdf140_pow$true_eft_mean, sdf140_pow$nl_mse_ncc_subj, type = "l", col = "red", main = "NCC subj Mean Square Error", sub = "NCC = -5.5")
lines(sdf140_pow$true_eft_mean, sdf140_pow$lin_mse_ncc_subj, col = "blue")

# Power scenarios (N=35)
plot(sdf35_pow$true_eft_mean, sdf35_pow$nl_mse_eft_subj, type = "l", col = "red", main = "EFT subj Mean Square Error", ylim = c(0,0.8))
lines(sdf35_pow$true_eft_mean, sdf35_pow$lin_mse_eft_subj, col = "blue")

plot(sdf35_pow$true_eft_mean, sdf35_pow$nl_mse_ncc_subj, type = "l", col = "red", main = "NCC subj Mean Square Error", ylim = c(0,0.8))
lines(sdf35_pow$true_eft_mean, sdf35_pow$lin_mse_ncc_subj, col = "blue")

# t1e scenarios
plot(sdf35_t1e$true_eft_mean, sdf35_t1e$nl_mse_eft_subj, type = "l", col = "red", main = "EFT subj Mean Square Error", ylim = c(0,0.8))
lines(sdf35_t1e$true_eft_mean, sdf35_t1e$lin_mse_eft_subj, col = "blue")

plot(sdf35_t1e$true_eft_mean, sdf35_t1e$nl_mse_ncc_subj, type = "l", col = "red", main = "NCC subj Mean Square Error", ylim = c(0,0.8))
lines(sdf35_t1e$true_eft_mean, sdf35_t1e$lin_mse_ncc_subj, col = "blue")



# Filtration proportion by rule check
plot(sdf140_t1e$true_eft_mean, (sdf140_t1e$prop_rc_lost_eft + sdf140_t1e$prop_rc_lost_ncc)/2, type="l", lwd=2,
     main = "Rule check failures by true ln(k)", xlab = "log(k)", ylab = "Rule check failure proportion")

# Get as a data frame, the filtration proportion
filt_prop <- sdf140_t1e %>%
  mutate(prop_fail_rc = (prop_rc_lost_eft+prop_rc_lost_ncc)/2,
         N = (n_eft+n_ncc)*sims) %>%
  rename(true_ln_k_mean = true_eft_mean) %>%
  ungroup() %>%
  select(true_ln_k_mean, sigma_sq, g, N, prop_fail_rc) %>%
  select(true_ln_k_mean, prop_fail_rc)  # Less complete but simpler version

# write.csv(filt_prop, "Plots/filt_prop.csv", row.names = FALSE)


t1e_n35 <- sdf35_t1e %>%  
  ungroup() %>%
  rename(true_ln_k_mean = true_eft_mean) %>%
  mutate(nl_param_mse = (nl_mse_eft_param+nl_mse_ncc_param)/2,
         lin_param_mse = (lin_mse_eft_param+lin_mse_ncc_param)/2,
         nli_param_mse = (nli_mse_eft_param+nli_mse_ncc_param)/2,
         nl_subj_mse = (nl_mse_eft_subj+nl_mse_ncc_subj)/2,
         lin_subj_mse = (lin_mse_eft_subj+lin_mse_ncc_subj)/2,
         nli_subj_mse = (nli_mse_eft_subj+nli_mse_ncc_subj)/2,
         nl_bias = (nl_mean_est_eft_param + nl_mean_est_ncc_param)/2 - true_ln_k_mean,
         lin_bias = (lin_mean_est_eft_param + lin_mean_est_ncc_param)/2 - true_ln_k_mean,
         nli_bias = (nli_mean_est_eft_param + nli_mean_est_ncc_param)/2 - true_ln_k_mean) %>%
  select(n_eft, n_ncc, true_ln_k_mean, sigma_sq, g,   # Parameters
         nl_power, lin_power, nli_power, lin_rc_power, nli_rc_power, # Type 1 Errors
         nl_param_mse, lin_param_mse, nli_param_mse,  # hyperparameter ln(k) MSE
         nl_subj_mse, lin_subj_mse, nli_subj_mse,  # subject ln(k) MSE
         nl_bias, lin_bias, nli_bias)  
write.csv(t1e_n35, "Plots/t1e_n35.csv", row.names = FALSE)

t1e_n70 <- sdf70_t1e %>%  
  ungroup() %>%
  rename(true_ln_k_mean = true_eft_mean) %>%
  mutate(nl_param_mse = (nl_mse_eft_param+nl_mse_ncc_param)/2,
         lin_param_mse = (lin_mse_eft_param+lin_mse_ncc_param)/2,
         nli_param_mse = (nli_mse_eft_param+nli_mse_ncc_param)/2,
         nl_subj_mse = (nl_mse_eft_subj+nl_mse_ncc_subj)/2,
         lin_subj_mse = (lin_mse_eft_subj+lin_mse_ncc_subj)/2,
         nli_subj_mse = (nli_mse_eft_subj+nli_mse_ncc_subj)/2,
         nl_bias = (nl_mean_est_eft_param + nl_mean_est_ncc_param)/2 - true_ln_k_mean,
         lin_bias = (lin_mean_est_eft_param + lin_mean_est_ncc_param)/2 - true_ln_k_mean,
         nli_bias = (nli_mean_est_eft_param + nli_mean_est_ncc_param)/2 - true_ln_k_mean) %>%
  select(n_eft, n_ncc, true_ln_k_mean, sigma_sq, g,   # Parameters
         nl_power, lin_power, nli_power, lin_rc_power, nli_rc_power, # Type 1 Errors
         nl_param_mse, lin_param_mse, nli_param_mse,  # hyperparameter ln(k) MSE
         nl_subj_mse, lin_subj_mse, nli_subj_mse,  # subject ln(k) MSE
         nl_bias, lin_bias, nli_bias)  
write.csv(t1e_n70, "Plots/t1e_n70.csv", row.names = FALSE)

t1e_n140 <- sdf140_t1e %>%  
  ungroup() %>%
  rename(true_ln_k_mean = true_eft_mean) %>%
  mutate(nl_param_mse = (nl_mse_eft_param+nl_mse_ncc_param)/2,
         lin_param_mse = (lin_mse_eft_param+lin_mse_ncc_param)/2,
         nli_param_mse = (nli_mse_eft_param+nli_mse_ncc_param)/2,
         nl_subj_mse = (nl_mse_eft_subj+nl_mse_ncc_subj)/2,
         lin_subj_mse = (lin_mse_eft_subj+lin_mse_ncc_subj)/2,
         nli_subj_mse = (nli_mse_eft_subj+nli_mse_ncc_subj)/2,
         nl_bias = (nl_mean_est_eft_param + nl_mean_est_ncc_param)/2 - true_ln_k_mean,
         lin_bias = (lin_mean_est_eft_param + lin_mean_est_ncc_param)/2 - true_ln_k_mean,
         nli_bias = (nli_mean_est_eft_param + nli_mean_est_ncc_param)/2 - true_ln_k_mean) %>%
  select(n_eft, n_ncc, true_ln_k_mean, sigma_sq, g,   # Parameters
         nl_power, lin_power, nli_power, lin_rc_power, nli_rc_power, # Type 1 Errors
         nl_param_mse, lin_param_mse, nli_param_mse,  # hyperparameter ln(k) MSE
         nl_subj_mse, lin_subj_mse, nli_subj_mse,  # subject ln(k) MSE
         nl_bias, lin_bias, nli_bias)  
write.csv(t1e_n140, "Plots/t1e_n140.csv", row.names = FALSE)
