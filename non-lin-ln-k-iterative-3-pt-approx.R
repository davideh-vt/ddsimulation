iter_ln_k_3_pt_approx_w_init <- function(dd_data, init_ln_k=-6, default_delta = 1, tol = 1e-4, max.iter = 20){
  # The input dataset dd_data needs to be sorted by group and then subject so that group_by does not change the order
  #  are together, in order of the time points
  
  # init_ln_k contains initial values for each 
  
  delta = default_delta
  
  # Get the time points from the dataset
  # TODO: Make absolutely sure this is sorted numerically smallest to largest
  time_points = as.numeric(levels(as.factor(dd_data$delay)))
  n_tp = length(time_points)
  
  # set up bound vectors
  ln_k_lb = rep(-Inf, length(init_ln_k))
  ln_k_ub = rep(Inf, length(init_ln_k))
  value_change = Inf
  middle_best = F
  
  consec_no_middle_best = -1
  
  bound_met = T
  
  iters = 0
  
  ### Loop will start here
  while(max(ln_k_ub - ln_k_lb) > tol & iters < max.iter){
  # TODO: Address max error not improving
  consec_no_middle_best = ifelse(middle_best, 0, 1+consec_no_middle_best)
    
  # check to make sure there is no conflict between the initial k and the bounds
  init_ln_k_to_bound = pmin(init_ln_k - ln_k_lb, ln_k_ub - init_ln_k)
  delta = ifelse(init_ln_k_to_bound/2 > default_delta, default_delta, init_ln_k_to_bound/2)
  delta = ifelse(bound_met & middle_best, pmin(value_change, delta), delta)
  
  # reset if too many not middle best
  do_reset = consec_no_middle_best > 3 & delta < default_delta & ln_k_lb > -Inf & ln_k_ub < Inf
  init_ln_k = ifelse(do_reset, (ln_k_lb + ln_k_ub)/2, init_ln_k)
  delta = ifelse(do_reset, (ln_k_ub - ln_k_lb)/4, delta)
  
  init_k = exp(init_ln_k) # Get the initial value of k. TODO: Is this unnecessary now?
  
  # set up matrices with modelled indifference points for 3 initial k per subject
  model_indiff_vals = matrix(data = 0, nrow = n_tp*length(init_ln_k), ncol = 3)
  model_indiff_vals[,2] = rep(exp(init_ln_k), each = n_tp)
  model_indiff_vals[,1] = rep(exp(init_ln_k - delta), each = n_tp)
  model_indiff_vals[,3] = rep(exp(init_ln_k + delta), each = n_tp)
  model_indiff_vals = model_indiff_vals*time_points  # at this point we have k*t
  model_indiff_vals = 1/(1+model_indiff_vals)
  
  # Calculate sse for each observation for each k
  dd_data$sse1 = (dd_data$indiff - model_indiff_vals[,1])^2
  dd_data$sse2 = (dd_data$indiff - model_indiff_vals[,2])^2
  dd_data$sse3 = (dd_data$indiff - model_indiff_vals[,3])^2
  
  # Get the sum square error for each subject
  subj_sse <- dd_data %>%
    group_by(study, group, subj) %>%
    summarise(sse1 = sum(sse1), sse2 = sum(sse2), sse3 = sum(sse3))
  
  subj_sse_mat = as.matrix(subj_sse[4:6])
  subj_sse_diff_mat = t(subj_sse_mat %*% matrix(c(1, -1, 0,
                                                  0, -1, 1), nrow = 3, ncol = 2))
  # we have in this matrix, the first column is sse1 - sse2; the second column is sse3 - sse2
  
  # establish bounds for the true solution based on the calculated sse
  
  left_pt_worse = subj_sse_diff_mat[1,] > 0  # figure out if sse1 is worse than sse2
  right_pt_worse = subj_sse_diff_mat[2,] > 0  # figure out if sse3 is worse than sse2
  ln_k_lb_1 = ifelse(left_pt_worse,
                     ifelse(right_pt_worse, init_ln_k - delta, init_ln_k),
                     ln_k_lb)
  ln_k_lb = ifelse(ln_k_lb_1 > ln_k_lb, ln_k_lb_1, ln_k_lb)
  ln_k_ub_1 = ifelse(right_pt_worse,
                     ifelse(left_pt_worse, init_ln_k + delta, init_ln_k),
                     ln_k_ub)
  ln_k_ub = ifelse(ln_k_ub_1 < ln_k_ub, ln_k_ub_1, ln_k_ub)
  middle_best = left_pt_worse & right_pt_worse
  
  
  
  # This can be simplified to only do 2 matrix multiplication operations, applied to t(subj_sse_mat)
  solver_mat = matrix(c(0.5, -0.5, 0.5, 0.5), ncol = 2, nrow = 2)
  coeff_mat = solver_mat %*% subj_sse_diff_mat
  ln_k_offsets = -coeff_mat[2,]/coeff_mat[1,]/2
  ln_k_approx = ln_k_offsets*delta + init_ln_k
  
  value_change = abs(init_ln_k - ln_k_approx)
  
  # Do not accept the approximation if it is not within the bounds
  bound_met = ln_k_approx > ln_k_lb & ln_k_approx < ln_k_ub
  bound_met = !is.na(bound_met) & bound_met
  ln_k_approx = ifelse(bound_met, ln_k_approx, 
                     ifelse(ln_k_lb == -Inf | ln_k_ub == Inf, 
                            ifelse(ln_k_lb == -Inf, ln_k_ub - 2*default_delta, ln_k_lb + 2*default_delta),
                            (ln_k_lb + ln_k_ub)/2))
  
  init_ln_k = ln_k_approx
  
  iters = iters + 1
  }
  ### Loop will end here
  
  
  # We might give up on solving some observations through the iterative approximation
  # In that case, pass them to the optimize method
  for(i in which(ln_k_ub - ln_k_lb > tol)){
    print(i)
    one_subj_data = dd_data %>%
      filter(study == subj_sse$study[i], group == subj_sse$group[i], subj == subj_sse$subj[i])
    ln_k_approx[i] = optimize.ln.k.subj(one_subj_data, c(max(ln_k_lb, -17), min(ln_k_ub, 2)))
  }
  
  # Once we have our approximation, add it to a dataframe identifying each estimate by subject and return
  subj_sse$nonlin_ln_k = ln_k_approx
  subj_est_ln_k_vals <- subj_sse %>%
    dplyr::select(study, group, subj, nonlin_ln_k)
  return(subj_est_ln_k_vals)
}


optimize.ln.k.subj <- function(one_subj_data, interval = c(-17,2)){
  # find optimal ln k for a single subject
  ln.k.hat <- optimize(f = optim.sse, 
                       interval = interval, 
                       time = one_subj_data$delay, 
                       discount = one_subj_data$indiff)$minimum
  return(ln.k.hat)
}

optim.sse <- function(ln_k, time, discount){
  # input:
  #  k: the discounting parameter
  #  time: time points
  #  discount: observed discount rate at that time
  # output: estimated ln_k
  
  k = exp(ln_k)
  pred_indiff = 1/(1+k*time)
  sse = sum((pred_indiff - discount)^2)
  return(sse)
}
