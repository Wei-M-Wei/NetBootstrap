rm(list = ls())
library(NetBootstrap)
library(doParallel)
library (parallel)
library(foreach)
library(alpaca)
library(fixest)
library(speedglm)
source('simulation data.R')
packages_to_export <- c("dplyr", 'NetBootstrap', 'alpaca')
num_cores <- detectCores()
cl <- makeCluster(num_cores-7)
registerDoParallel(cl)

# probit model network
N = 30
N_seq = c(N)

# bootstrap times
bootstrap_time = 499

# repitition times
mle_num = 100

for(design in c(1,2,3,4)){

  # empty matrix
  Estimate_bias = matrix(0, length(N_seq), 6)
  Estimate_se = matrix(0, length(N_seq), 6)
  Estimate_bias_median = matrix(0, length(N_seq), 6)
  percentile = matrix(0, length(N_seq), 3)
  p_reject_ratio = matrix(0, length(N_seq), 3)
  p_rej_ratio_without = matrix(0, length(N_seq), 3)
  APE_estimate = matrix(0, length(N_seq), 3)
  Estimate_deviation_boot = matrix(0, length(N_seq), 1)
  Estimate_deviation = matrix(0, length(N_seq), 6)
  APE_deviation = matrix(0, length(N_seq), 3)
  Hessian_inv = matrix(0, length(N_seq), 1)
  cov_MLE = matrix(0, length(N_seq), 1)
  cov_jack = matrix(0, length(N_seq), 1)
  boost_variance = matrix(0, length(N_seq), 1)
  ratio_MLE = matrix(0, length(N_seq), 1)
  ratio_b = matrix(0, length(N_seq), 1)

  # parameters
  beta = c(1)
  beta_NULL = c(1)
  t_index = 1
  K = length(beta)
  design = 1
  for (N in N_seq) {
    result = foreach(k = 1:mle_num, .errorhandling = 'remove', .packages = packages_to_export) %dopar% {
      set.seed(k)
      cof = NULL
      cof_j = NULL
      cof_constrain = NULL
      cof_constrain_j = NULL
      estimate_jack = NULL

      # prepare the data
      DGP = data_generation(N = N, beta = beta, design = design )
      y = DGP$y
      X_design = DGP$X
      data_in = as.data.frame(DGP$data)
      index_name = colnames(data_in)[(K+2):(K+3)]

      # estimation
      fit = network_bootstrap(y = y, X = X_design, N = N, bootstrap_time = bootstrap_time, index = index_name, data = data_in, link = 'probit', beta_NULL = beta_NULL)
      fit_jack = split_jackknife(y = y, X = X_design, N = N, index = index_name, data = data_in, link = 'probit', beta_NULL = NULL)
      fit_analytical = analytical_Amrei(y = y, X = X_design, N = N, index = index_name, data = data_in, link = 'probit', L = 1, beta_NULL = NULL)
      fit_analytical_own = analytical_corrected(y = y, X = X_design, N = N, index = index_name, data = data_in, link = 'probit', L = 1, beta_NULL = NULL)

      # all the estimates
      MLE_estimate = fit$cof_MLE[1]
      Medain_bootstrap_estimate = fit$cof_median[1]
      Mean_bootstrap_estimate = fit$cof_mean[1]
      jack_estimate = fit_jack$cof_jack[1]
      analytical_estimate = fit_analytical$cm[1]
      analytical_estimate_own = fit_analytical_own$est
      bootstrap_estimate = fit$cof_bootstrap_all[,1]
      mean_bias = mean(fit$cof_bootstrap_all[,1]) - fit$cof_MLE[1]
      se_MLE = sqrt(fit$Hessian_MLE[1,1])
      se_bootstrap = fit$sd
      se_jack = fit_jack$se_weidner
      se_analytical = fit_analytical$cm[2]
      se_analytical_own = fit_analytical_own$se


      # several evaluations
      if (MLE_estimate  - quantile(fit$cof_bootstrap_all[,1] - MLE_estimate, 0.975) <= beta & beta <= MLE_estimate  - quantile(fit$cof_bootstrap_all[,1] - MLE_estimate, 0.025)){
        p_cover_bootstrap  =  1
      }else{
        p_cover_bootstrap  =  0
      }

      if (MLE_estimate  - 1.96*se_MLE <= beta & beta <= MLE_estimate  + 1.96*se_MLE ){
        p_cover_MLE  =  1
      }else{
        p_cover_MLE  =  0
      }

      if (jack_estimate  - 1.96*se_jack <= beta & beta <= jack_estimate  + 1.96*se_jack ){
        p_cover_jack  =  1
      }else{
        p_cover_jack  =  0
      }

      if (analytical_estimate  - 1.96*se_analytical <= beta & beta <= analytical_estimate  + 1.96*se_analytical ){
        p_cover_analytical  =  1
      }else{
        p_cover_analytical  =  0
      }

      if (analytical_estimate_own  - 1.96*se_analytical_own <= beta & beta <= analytical_estimate_own  + 1.96*se_analytical_own ){
        p_cover_analytical_own  =  1
      }else{
        p_cover_analytical_own  =  0
      }



      results = list( MLE_estimate = MLE_estimate, Mean_bootstrap_estimate  = Mean_bootstrap_estimate , Medain_bootstrap_estimate = Medain_bootstrap_estimate, bootstrap_estimate = bootstrap_estimate,
                      jack_estimate = jack_estimate, analytical_estimate = analytical_estimate, analytical_estimate_own = analytical_estimate_own,
                      bootstrap_estimate_low = quantile(bootstrap_estimate - MLE_estimate, 0.025), bootstrap_estimate_up = quantile(bootstrap_estimate - MLE_estimate, 0.975),
                      cov_var_MLE = fit$Hessian_MLE[1,1], cov_var_jack = fit$Hessian_MLE[1,1],
                      se_MLE = se_MLE, se_bootstrap = se_bootstrap, se_jack = se_jack, se_analytical = se_analytical, se_analytical_own = se_analytical_own,
                      p_cover_MLE = p_cover_MLE, p_cover_bootstrap = p_cover_bootstrap, p_cover_jack = p_cover_jack,
                      p_cover_analytical = p_cover_analytical,p_cover_analytical_own = p_cover_analytical_own)
      return(results)
    }



    # svae the results
    result <- Filter(Negate(is.null), result)

    if (length(result) == 0){
      print('run again')
      t_index = t_index + 1
    } else {
      all_estimator = matrix(0, length(result), 6)
      Boot_estimate = matrix(0, length(result), 2)
      se_estimator = matrix(0, length(result), 6)
      p_cover_bootstrap = 0
      p_cover_MLE = 0
      p_cover_jack = 0
      p_cover_analytical = 0
      p_cover_analytical_own = 0

      for (i in seq(length(result))) {
        all_estimator[i, ] = c(result[[i]]$MLE_estimate, result[[i]]$Mean_bootstrap_estimate, result[[i]]$Medain_bootstrap_estimate, result[[i]]$jack_estimate, result[[i]]$analytical_estimate, result[[i]]$analytical_estimate_own)
        se_estimator[i, ] = c(result[[i]]$se_MLE, result[[i]]$se_bootstrap,result[[i]]$se_bootstrap, result[[i]]$se_jack, result[[i]]$se_analytical, result[[i]]$se_analytical_own)
        p_cover_MLE = p_cover_MLE + result[[i]]$p_cover_MLE
        p_cover_jack = p_cover_jack + result[[i]]$p_cover_jack
        p_cover_analytical = p_cover_analytical + result[[i]]$p_cover_analytical
        p_cover_analytical_own = p_cover_analytical_own + result[[i]]$p_cover_analytical_own
        p_cover_bootstrap = p_cover_bootstrap + result[[i]]$p_cover_bootstrap

      }

      # beta bias
      Estimate_bias[t_index, ] = apply(all_estimator - beta, 2, mean)
      Estimate_se[t_index, ] = apply(se_estimator, 2, mean)
      Estimate_deviation[t_index, ] = apply(all_estimator , 2, sd)

      # coverage rate for beta
      cover_MLE = p_cover_MLE/length(result)
      cover_jack = p_cover_jack/length(result)
      cover_analytical = p_cover_analytical/length(result)
      cover_analytical_own = p_cover_analytical_own/length(result)
      cover_bootstrap = p_cover_bootstrap/length(result)

      }
    t_index = t_index + 1
  }

  # save the table
  table_name1 <- paste0('setting_', design, '_bias_estimation', ".csv")
  table_name2 <- paste0('setting_', design, '_cover_rate', ".csv")
  table_name3 <- paste0('setting_', design, '_se', ".csv")
  table_name5 <- paste0('setting_', design, 'all_estimator', ".csv")
  write.csv(rbind( Estimate_bias, Estimate_deviation ) , table_name1)
  write.csv(rbind(cover_MLE, cover_bootstrap, cover_jack, cover_analytical,cover_analytical_own) , table_name2)
  write.csv(Estimate_se , table_name3)
  write.csv(all_estimator , table_name5)

}

