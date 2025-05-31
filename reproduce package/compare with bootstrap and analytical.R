rm(list = ls())
library(NetBootstrap)
library(doParallel)
library (parallel)
library(foreach)
library(alpaca)
library(fixest)
library(speedglm)
library(margins)
library(MASS)

# functions needed

source('simulation data.R')
packages_to_export <- c("dplyr", 'NetBootstrap', 'alpaca', 'speedglm', 'fixest', 'margins')
num_cores <- detectCores()
cl <- makeCluster(num_cores-25)
registerDoParallel(cl)

# probit model network
N = 40
N_seq = c(N)

# bootstrap times
bootstrap_time = 19

# repitition times
mle_num = 300
design = 3

# parameters
beta = c(1)
beta_NULL = c(1, 2)
link = 'probit'
K = length(beta)

for(design in c(3)){
  start_time <- Sys.time()
  t_index = 1
  # empty matrix
  Estimate_bias = matrix(0, length(N_seq), 5*K)
  Estimate_se = matrix(0, length(N_seq), 5*K)
  Estimate_bias_median = matrix(0, length(N_seq), 5*K)
  percentile = matrix(0, length(N_seq), 3*K)
  p_reject_ratio = matrix(0, length(N_seq), 3*K)
  p_rej_ratio_without = matrix(0, length(N_seq), 3*K)
  APE_estimate = matrix(0, length(N_seq), 5*K)
  Estimate_deviation_boot = matrix(0, length(N_seq), 1*K)
  Estimate_deviation = matrix(0, length(N_seq), 5*K)
  APE_deviation = matrix(0, length(N_seq), 5*K)
  Hessian_inv = matrix(0, length(N_seq), 1*K)
  cov_MLE = matrix(0, length(N_seq), 1*K)
  cov_jack = matrix(0, length(N_seq), 1*K)
  boost_variance = matrix(0, length(N_seq), 1*K)
  ratio_MLE = matrix(0, length(N_seq), 1*K)
  ratio_b = matrix(0, length(N_seq), 1*K)
  for (N in N_seq) {
    result = foreach(k = 1:mle_num, .errorhandling = 'remove', .packages = packages_to_export) %dopar% {
      set.seed(k)
      cof = NULL
      cof_j = NULL
      cof_constrain = NULL
      cof_constrain_j = NULL
      estimate_jack = NULL

      # prepare the data
      DGP = data_generation( N = N, beta = beta, design = design, model = link)
      y = DGP$y
      X_design = DGP$X
      data_in = as.data.frame(DGP$data)
      index_name = colnames(data_in)[(K+2):(K+3)]

      # estimation
      fit = network_bootstrap(y = y, X = X_design, N = N, bootstrap_time = bootstrap_time, index = index_name, data = data_in, link = link, boot_repeat = 1)
      fit_analytical = analytical_Amrei(y = y, X = X_design, N = N, index = index_name, data = data_in, link = link, L = 1, beta_NULL = NULL)
      fit_analytical_own = analytical_corrected(y = y, X = X_design, N = N, index = index_name, data = data_in, link = link, L = 1, beta_NULL = NULL)

      # vector_to_matrix(fit_analytical_own$eta, N, data_in$index.1, data_in$index.2)[1,]
      # fit_analytical$est$eta
      # all the estimates
      MLE_estimate = fit$est_MLE[1:K]
      Medain_bootstrap_estimate = fit$est_median[1:K]
      Mean_bootstrap_estimate = fit$est_mean[1:K]
      analytical_estimate = summary(fit_analytical$est)$cm[1:K,1]
      analytical_estimate_own = fit_analytical_own$est[1:K]
      bootstrap_estimate = fit$est_bootstrap_all[,1:K]
      mean_bias = apply(as.matrix(fit$est_bootstrap_all[,1:K]), 2 , mean) - fit$est_MLE[1:K]
      se_MLE = fit_analytical_own$se_MLE
      se_bootstrap = fit$sd
      se_bootstrap_critical = fit$sd_critical
      se_analytical = summary(fit_analytical$est)$cm[1:K,2]
      se_analytical_own = fit_analytical_own$se

      low_alpha = find_critical_low(fit$est_bootstrap_all_critical[,1], fit$est_bootstrap_all[,1] - fit$est_critical[1:K], MLE_estimate, 0.025 )
      up_alpha = find_critical_up(fit$est_bootstrap_all_critical[,1], fit$est_bootstrap_all[,1] - fit$est_critical[1:K], MLE_estimate, 0.025 )
      p_cover_bootstrap <- as.integer(
        (MLE_estimate - quantile(fit$est_bootstrap_all_critical[,1] - MLE_estimate, up_alpha) <= beta) &
          (beta <= MLE_estimate - quantile(fit$est_bootstrap_all_critical[,1] - MLE_estimate, low_alpha))
      )

      p_cover_bootstrap_basic <- as.integer(
        (MLE_estimate - quantile(fit$est_bootstrap_all_critical[,1] - MLE_estimate, 0.975) <= beta) &
          (beta <= MLE_estimate - quantile(fit$est_bootstrap_all_critical[,1] - MLE_estimate, 0.025))
      )

      p_cover_MLE <- as.integer(
        (MLE_estimate - 1.96*se_MLE <= beta) & (beta <= MLE_estimate + 1.96*se_MLE)
      )


      p_cover_analytical <- as.integer(
        (analytical_estimate - 1.96*se_analytical <= beta) & (beta <= analytical_estimate + 1.96*se_analytical)
      )

      p_cover_analytical_own <- as.integer(
        (analytical_estimate_own - 1.96*se_analytical_own <= beta) & (beta <= analytical_estimate_own + 1.96*se_analytical_own)
      )


      results = list( MLE_estimate = MLE_estimate, Mean_bootstrap_estimate  = fit$est_critical[1:K] , Medain_bootstrap_estimate = Medain_bootstrap_estimate, bootstrap_estimate = bootstrap_estimate,
                      analytical_estimate = analytical_estimate, analytical_estimate_own = analytical_estimate_own,
                      bootstrap_estimate_low = quantile(bootstrap_estimate - MLE_estimate, 0.025), bootstrap_estimate_up = quantile(bootstrap_estimate - MLE_estimate, 0.975),
                      cov_var_MLE = fit$Hessian_MLE[1,1], cov_var_jack = fit$Hessian_MLE[1,1],
                      se_MLE = se_MLE, se_bootstrap = se_bootstrap, se_bootstrap_critical = se_bootstrap_critical, se_analytical = se_analytical, se_analytical_own = se_analytical_own,
                      p_cover_MLE = p_cover_MLE, p_cover_bootstrap = p_cover_bootstrap, p_cover_bootstrap_basic = p_cover_bootstrap_basic,
                      p_cover_analytical = p_cover_analytical, p_cover_analytical_own = p_cover_analytical_own
      )
      return(results)
    }



    # save the results
    result <- Filter(Negate(is.null), result)

    if (length(result) == 0){
      print('run again')
      t_index = t_index + 1
    } else {
      all_estimator = matrix(0, length(result), 5*K)
      Boot_estimate_all = c()
      se_estimator = matrix(0, length(result), 5*K)
      p_cover_bootstrap = rep(0, K)
      p_cover_bootstrap_basic = rep(0, K)
      p_cover_MLE = rep(0, K)
      p_cover_jack = rep(0, K)
      p_cover_analytical = rep(0, K)
      p_cover_analytical_own = rep(0, K)
      p_cover_bootstrap_pre_pivoting = rep(0, K)

      for (i in seq(length(result))) {

        all_estimator[i, ] = c(result[[i]]$MLE_estimate, result[[i]]$Mean_bootstrap_estimate, result[[i]]$Medain_bootstrap_estimate, result[[i]]$analytical_estimate, result[[i]]$analytical_estimate_own)
        se_estimator[i, ] = c(result[[i]]$se_MLE, result[[i]]$se_bootstrap_critical, result[[i]]$se_bootstrap, result[[i]]$se_analytical, result[[i]]$se_analytical_own)
        p_cover_MLE = p_cover_MLE + result[[i]]$p_cover_MLE
        p_cover_analytical = p_cover_analytical + result[[i]]$p_cover_analytical
        p_cover_analytical_own = p_cover_analytical_own + result[[i]]$p_cover_analytical_own
        p_cover_bootstrap = p_cover_bootstrap + result[[i]]$p_cover_bootstrap
        p_cover_bootstrap_basic = p_cover_bootstrap_basic + result[[i]]$p_cover_bootstrap_basic
        # p_cover_bootstrap_pre_pivoting = p_cover_bootstrap_pre_pivoting + result[[i]]$p_cover_bootstrap_pre_pivoting
        # all bootstrap estimate
        Boot_estimate_all = rbind(Boot_estimate_all, result[[i]]$bootstrap_estimate)
      }

      # beta bias
      Estimate_bias[t_index, ] = apply(all_estimator, 2, mean) - rep(beta, 5)
      Estimate_se[t_index, ] = apply(se_estimator, 2, mean)
      Estimate_deviation[t_index, ] = apply(all_estimator , 2, sd)

      # coverage rate for beta
      cover_MLE = p_cover_MLE/length(result)
      cover_analytical = p_cover_analytical/length(result)
      cover_analytical_own = p_cover_analytical_own/length(result)
      cover_bootstrap = p_cover_bootstrap/length(result)
    }
    t_index = t_index + 1
  }
  end_time <- Sys.time()
  running_time = end_time - start_time

  # save the table
  table_name1 <- paste0('setting_', design, '_bias_estimation_time', running_time ,  ".csv")
  table_name2 <- paste0('setting_', design, '_cover_rate', running_time , ".csv")
  table_name3 <- paste0('setting_', design, 'se', running_time , ".csv")
  table_name4 <- paste0('setting_', design, 'all_estimator', running_time , ".csv")
  table_name5 <- paste0('setting_', design, '_cover_rate_APE', running_time , ".csv")
  table_name6 <- paste0('setting_', design, 'all_APE', running_time , ".csv")
  table_name7 <- paste0('setting_', design, 'all_APE_se', running_time , ".csv")
  table_name8 <- paste0('setting_', design, 'APE_bias_estimation', running_time , ".csv")
  table_name9 <- paste0('setting_', design, 'all_APE_se_own', running_time, ".csv")
  table_name10 <- paste0('setting_', design, 'all_se', running_time, ".csv")
  table_name11 <- paste0('setting_', design, 'boot_estimate', running_time, ".csv")
  write.csv(rbind( Estimate_bias, Estimate_deviation ) , table_name1)
  write.csv(rbind(cover_MLE, cover_bootstrap_basic, cover_bootstrap, cover_analytical, cover_analytical_own) , table_name2)
  write.csv(Estimate_se , table_name3)
  write.csv(all_estimator , table_name4)
  # write.csv(cover_APE, table_name5)
  # write.csv(all_APE , table_name6)
  # write.csv(APE_se_formula, table_name7)
  # write.csv(rbind( APE_estimate, APE_deviation) , table_name8)
  # write.csv(APE_se_formula_own, table_name9)
  # write.csv(se_estimator, table_name10)
  write.csv(Boot_estimate_all, table_name11)

}

