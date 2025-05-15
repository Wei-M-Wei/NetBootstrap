rm(list = ls())
library(NetBootstrap)
library(doParallel)
library (parallel)
library(foreach)
library(alpaca)
library(fixest)
library(speedglm)
library(margins)
source('simulation data.R')
packages_to_export <- c("dplyr", 'NetBootstrap', 'alpaca', 'speedglm', 'fixest', 'margins')
num_cores <- detectCores()
cl <- makeCluster(num_cores-5)
registerDoParallel(cl)

# probit model network
N = 20
N_seq = c(N)

# bootstrap times
bootstrap_time = 5

# repitition times
mle_num = 3

for(design in c(1,2,3,4)){
  start_time <- Sys.time()
  # empty matrix
  Estimate_bias = matrix(0, length(N_seq), 6)
  Estimate_se = matrix(0, length(N_seq), 6)
  Estimate_bias_median = matrix(0, length(N_seq), 6)
  percentile = matrix(0, length(N_seq), 3)
  p_reject_ratio = matrix(0, length(N_seq), 3)
  p_rej_ratio_without = matrix(0, length(N_seq), 3)
  APE_estimate = matrix(0, length(N_seq), 5)
  Estimate_deviation_boot = matrix(0, length(N_seq), 1)
  Estimate_deviation = matrix(0, length(N_seq), 6)
  APE_deviation = matrix(0, length(N_seq), 5)
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
  for (N in N_seq) {
    result = foreach(k = 1:mle_num, .errorhandling = 'remove', .packages = packages_to_export) %dopar% {
      set.seed(k)
      cof = NULL
      cof_j = NULL
      cof_constrain = NULL
      cof_constrain_j = NULL
      estimate_jack = NULL

      # prepare the data
      DGP = data_generation( N = N, beta = beta, design = design )
      y = DGP$y
      X_design = DGP$X
      data_in = as.data.frame(DGP$data)
      index_name = colnames(data_in)[(K+2):(K+3)]

      # estimation
      fit = network_bootstrap(y = y, X = X_design, N = N, bootstrap_time = bootstrap_time, index = index_name, data = data_in, link = 'probit', beta_NULL = NULL)
      fit_jack = split_jackknife(y = y, X = X_design, N = N, index = index_name, data = data_in, link = 'probit',  beta_NULL = NULL)
      fit_analytical = analytical_Amrei(y = y, X = X_design, N = N, index = index_name, data = data_in, link = 'probit', L = 1, beta_NULL = NULL)
      fit_analytical_own = analytical_corrected(y = y, X = X_design, N = N, index = index_name, data = data_in, link = 'probit', L = 1, beta_NULL = NULL)

      # all the estimates
      MLE_estimate = fit$est_MLE[1]
      Medain_bootstrap_estimate = fit$est_median[1]
      Mean_bootstrap_estimate = fit$est_mean[1]
      jack_estimate = fit_jack$est[1]
      analytical_estimate = summary(fit_analytical$est)$cm[1]
      analytical_estimate_own = fit_analytical_own$est[1]
      bootstrap_estimate = fit$est_bootstrap_all[,1]
      mean_bias = mean(fit$est_bootstrap_all[,1]) - fit$est_MLE[1]
      se_MLE = sqrt(fit$Hessian_MLE[1,1])
      se_bootstrap = fit$sd
      se_jack = fit_jack$se
      se_analytical = summary(fit_analytical$est)$cm[2]
      se_analytical_own = fit_analytical_own$se

      pre_pivoting_bootstrap = (Medain_bootstrap_estimate - beta)/se_analytical_own

      if ( pre_pivoting_bootstrap <= qnorm(0.975) & pre_pivoting_bootstrap >= qnorm(0.025)){
        p_cover_bootstrap_pre_pivoting  =  1
      }else{
        p_cover_bootstrap_pre_pivoting  =  0
      }


      # several evaluations
      if (MLE_estimate  - quantile(fit$est_bootstrap_all[,1] - MLE_estimate, 0.975) <= beta & beta <= MLE_estimate  - quantile(fit$est_bootstrap_all[,1] - MLE_estimate, 0.025)){
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

      APE_true = matrix(0, nrow = N, ncol = N)
      alpha = DGP$alpha
      gamma = DGP$gamma
      for (i in seq(N)){
        for ( j in seq(N)){
          if (i!=j)
            APE_true[i, j] = (pnorm(beta + alpha[i] + gamma[j]) - pnorm(-beta + alpha[i] + gamma[j]))/2
        }
      }

      # APE_true = matrix(0, nrow = N, ncol = N)
      # alpha = DGP$alpha
      # gamma = DGP$gamma
      # for (i in seq(N)){
      #   for ( j in seq(N)){
      #     if (i!=j)
      #       APE_true[i, j] = beta * dnorm(beta + alpha[i] + gamma[j])
      #   }
      # }

      APE_true = sum(APE_true)/(N*(N-1))
      APE_bootstrap = get_APE_bootstrap(y = y, X = X_design, N = N, data = data_in, fit = fit)
      APE_jackknife = get_APE_jackknife(y = y, X = X_design, N = N, data = data_in, index = index_name, fit = fit_jack)
      APE_analytical = get_APE_analytical(y = y, X = X_design, N = N, data = data_in, index = index_name, fit = fit_analytical_own, L = 1)
      APE_MLE = getAPEs(fit_analytical$est_MLE, panel.structure = 'network')
      APE_MLE_estimate = summary(APE_MLE)$cm[1]
      APE_MLE_se = summary(APE_MLE)$cm[2]
      APE_corrected = getAPEs(fit_analytical$est, panel.structure = 'network')
      APE_corrected_estimate = summary(APE_corrected)$cm[1]
      APE_corrected_se = summary(APE_corrected)$cm[2]
      APE_se_formula = c(APE_MLE_se, APE_bootstrap$se, APE_bootstrap$se, APE_corrected_se, APE_corrected_se)
      APE_se_formula_own = c(APE_MLE_se, APE_bootstrap$se, APE_bootstrap$se, APE_analytical$se, APE_analytical$se)

      results = list( MLE_estimate = MLE_estimate, Mean_bootstrap_estimate  = Mean_bootstrap_estimate , Medain_bootstrap_estimate = Medain_bootstrap_estimate, bootstrap_estimate = bootstrap_estimate,
                      jack_estimate = jack_estimate, analytical_estimate = analytical_estimate, analytical_estimate_own = analytical_estimate_own,
                      bootstrap_estimate_low = quantile(bootstrap_estimate - MLE_estimate, 0.025), bootstrap_estimate_up = quantile(bootstrap_estimate - MLE_estimate, 0.975),
                      cov_var_MLE = fit$Hessian_MLE[1,1], cov_var_jack = fit$Hessian_MLE[1,1],
                      se_MLE = se_MLE, se_bootstrap = se_bootstrap, se_jack = se_jack, se_analytical = se_analytical, se_analytical_own = se_analytical_own,
                      p_cover_MLE = p_cover_MLE, p_cover_bootstrap = p_cover_bootstrap, p_cover_jack = p_cover_jack,
                      p_cover_analytical = p_cover_analytical,p_cover_analytical_own = p_cover_analytical_own, p_cover_bootstrap_pre_pivoting = p_cover_bootstrap_pre_pivoting,
                      APE_true = APE_true, APE_MLE = fit$APE_MLE_estimate, APE_mean = APE_bootstrap$APE_mean_bootstrap, APE_median = APE_bootstrap$APE_median_bootstrap, APE_jackknife = APE_jackknife$APE,  APE_analytical =  APE_analytical$APE,  APE_corrected_estimate = APE_corrected_estimate,
                      APE_se_formula = APE_se_formula, APE_se_formula_own = APE_se_formula_own)
      return(results)
    }



    # svae the results
    result <- Filter(Negate(is.null), result)

    if (length(result) == 0){
      print('run again')
      t_index = t_index + 1
    } else {
      all_estimator = matrix(0, length(result), 6)
      Boot_estimate_all = matrix(0, length(result), bootstrap_time)
      se_estimator = matrix(0, length(result), 6)
      p_cover_bootstrap = 0
      p_cover_MLE = 0
      p_cover_jack = 0
      p_cover_analytical = 0
      p_cover_analytical_own = 0
      p_cover_bootstrap_pre_pivoting = 0
      p_pivoting = 0
      all_APE = matrix(0, length(result), 5)
      APE_se_formula = matrix(0, length(result), 5)
      APE_se_formula_own = matrix(0, length(result), 5)

      for (i in seq(length(result))) {
        all_estimator[i, ] = c(result[[i]]$MLE_estimate, result[[i]]$Mean_bootstrap_estimate, result[[i]]$Medain_bootstrap_estimate, result[[i]]$jack_estimate, result[[i]]$analytical_estimate, result[[i]]$analytical_estimate_own)
        se_estimator[i, ] = c(result[[i]]$se_MLE, result[[i]]$se_bootstrap,result[[i]]$se_bootstrap, result[[i]]$se_jack, result[[i]]$se_analytical, result[[i]]$se_analytical_own)
        p_cover_MLE = p_cover_MLE + result[[i]]$p_cover_MLE
        p_cover_jack = p_cover_jack + result[[i]]$p_cover_jack
        p_cover_analytical = p_cover_analytical + result[[i]]$p_cover_analytical
        p_cover_analytical_own = p_cover_analytical_own + result[[i]]$p_cover_analytical_own
        p_cover_bootstrap = p_cover_bootstrap + result[[i]]$p_cover_bootstrap
        p_cover_bootstrap_pre_pivoting = p_cover_bootstrap_pre_pivoting + result[[i]]$p_cover_bootstrap_pre_pivoting
        # APE
        all_APE[i, ] = c(result[[i]]$APE_MLE - result[[i]]$APE_true, result[[i]]$APE_mean - result[[i]]$APE_true, result[[i]]$APE_median - result[[i]]$APE_true, result[[i]]$APE_jackknife - result[[i]]$APE_true, result[[i]]$APE_analytical - result[[i]]$APE_true)
        APE_se_formula[i, ] = c(result[[i]]$APE_se_formula)
        APE_se_formula_own[i, ] = c(result[[i]]$APE_se_formula_own)
        # all bootstrap estimate
        Boot_estimate_all[i,] = result[[i]]$bootstrap_estimate
        T_boot = Boot_estimate_all[i,]- (mean(Boot_estimate_all[i,]) - result[[i]]$MLE_estimate) - result[[i]]$MLE_estimate
        F_hat <- ecdf(T_boot)
        # 4. Prepivot bootstrap replicates
        U_boot <- sapply(T_boot, F_hat)
        # 5. Get CI bounds in uniform space
        U_q <- quantile(U_boot, probs = c(0.025, 0.975))
        CI <- quantile(T_boot, probs = U_q)
        if (CI[1] <= beta && beta <= CI[2]){
          p_pivoting = p_pivoting + 1
        }

        # T_boot <- Boot_estimate_all[i,]- (mean(Boot_estimate_all[i,]) - result[[i]]$MLE_estimate) - result[[i]]$MLE_estimate
        # u_obs <- mean(T_boot <= result[[i]]$MLE_estimate)
        # U_boot <- sapply(1:bootstrap_time, function(i) mean(T_boot <= T_boot[i]))
        # u_bounds <- quantile(U_boot, probs = c(0.025, 0.975))
        # lower <- quantile(T_boot, probs = u_bounds[1])
        # upper <- quantile(T_boot, probs = u_bounds[2])
        # CI <- c(lower, upper)
        # if (CI[1] <= beta && beta <= CI[2]){
        #   p_pivoting = p_pivoting + 1
        # }
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
      cover_pivoting = p_cover_bootstrap_pre_pivoting/length(result)
      cover_pivoting_another = p_pivoting/length(result)
      # APE
      APE_estimate[t_index, ] = apply(all_APE / result[[1]]$APE_true, 2, mean)
      APE_deviation[t_index, ] = apply(all_APE/result[[1]]$APE_true, 2, sd)
      cover_APE = apply( all_APE + result[[1]]$APE_true - 1.96*APE_se_formula_own <= result[[1]]$APE_true & result[[1]]$APE_true <= all_APE + result[[1]]$APE_true + 1.96*APE_se_formula_own, 2, function(x) formatC(mean(x), format = "f", digits = 3))


    }
    t_index = t_index + 1
  }
  end_time <- Sys.time()
  running_time = end_time - start_time

  # save the table
  table_name1 <- paste0('setting_', design, '_bias_estimation_time', running_time ,  ".csv")
  table_name2 <- paste0('setting_', design, '_cover_rate', running_time , ".csv")
  table_name3 <- paste0('setting_', design, 'all_se', running_time , ".csv")
  table_name4 <- paste0('setting_', design, 'all_estimator', running_time , ".csv")
  table_name5 <- paste0('setting_', design, '_cover_rate_APE', running_time , ".csv")
  table_name6 <- paste0('setting_', design, 'all_APE', running_time , ".csv")
  table_name7 <- paste0('setting_', design, 'all_APE_se', running_time , ".csv")
  table_name8 <- paste0('setting_', design, 'APE_bias_estimation', running_time , ".csv")
  table_name9 <- paste0('setting_', design, 'all_APE_se_own', running_time, ".csv")
  write.csv(rbind( Estimate_bias, Estimate_deviation ) , table_name1)
  write.csv(rbind(cover_MLE, cover_bootstrap, cover_jack, cover_analytical,cover_analytical_own, cover_pivoting, cover_pivoting_another) , table_name2)
  write.csv(Estimate_se , table_name3)
  write.csv(all_estimator , table_name4)
  write.csv(cover_APE, table_name5)
  write.csv(all_APE , table_name6)
  write.csv(APE_se_formula, table_name7)
  write.csv(rbind( APE_estimate, APE_deviation) , table_name8)
  write.csv(APE_se_formula_own, table_name9)

}

