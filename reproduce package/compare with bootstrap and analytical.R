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
source('estimation functions.R')
source('necessary functions.R')
source('APE functions for netbootstrap.R')
packages_to_export <- c("dplyr", 'NetBootstrap', 'alpaca', 'speedglm', 'fixest', 'margins')
num_cores <- detectCores()
cl <- makeCluster(num_cores-5)
registerDoParallel(cl)

# probit model network
N = 30
N_seq = c(N)

# bootstrap times
bootstrap_time = 2

# repitition times
mle_num = 2
design = 5

# parameters
beta = c(1, 1)
beta_NULL = c(1, 2)
link = 'probit'
K = length(beta)

for(design in c(1,2,3,4)){
  start_time <- Sys.time()
  t_index = 1
  # empty matrix
  Estimate_bias = matrix(0, length(N_seq), 4*K)
  Estimate_se = matrix(0, length(N_seq), 4*K)
  Estimate_bias_median = matrix(0, length(N_seq), 4*K)
  percentile = matrix(0, length(N_seq), 3*K)
  p_reject_ratio = matrix(0, length(N_seq), 3*K)
  p_rej_ratio_without = matrix(0, length(N_seq), 3*K)
  APE_estimate = matrix(0, length(N_seq), 4*K)
  Estimate_deviation = matrix(0, length(N_seq), 4*K)
  APE_deviation = matrix(0, length(N_seq), 4*K)
  APE_deviation_no_scale = matrix(0, length(N_seq), 4*K)
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
      DGP = data_generation( N = N, beta = beta, design = design, model = link )
      y = DGP$y
      sum(y)
      X_design = DGP$X
      data_in = as.data.frame(DGP$data)
      DGP$eta
      DGP$alpha
      DGP$alpha + DGP$gamma
      pnorm(DGP$eta)
      index_name = colnames(data_in)[(K+2):(K+3)]


      # esalpha# estimation
      fit = network_bootstrap(y = y, X = X_design, N = N, bootstrap_time = bootstrap_time, index = index_name, data = data_in, link = link, boot_repeat = 0)
      #fit_advanced = network_bootstrap_advanced(y = y, X = X_design, N = N, bootstrap_time = bootstrap_time, index = index_name, data = data_in, link = link, boot_repeat = 1)
      fit_analytical = analytical_Amrei(y = y, X = X_design, N = N, index = index_name, data = data_in, link = link, L = 0, beta_NULL = NULL)
      #fit_analyticest_median#fit_analytical = analytical_corrected(y = y, X = X_design, N = N, index = index_name, data = data_in, link = link, L = 0, beta_NULL = NULL)
      #dim(fit_advanced$model_MLE$data)
      #dim(fit_analytical$est_MLE$data)
      # fit$est_critical
      # fit_advanced$est_critical
      # all the estimates
      MLE_estimate = fit$est_MLE[1:K]
      Medain_bootstrap_estimate = fit$est_median[1:K]
      Mean_bootstrap_estimate = fit$est_mean[1:K]
      analytical_estimate = summary(fit_analytical$est)$cm[1:K,1]
      #analytical_estimate = fit_analytical$est[1:K]
      bootstrap_estimate = fit$est_bootstrap_all[,1:K]
      mean_bias = apply(as.matrix(fit$est_bootstrap_all[,1:K]), 2 , mean) - fit$est_MLE[1:K]
      se_MLE = summary(fit_analytical$est_MLE)$cm[1:K,2]
      #se_MLE = fit_analytical$se_MLE
      se_bootstrap = fit$sd
      se_bootstrap_critical = fit$sd_critical
      #se_analytical = fit_analytical$se
      se_analytical = summary(fit_analytical$est)$cm[1:K,2]

      p_cover_bootstrap <- as.integer(
        (MLE_estimate - quantile(fit$est_bootstrap_all[,1:K] - fit$est_critical[1:K], 0.975) <= beta) &
          (beta <= MLE_estimate - quantile(fit$est_bootstrap_all[,1:K] - fit$est_critical[1:K], 0.025))
      )

      p_cover_bootstrap_basic <- as.integer(
        (MLE_estimate - quantile(fit$est_bootstrap_all_critical[,1:K] - MLE_estimate, 0.975) <= beta) &
          (beta <= MLE_estimate - quantile(fit$est_bootstrap_all_critical[,1:K] - MLE_estimate, 0.025))
      )


      p_cover_MLE <- as.integer(
        (MLE_estimate - 1.96*se_MLE <= beta) & (beta <= MLE_estimate + 1.96*se_MLE)
      )


      p_cover_analytical <- as.integer(
        (analytical_estimate - 1.96*se_analytical <= beta) & (beta <= analytical_estimate + 1.96*se_analytical)
      )

      APE_true = array(0, dim = c(N, N, K))
      alpha = DGP$alpha
      gamma = DGP$gamma
      eta = DGP$eta
      X_original = DGP$X_original
      for (l in 1:K){
        for (i in seq(N)){
          for ( j in seq(N)){

            if (i!=j){

              if (length(unique(X_original[,1,l])) == 2){
                max_X = max(unique(X_original[,1,l]))
                min_X = min(unique(X_original[,1,l]))
                if (link == 'probit'){
                  APE_true[i, j, l] = (pnorm(beta[l]*max_X + eta[i,j] - beta[l] * X_original[i,j,l]) - pnorm(beta[l]*min_X + eta[i,j] - beta[l] * X_original[i,j,l]))
                }else{
                  APE_true[i, j, l] = (plogis(beta[l]*max_X + eta[i,j] - beta[l] * X_original[i,j,l]) - plogis(beta[l]*min_X + eta[i,j] - beta[l] * X_original[i,j,l]))
                }
              }else{
                if (link == 'probit'){
                  APE_true[i, j, l] = beta[l] * dnorm(eta[i,j])
                }else{
                  APE_true[i, j, l] = beta[l] * dlogis(eta[i,j])

                }

              }
            }

          }
        }

      }


      APE_true = apply(APE_true, 3, sum )  /(N*(N-1))
      APE_MLE = get_APE_MLE(y = y, X = X_design, N = N, data = data_in, index = index_name, fit = fit, model = link)
      APE_bootstrap = get_APE_bootstrap(y = y, X = X_design, N = N, data = data_in, fit = fit, model = link)
      #APE_bootstrap_advanced = get_APE_bootstrap(y = y, X = X_design, N = N, data = data_in, fit = fit_advanced, model = link)
      #APE_bootstrap_iter = get_APE_bootstrap_iter(y = y, X = X_design, N = N, data = data_in, index = index_name, fit = fit, model = link)
      #APE_bootstrap_iter = get_APE_bootstrap_iter(y = y, X = X_design, N = N, data = data_in, index = index_name, fit = fit_advanced, model = link)
      #APE_bootstrap$APE_median_bootstrap
      #APE_bootstrap_iter$APE_median_bootstrap
      # APE_bootstrap$se
      # APE_analytical = get_APE_analytical_ameri(y = y, X = X_design, N = N, data = data_in, index = index_name, fit = fit_analytical, L = 1, model = link)
      # APE_analytical = get_APE_analytical(y = y, X = X_design, N = N, data = data_in, index = index_name, fit = fit_analytical, L = 1, model = link)
      APE_analytical = getAPEs(fit_analytical$est, panel.structure = 'classic')
      APE_analytical$APE = summary(APE_analytical)$cm[,1]
      APE_analytical$se = summary(APE_analytical)$cm[,2]
      lower <- APE_MLE$APE - apply(APE_bootstrap$APE_bootstrap_all - APE_MLE$APE, 2, quantile, 0.975)
      upper <- APE_MLE$APE - apply(APE_bootstrap$APE_bootstrap_all - APE_MLE$APE, 2, quantile, 0.025)
      cover_APE_bootstrap <- as.integer((lower <= APE_true) & (APE_true <= upper))

      APE_se_formula = c(APE_MLE$se, APE_bootstrap$se, APE_bootstrap$se, APE_analytical$se)


      results = list( MLE_estimate = MLE_estimate, Mean_bootstrap_estimate  = Mean_bootstrap_estimate,
                      Medain_bootstrap_estimate = Medain_bootstrap_estimate,
                      bootstrap_estimate = bootstrap_estimate, analytical_estimate = analytical_estimate,
                      se_MLE = se_MLE, se_bootstrap = se_bootstrap, se_bootstrap_critical = se_bootstrap_critical,
                      se_analytical = se_analytical,
                      p_cover_MLE = p_cover_MLE, p_cover_bootstrap = p_cover_bootstrap,
                      p_cover_bootstrap_basic = p_cover_bootstrap_basic,
                      p_cover_analytical = p_cover_analytical,
                      APE_true = APE_true, APE_MLE = APE_MLE$APE, APE_mean = APE_bootstrap$APE_median_bootstrap,
                      APE_median = APE_bootstrap$APE_median_bootstrap, APE_analytical =  APE_analytical$APE,
                      APE_se_formula = APE_se_formula, cover_APE_bootstrap = cover_APE_bootstrap
      )
      return(results)
    }



    # save the results
    result <- Filter(Negate(is.null), result)

    if (length(result) == 0){
      print('run again')
      t_index = t_index + 1
    } else {
      all_estimator = matrix(0, length(result), 4*K)
      Boot_estimate_all = c()
      se_estimator = matrix(0, length(result), 4*K)
      p_cover_bootstrap = rep(0, K)
      p_cover_bootstrap_basic = rep(0, K)
      p_cover_MLE = rep(0, K)
      p_cover_jack = rep(0, K)
      p_cover_analytical = rep(0, K)
      p_cover_bootstrap_pre_pivoting = rep(0, K)
      cover_APE_bootstrap = rep(0, K)
      all_APE = matrix(0, length(result), 4*K)
      APE_no_scale = matrix(0, length(result), 4*K)
      APE_se_formula = matrix(0, length(result), 4*K)

      for (i in seq(length(result))) {

        all_estimator[i, ] = c(result[[i]]$MLE_estimate, result[[i]]$Mean_bootstrap_estimate, result[[i]]$Medain_bootstrap_estimate, result[[i]]$analytical_estimate)
        se_estimator[i, ] = c(result[[i]]$se_MLE, result[[i]]$se_bootstrap_critical, result[[i]]$se_bootstrap, result[[i]]$se_analytical)
        p_cover_MLE = p_cover_MLE + result[[i]]$p_cover_MLE
        p_cover_analytical = p_cover_analytical + result[[i]]$p_cover_analytical
        p_cover_bootstrap = p_cover_bootstrap + result[[i]]$p_cover_bootstrap
        p_cover_bootstrap_basic = p_cover_bootstrap_basic + result[[i]]$p_cover_bootstrap_basic

        # APE
        all_APE[i, ] = c(result[[i]]$APE_MLE - result[[i]]$APE_true, result[[i]]$APE_mean - result[[i]]$APE_true, result[[i]]$APE_median - result[[i]]$APE_true, result[[i]]$APE_analytical - result[[i]]$APE_true)
        APE_no_scale[i, ] = c(result[[i]]$APE_MLE, result[[i]]$APE_mean, result[[i]]$APE_median, result[[i]]$APE_analytical)
        APE_se_formula[i, ] = c(result[[i]]$APE_se_formula)
        cover_APE_bootstrap = cover_APE_bootstrap + result[[i]]$cover_APE_bootstrap

        # all bootstrap estimate
        Boot_estimate_all = rbind(Boot_estimate_all, result[[i]]$bootstrap_estimate)
      }

      # beta bias
      Estimate_bias[t_index, ] = apply(all_estimator, 2, mean) - rep(beta, 4)
      Estimate_se[t_index, ] = apply(se_estimator, 2, mean)
      Estimate_deviation[t_index, ] = apply(all_estimator , 2, sd)

      # coverage rate for beta
      cover_MLE = p_cover_MLE/length(result)
      cover_analytical = p_cover_analytical/length(result)
      cover_bootstrap = p_cover_bootstrap/length(result)
      cover_bootstrap_basic = p_cover_bootstrap_basic/length(result)

      # APE
      APE_estimate[t_index, ] = apply(all_APE / rep(result[[1]]$APE_true, 4), 2, mean)
      APE_deviation[t_index, ] = apply(all_APE/ rep(result[[1]]$APE_true, 4), 2, sd)
      APE_deviation_no_scale[t_index, ] = apply(APE_no_scale, 2, sd)
      coverage_matrix <- (all_APE - 1.96 * APE_se_formula <= 0) & (0 <= all_APE + 1.96 * APE_se_formula)
      coverage_rate <- apply(coverage_matrix, 2, mean)  # coverage per column
      cover_APE <- formatC(coverage_rate, format = "f", digits = 3)
      cover_APE[(K+1):(2*K)] = cover_APE_bootstrap/length(result) # bootstrap
      cover_APE[(2*K+1):(3*K)] = cover_APE_bootstrap/length(result) # bootstrap
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
  cat("=== Estimate_bias & Estimate_deviation ===\n")
  print(rbind(Estimate_bias, Estimate_deviation))
  write.csv(rbind( Estimate_bias, Estimate_deviation ) , table_name1)
  cat("\n=== cover_MLE, cover_bootstrap_basic, cover_bootstrap, cover_analytical ===\n")
  print(rbind(cover_MLE, cover_bootstrap_basic, cover_bootstrap, cover_analytical))
  write.csv(rbind(cover_MLE, cover_bootstrap_basic, cover_bootstrap, cover_analytical) , table_name2)
  write.csv(Estimate_se , table_name3)
  write.csv(all_estimator , table_name4)
  write.csv(cover_APE, table_name5)
  write.csv(all_APE , table_name6)
  write.csv(APE_se_formula, table_name7)
  cat("\n=== APE_estimate & APE_deviation ===\n")
  print(rbind(APE_estimate, APE_deviation))
  cat("\n=== cover_APE ===\n")
  print(cover_APE)
  write.csv(rbind( APE_estimate, APE_deviation) , table_name8)
  write.csv( rbind(APE_deviation_no_scale, apply(APE_se_formula, 2, mean)) , table_name9)
  write.csv(se_estimator, table_name10)
  write.csv(Boot_estimate_all, table_name11)

}

