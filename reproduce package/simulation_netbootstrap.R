rm(list = ls())
library(NetBootstrap)
library(doParallel)
library (parallel)
library(foreach)
source('simulation data.R')
packages_to_export <- c("dplyr", 'NetBootstrap')
num_cores <- detectCores()
cl <- makeCluster(num_cores-7)
registerDoParallel(cl)

# probit model network
N = 50
N_seq = c(N)

# bootstrap times
bootstrap_time = 499

# repitition times
mle_num = 500

for(design in c(1,2,3,4)){

  # empty matrix
  Estimate_bias = matrix(0, length(N_seq), 3)
  Estimate_bias_median = matrix(0, length(N_seq), 3)
  percentile = matrix(0, length(N_seq), 3)
  p_reject_ratio = matrix(0, length(N_seq), 3)
  p_rej_ratio_without = matrix(0, length(N_seq), 3)
  APE_estimate = matrix(0, length(N_seq), 3)
  Estimate_deviation_boot = matrix(0, length(N_seq), 1)
  Estimate_deviation = matrix(0, length(N_seq), 3)
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

for (N in N_seq) {
  result = foreach(k = 1:mle_num,  .packages = packages_to_export) %dopar% {
    set.seed(k)
    cof = NULL
    cof_j = NULL
    cof_constrain = NULL
    cof_constrain_j = NULL
    estimate_jack = NULL

    # prepare the data
    DGP = data_generation(N = N, beta = beta, design = design)
    y = DGP$y
    X_design = DGP$X
    data_in = DGP$data
    index_name = colnames(data_in)[(K+2):(K+3)]

    # estimation
    fit = network_bootstrap(y = y, X = X_design, N = N, bootstrap_time = bootstrap_time, index = index_name, data = data_in, link = 'probit', beta_NULL = beta_NULL)

    # all the estimates
    MLE_estimate = fit$cof_MLE[1]
    Medain_bootstrap_estimate = fit$cof_median[1]
    Mean_bootstrap_estimate = fit$cof_mean[1]
    bootstrap_estimate = fit$cof_bootstrap_all[,1]
    mean_bias = mean(fit$cof_bootstrap_all[,1]) - fit$cof_MLE[1]
    se_MLE = sqrt(fit$Hessian_MLE[1,1])
    se_bootstrap = fit$sd
    ratio_MLE = 2 * (fit$log_likelihood_MLE - fit$log_likelihood_MLE_NULL)
    ratio_bootstrap = 2 * (fit$log_likelihood_Bootstrap - fit$log_likelihood_Bootstrap_NULL)


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

    # likelihood ratio test
    if (ratio_MLE >= quantile(ratio_bootstrap, 0.95) ){
      p_rej_ratio_bootstrap = 1
    }else{
      p_rej_ratio_bootstrap = 0
    }

    non_central_parameter = (mean_bias) * (se_bootstrap^2)^(-1) * (mean_bias)
    quantile_value <- qchisq(0.95, 1, ncp = non_central_parameter)
    if (ratio_MLE >= quantile_value ){
      p_rej_ratio_without = 1
    }else{
      p_rej_ratio_without = 0
    }

    # APE
    APE_true = matrix(0, nrow = N, ncol = N)
    APE_jack = matrix(0, nrow = N, ncol = N)
    APE_MLE = matrix(0, nrow = N, ncol = N)
    APE_mean = matrix(0, nrow = N, ncol = N)
    APE_median = matrix(0, nrow = N, ncol = N)
    alpha = DGP$alpha
    gamma = DGP$gamma

    # calculate APE
    for (i in seq(N)){
      for ( j in seq(T)){
        APE_true[i, j] = pnorm(beta + alpha[i] + gamma[j]) - pnorm(-beta + alpha[i] + gamma[j])/2
        APE_MLE[i, j] = pnorm(MLE_estimate  + alpha[i] + gamma[j]) - pnorm(-MLE_estimate  + alpha[i] + gamma[j])/2
        APE_mean[i, j] = pnorm(Mean_bootstrap_estimate  + alpha[i] + gamma[j]) - pnorm(-Mean_bootstrap_estimate  + alpha[i] + gamma[j])/2
        APE_median[i, j] = pnorm(Medain_bootstrap_estimate  + alpha[i] + gamma[j]) - pnorm(-Medain_bootstrap_estimate  + alpha[i] + gamma[j])/2
      }
    }


    APE_se_formula = c()

    Z = DGP$X_original
    X_in = Z[,,1]


    # calculate the standard error of the APE based on formula
    APE_se_formula = cbind(APE_se_formula, APE_se( fit = fit, est = fit$cof_MLE, N = N, X = X_in, y = DGP$Y_original, APE =  APE_MLE ))

    APE_se_formula = cbind(APE_se_formula, APE_se( fit = fit, est = fit$cof_mean, N = N, X = X_in, y = DGP$Y_original, APE =  APE_mean ))

    APE_se_formula = cbind(APE_se_formula, APE_se( fit = fit, est = fit$cof_median, N = N, X = X_in, y = DGP$Y_original, APE =  APE_median ))



    results = list( MLE_estimate = MLE_estimate, Mean_bootstrap_estimate  = Mean_bootstrap_estimate , Medain_bootstrap_estimate = Medain_bootstrap_estimate, bootstrap_estimate = bootstrap_estimate,
                    bootstrap_estimate_low = quantile(bootstrap_estimate - MLE_estimate, 0.025), bootstrap_estimate_up = quantile(bootstrap_estimate - MLE_estimate, 0.975),
                    cov_var_MLE = fit$Hessian_MLE[1,1], cov_var_jack = fit$Hessian_MLE[1,1],
                    APE_true = mean(APE_true[!row(APE_true) == col(APE_true)]), APE_jack = mean(APE_jack[!row(APE_jack) == col(APE_jack)]), APE_MLE = mean(APE_MLE[!row(APE_MLE) == col(APE_MLE)]), APE_mean = mean(APE_mean[!row(APE_mean) == col(APE_mean)]), APE_median = mean(APE_median[!row(APE_median) == col(APE_median)]),
                    ratio_MLE = ratio_MLE, ratio_bootstrap = ratio_bootstrap,
                    p_cover_MLE = p_cover_MLE, p_cover_bootstrap = p_cover_bootstrap, p_cover_jack = p_cover_MLE,
                    p_rej_ratio_bootstrap = p_rej_ratio_bootstrap, p_rej_ratio_without = p_rej_ratio_without,
                    APE_se_formula = APE_se_formula)
    return(results)
  }



  # svae the results
  result <- Filter(Negate(is.null), result)

  if (length(result) == 0){
    print('run again')
    t_index = t_index + 1
  } else {
    all_estimator = matrix(0, length(result), 3)
    all_APE = matrix(0, length(result), 3)
    Boot_estimate = matrix(0, length(result), 2)
    cov_var_MLE = matrix(0, length(result), 1)
    p_cover_bootstrap = 0
    p_cover_MLE = 0
    p_rej_ratio_bootstrap = 0
    p_rej_ratio_without = 0
    APE_se_formula = matrix(0, length(result), 3)
    APE_MLE = matrix(0, length(result), 1)

    for (i in seq(length(result))) {
      all_estimator[i, ] = c(result[[i]]$MLE_estimate, result[[i]]$Mean_bootstrap_estimate, result[[i]]$Medain_bootstrap_estimate)
      all_APE[i, ] = c(result[[i]]$APE_MLE - result[[i]]$APE_true, result[[i]]$APE_mean - result[[i]]$APE_true, result[[i]]$APE_median - result[[i]]$APE_true)
      cov_var_MLE[i, ] = c(result[[i]]$cov_var_MLE)
      p_cover_MLE = p_cover_MLE + result[[i]]$p_cover_MLE
      p_cover_bootstrap = p_cover_bootstrap + result[[i]]$p_cover_bootstrap
      p_rej_ratio_bootstrap = p_rej_ratio_bootstrap + result[[i]]$p_rej_ratio_bootstrap
      p_rej_ratio_without = p_rej_ratio_without + result[[i]]$p_rej_ratio_without
      APE_se_formula[i,] = c(result[[i]]$APE_se_formula)
      APE_MLE[i] = result[[i]]$APE_MLE

    }

    # beta bias
    Estimate_bias[t_index, ] = apply(all_estimator - beta, 2, mean)
    APE_estimate[t_index, ] = apply(all_APE/ result[[1]]$APE_true, 2, mean)
    Estimate_deviation[t_index, ] = apply(all_estimator , 2, sd)
    APE_deviation[t_index, ] = apply(all_APE/result[[1]]$APE_true, 2, sd)
    APE_se_formula = colMeans(APE_se_formula)
    APE_sd = sd(APE_MLE)
    cov_MLE[t_index, ] = apply(cov_var_MLE, 2, mean)
    cov_jack[t_index, ] = apply(cov_var_jack, 2, mean)

    # coverage rate for beta
    cover_MLE = p_cover_MLE/length(result)
    cover_jack = p_cover_jack/length(result)
    cover_bootstrap = p_cover_bootstrap/length(result)

    # rejection rate for APE
    cover_APE = apply( all_APE + result[[1]]$APE_true - rep(1.96*APE_se_formula, length(result)) <= result[[1]]$APE_true & result[[1]]$APE_true <= all_APE + result[[1]]$APE_true + rep(1.96*APE_se_formula, length(result)), 2, function(x) formatC(mean(x), format = "f", digits = 3))

    # LR ratio by parametric bootstrap
    p_reject_ratio = p_rej_ratio_bootstrap/length(result)
    p_rej_ratio_without = p_rej_ratio_without/length(result)
  }
  t_index = t_index + 1
}

# hist plot
hist(all_estimator[,1])

# print out to check
print(Estimate_bias)
print(Estimate_deviation)
print(cov_MLE)
print(APE_estimate)
print(cover_APE)

# save the table
table_name1 <- paste0('setting_', design, '_bias_estimation', ".csv")
table_name2 <- paste0('setting_', design, '_cover_rate', ".csv")
table_name3 <- paste0('setting_', design, '_test_power', ".csv")
table_name4 <- paste0('setting_', design, '_cover_rate_APE', ".csv")
write.csv(rbind( Estimate_bias, Estimate_deviation, APE_estimate, APE_deviation ) , table_name1)
write.csv(rbind(cover_MLE, cover_bootstrap) , table_name2)
write.csv(rbind(p_reject_ratio, p_rej_ratio_without), table_name3)
write.csv(cover_APE, table_name4)
}


