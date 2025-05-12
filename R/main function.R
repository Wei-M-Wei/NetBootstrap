#' Parametric bootstrap inference in network formation models
#'@title Parametric bootstrap inference in network formation models
#'
#' @description R package 'NetBootstrap' is dedicated to do parametric bootstrap inference in network formation models.
#'
#' @param y outcome variabe
#' @param X dependent variables
#' @param N size of the individuals
#' @param bootstrap_time number of the bootstrap times
#' @param data data which you want to apply
#' @param index name of the import side and export side
#' @param link probit or logit model
#' @param beta_NULL Null hyphothesis of the first covariate
#'
#' @returns A list of fitted results is returned.
#' Within this outputted list, the following elements can be found:
#'     \item{cof_mean}{estimated coefficients corrected by the mean of the bootstrap estimates.}
#'     \item{cof_median}{estimated coefficients corrected by the median of the bootstrap estimates.}
#'     \item{sd}{stand deviation form the bootstrap procedure.}
#'     \item{cof_MLE}{estimated coefficients from Maximum Likelihood Estimates (MLE).}
#'     \item{cof_bootstrap_all}{all bootstrap estimates.}
#'     \item{cof_MLE_NULL}{Maximum Likelihood Estimates (MLE) of coefficients under the null hyphothesis.}
#'     \item{cof_bootstrap_NULL}{Coefficient estimates from the bootstrap under the null hyphothesis.}
#'     \item{log_likelihood_MLE}{Log-likelihood of the full model based on Maximum Likelihood Estimates (MLE).}
#'     \item{log_likelihood_Bootstrap}{Log-likelihood of the full model using bootstrap estimates.}
#'     \item{log_likelihood_MLE_NULL}{Log-likelihood of under the null hyphothesis based on MLE.}
#'     \item{log_likelihood_Bootstrap_NULL}{Log-likelihood under the null hyphothesis using bootstrap estimates.}
#'     \item{Hessian_MLE}{Inverse Hessian matrix for the full model (used for variance estimation).}
#'     \item{Hessian_MLE_NULL}{Inverse Hessian matrix under the null hyphothesis.}
#'     \item{X_origin}{Original design matrix used in the model, as a numeric matrix, which contains the dummy variable.}
#'
#' @import speedglm
#' @export
#'
#'@examples
#' rm(list = ls())
#' library(NetBootstrap)

#' generate the simulated data
#' N = 30
#' K = 1
#' beta = 1
#' alpha = rep(0, N)
#' gamma = rep(0, N)
#' dens_N_up = log(log(N))
#' dens_N_low = -log(log(N))
#' dens_T_up = log(log(N))
#' dens_T_low = -log(log(N))
#' for (i in seq(N)){
#'   alpha[i] = dens_N_low + ((i - 1) / ( N - 1)) * (dens_N_up - dens_N_low)
#' }
#' for (i in seq(N)){
#'   gamma[i] = dens_T_low + ((i - 1) / ( N - 1)) * (dens_T_up - dens_T_low)
#' }
#'
#' Z = array(0, dim = c(N, N, K))
#' Y = matrix(0, nrow = N, ncol = N)
#' for (t in seq(N)) {
#'   for (i in seq(N)) {
#'     v_it = rnorm(1, 0, 1 / sqrt(2))
#'     epsi_it = rnorm(1, 0, 1)
#'     if ( t %% 2 == 0 & i %% 2 == 0){
#'       Z[i, t, ] = (1 - 2 * 1*(0)) * (1 - 2 * 1*(0))
#'     } else if(t %% 2 == 0 & i %% 2 == 1){
#'       Z[i, t, ] = (1 - 2 * 1*(1)) * (1 - 2 * 1*(0))
#'     } else if(t %% 2 == 1 & i %% 2 == 0){
#'       Z[i, t, ] = (1 - 2 * 1*(0)) * (1 - 2 * 1*(1))
#'     } else{
#'       Z[i, t, ] = (1 - 2 * 1*(1)) * (1 - 2 * 1*(1))
#'     }
#'     Y[i, t] = (beta * Z[i, t, ] + alpha[i] + gamma[t] > epsi_it)
#'
#'   }
#' }
#'
#'
#' # final data
#' y = c(Y[!row(Y) == col(Y)])
#' X_in = matrix(0, N * N, K )
#' index_matrix = matrix(0, N * N, 2)
#' for (t in seq(N)) {
#'   for (i in seq(N)) {
#'     alpha_in = rep(0, N)
#'     gamma_in = rep(0, N)
#'     alpha_in[i] = 1
#'     gamma_in[t] = 1
#'     X_in[i + (t - 1) * N,] = c(Z[i,t,])
#'     index_matrix[i + (t - 1) * N,] = c(i,t)
#'   }
#' }
#' drop_index = NULL
#' for (i in seq(N)){
#'   drop_index = cbind(drop_index, i + (i - 1) * N)
#' }
#' X_design = as.matrix(X_in[-drop_index,])
#' index_in = index_matrix[-drop_index,]
#'
#' # estimation procedure
#' data_in = data.frame(y = y, X = X_design, index = index_in)
#' index_name = colnames(data_in)[(K+2):(K+3)]
#' fit = network_bootstrap(y, X = X_design, N, bootstrap_time = 10, index = c('index.1', 'index.2'), data = data_in, link = 'probit', beta_NULL = NULL)
#'
#' # get the results
#' est_MLE = fit$cof_MLE[1]
#' est_corrected = fit$cof[1]
#' sd = fit$sd

#'
network_bootstrap = function(y, X, N, bootstrap_time, index, data, link = 'probit', beta_NULL = NULL){

  data = data.frame(y = y, X = X, data[,index[1]], data[,index[2]])
  K = dim(X)[2]
  # order the data
  if(is.null(colnames(X)) == 1){
    colnames(data)[(dim(data)[2]-1):dim(data)[2]] = index
    data = data[order(data[,index[2]], data[,index[1]]),]
    p = dim(X)[2]
    y = data$y
    if (p==1){
      X = data[,2]
    }else{
      X = data[,(2:p+1)]
    }
  }else{
    colnames(data) = c('y', colnames(X), index)
    data = data[order(data[,index[2]], data[,index[1]]),]
    p = dim(X)[2]

    y = data$y
    X = data[,colnames(X)]
  }

  # add dummy variable
  fix_effect = matrix(0, N*N, N + N )
  for (t in seq(N)) {
    for (i in seq(N)) {
      alpha_in = rep(0, N)
      gamma_in = rep(0, N)
      alpha_in[i] = 1
      gamma_in[t] = 1
      fix_effect[i + (t - 1) * N,] = c(alpha_in - alpha_in[1], gamma_in + alpha_in[1])
    }
  }
  drop_index = NULL
  for (i in seq(N)){
    drop_index = cbind(drop_index, i + (i - 1) * N)
  }
  fix = fix_effect[-drop_index,]

  # prepare teh final data
  X_design = cbind(X, fix)
  X_design = apply(X_design, 2, as.numeric)
  data_in <- data.frame(y = y, X = X_design)

  if (sum(y) <= N*N*0.01){
    return(NULL)
  }

  # MLE
  model <-
    speedglm(y ~ . - 1, data = data_in, family = binomial(link = link))
  fit = summary(model)
  Hessian_inv = vcov(model)
  cof = unlist(as.list(fit$coefficients[, 1]))
  log_likelihood_estimate <- logLik(model)
  cof[K+1] = sum(cof[(N+K+1):(N+N+K)]) - sum(cof[(K+2):(N+K)])

  # constrained MLE
  data_2 = data_in
  if(is.null(beta_NULL) != 1){

    formula <- as.formula( paste("y ~ -1 +", paste(colnames(data_2[,-2])[-1], collapse = " + "), "+ offset(offset_term)"))
    data_2$offset_term <- beta_NULL * X_design[,1]
    model_NULL <-
      glm(formula = formula, data = data_2, family=binomial(link = link))
    fit_NULL = summary(model_NULL)
    Hessian_inv_NULL = vcov(model_NULL)
    cof_NULL =  c(beta_NULL, NA, unlist(as.list(fit_NULL$coefficients[, 1])))
    log_likelihood_estimate_NULL <- logLik(model_NULL)

  }

  # prepare for bootstrap
  cof_B = NULL
  log_likelihood_estimate_B = NULL
  cof_B_NULL = NULL
  log_likelihood_estimate_B_NULL = NULL
  X_design = apply(X_design, 2, as.numeric)

  # bootstrap procedure
  for (k in seq(bootstrap_time)) {
    set.seed(k)
    epsi_it = rnorm(length(y), 0, 1)
    Y = as.numeric(X_design %*% cof > epsi_it)
    data_boostrap <- data.frame(y = Y, X = X_design)
    model_B <-
      speedglm(y ~ . - 1, data = data_boostrap, family = binomial(link = link))
    fit_B = summary(model_B)
    cof_B = rbind(cof_B, unlist(as.list(fit_B$coefficients[, 1])))
    log_likelihood_estimate_B <- rbind(log_likelihood_estimate_B, logLik(model_B))


    # constrained
    data_2 = data_boostrap
    if(is.null(beta_NULL) != 1){
      formula <- as.formula( paste("y ~ -1 +", paste(colnames(data_2[,-2])[-1], collapse = " + "), "+ offset(offset_term)"))
      data_2$offset_term <- cof[1] * X_design[,1]
      model_B_NULL <-
        speedglm(formula = formula, data = data_2, family=binomial(link = link))
      fit_B_NULL = summary(model_B_NULL)
      cof_B_NULL = rbind(cof_B_NULL, c(beta_NULL, NA, unlist(as.list(fit_B_NULL$coefficients[, 1]))))
      log_likelihood_estimate_B_NULL <- rbind(log_likelihood_estimate_B_NULL, logLik(model_B_NULL))
    }
  }

  # get the final results
  cof_B[,K+1] = apply(cof_B[,(N+K+1):(N+N+K)],1,sum) - apply(cof_B[,(K+2):(N+K)],1,sum)
  eta = cof_B %*% t(X_design)
  eta_MLE = cof %*% t(X_design)
  cof_boost_mean = apply(cof_B, 2, mean)
  cof_boost_median = apply(cof_B, 2, median)
  est_correct_mean = cof - (cof_boost_mean - cof)
  est_correct_median = cof - (cof_boost_median - cof)
  est_correct_mean[K+1] = sum(est_correct_mean[(N+K+1):(N+N+K)]) - sum(est_correct_mean[(K+2):(N+K)])
  est_correct_median[K+1] = sum(est_correct_median[(N+K+1):(N+N+K)]) - sum(est_correct_median[(K+2):(N+K)])
  if (K == 1){
    boostrap_sd = sd(cof_B[,1:K])
  }
  else{
    boostrap_sd = apply(cof_B[,1:K],2,sd)
  }


  if(is.null(beta_NULL) != 1){
    res = list(est_MLE = cof, est_mean = est_correct_mean, est_median = est_correct_median, sd = boostrap_sd,
               est_bootstrap_all = cof_B, cof_MLE_NULL = cof_NULL, cof_bootstrap_NULL = cof_B_NULL,
               log_likelihood_MLE = log_likelihood_estimate, log_likelihood_Bootstrap = log_likelihood_estimate_B,
               log_likelihood_MLE_NULL = log_likelihood_estimate_NULL, log_likelihood_Bootstrap_NULL = log_likelihood_estimate_B_NULL,
               Hessian_MLE = Hessian_inv, Hessian_MLE_NULL = Hessian_inv_NULL, X_origin = as.matrix(X_design), eta = eta, eta_MLE = eta_MLE, data = data
    )
  }
  else{
    res = list(est_MLE = cof, est_mean = est_correct_mean, est_median = est_correct_median, sd = boostrap_sd, est_bootstrap_all = cof_B,
               log_likelihood_MLE = log_likelihood_estimate, log_likelihood_Bootstrap = log_likelihood_estimate_B,
               Hessian_MLE = Hessian_inv, X_origin = as.matrix(X_design), eta = eta, eta_MLE = eta_MLE, data = data
    )
  }
  return(res)
}


#' @export
split_jackknife = function(y, X, N, index, data, link = 'probit', beta_NULL = NULL){

  ##################
  # part of estimate
  ##################

  data_j = data.frame(y = y, X = X, data[,index[1]], data[,index[2]])
  K = dim(X)[2]
  # order the data
  if(is.null(colnames(X)) == 1){
    colnames(data_j)[(dim(data_j)[2]-1):dim(data_j)[2]] = index
    data_j = data_j[order(data_j[,index[2]], data_j[,index[1]]),]
    p = dim(X)[2]
    y = data_j$y
    if (p==1){
      X = data_j[,2]
    }else{
      X = data_j[,(2:p+1)]
    }
  }else{
    colnames(data_j) = c('y', colnames(X), index)
    data_j = data_j[order(data_j[,index[2]], data_j[,index[1]]),]
    p = dim(X)[2]

    y = data_j$y
    X = data_j[,colnames(X)]
  }

  # add dummy variable
  fix_effect = matrix(0, N*N, N + N )
  for (t in seq(N)) {
    for (i in seq(N)) {
      alpha_in = rep(0, N)
      gamma_in = rep(0, N)
      alpha_in[i] = 1
      gamma_in[t] = 1
      fix_effect[i + (t - 1) * N,] = c(alpha_in - alpha_in[1], gamma_in + alpha_in[1])
    }
  }
  drop_index = NULL
  for (i in seq(N)){
    drop_index = cbind(drop_index, i + (i - 1) * N)
  }
  fix = fix_effect[-drop_index,]
  # prepare teh final data
  X_design = cbind(X, fix)
  X_design = apply(X_design, 2, as.numeric)
  X_save = X_design
  data_in <- data.frame(y = y, X = X_design)

  # MLE
  model <-
    speedglm(y ~ . - 1, data = data_in, family = binomial(link = link))
  fit = summary(model)
  Hessian_inv = vcov(model)
  cof = unlist(as.list(fit$coefficients[, 1]))
  log_likelihood_estimate <- logLik(model)
  cof[K+1] = sum(cof[(N+K+1):(N+N+K)]) - sum(cof[(K+2):(N+K)])


  cof_j = NULL
  cof_constrain_j = NULL
  log_likelihood_j = NULL
  log_likelihood_constrain_j = NULL

  index1 = data_j[,index[1]]
  index2 = data_j[,index[2]]
  new_index1 <- match(index1, sort(unique(index1)))
  new_index2 <- match(index2, sort(unique(index2)))

  data_j[,index[1]] =  new_index1
  data_j[,index[2]] =  new_index2
  data = data_j

  for (k in seq(N-1) ){
    index_jack = list()
    for (i in seq(N-k)){
      index_jack[[i]] =   c(i,i+k)
    }
    for (j in seq(N-k+1,N)){
      index_jack[[j]] = c(j,j-N+k)
    }

    match_idx <- sapply(index_jack, function(pair) {
      which(data[,index[2]] == pair[1] & data[,index[1]] == pair[2])
    })

    y_split = y[-match_idx]
    X_split = as.matrix(X)[-match_idx,]

    # add dummy variable
    fix_effect = matrix(0, N*N, N + N )
    for (t in seq(N)) {
      for (i in seq(N)) {
        alpha_in = rep(0, N)
        gamma_in = rep(0, N)
        alpha_in[i] = 1
        gamma_in[t] = 1
        fix_effect[i + (t - 1) * N,] = c(alpha_in - alpha_in[1], gamma_in + alpha_in[1])
      }
    }
    drop_index = NULL
    for (i in seq(N)){
      drop_index = cbind(drop_index, i + (i - 1) * N)
    }
    fix = fix_effect[-drop_index,]

    # prepare the final data
    X_design = cbind(X_split, fix[-match_idx,])
    X_design = apply(X_design, 2, as.numeric)

    data_in <- data.frame(y = y_split, X = X_design)


    model_j <-
      speedglm(y ~ . - 1, data = data_in, family = binomial(link = "probit"))
    fit_j = summary(model_j)
    Hessian_inv_1_j = vcov(model_j)
    cof_j = rbind(cof_j, unlist(as.list(fit_j$coefficients[, 1])))
    log_likelihood_j <- rbind( log_likelihood_j,logLik(model_j))

    if(is.null(beta_NULL) != 1){
      ## Constrainted MLE
      data_2 <- data_in
      formula <- as.formula( paste("y ~ -1 +", paste(colnames(data_2[,-2])[-1], collapse = " + "), "+ offset(offset_term)"))
      data_2$offset_term <- beta_NULL * X_design[,1]
      model_j_2 <-
        glm(formula = formula, data = data_2, family=binomial(link = 'probit'))
      fit_j_2 = summary(model_j_2)
      Hessian_inv_2_constrain_j = vcov(model_j_2)
      cof_constrain_j = rbind(cof_constrain_j, c(beta_NULL, NA, unlist(as.list(fit_j_2$coefficients[, 1]))))
      log_likelihood_constrain_j <- rbind( log_likelihood_constrain_j  , logLik(model_j_2))
    }
  }

  estimate_jack = (N-1)*cof - (N-2) * apply(cof_j, 2, mean)
  estimate_jack[K+1] = sum(estimate_jack[(N+K+1):(N+N+K)]) - sum(estimate_jack[(K+2):(N+K)])
  cof_j[,K+1] = apply(cof_j[,(N+K+1):(N+N+K)],1,sum) - apply(cof_j[,(K+2):(N+K)],1,sum)
  cof_jack_all = cof_j


  eta = cof_jack_all %*% t(X_save)
  eta_MLE = cof %*% t(X_save)


  ################
  # standard error
  ################

  # use the MLE estimator
  # the use of bias-corrected estimator is similar to the biased one
  # calculate the X'beta + pi, formed as a matrix
  cov_sum_1 = X_save[,1] * cof[1]
  cov_sum_2 = X_save[,-1] %*% cof[-1]
  cov_sum = matrix(cov_sum_1 + cov_sum_2, N-1, N)
  cov_sum = shift_lower_triangle_and_add_zero_diag(cov_sum)
  cov_sum[row(cov_sum) == col(cov_sum)] = 0
  X_to = as.matrix(X)
  y_to = y
  se = rep(0,K)
  se_huges = rep(0,K)
  se_no_MLE = rep(0, K)
  se_huges_no_MLE = rep(0, K)
  est_another = rep(0, K)
  se_corrected = rep(0, K)

  for (index_covariate in 1:K) {

    # X is a matrix of a single covariate, and the same for y
    X = vector_to_matrix(X_to[,index_covariate], N, ind1 = data_j[,index[1]], ind2 = data_j[,index[2]])
    y = vector_to_matrix(y_to, N, ind1 = data_j[,index[1]], ind2 = data_j[,index[2]])
    y[row(y) == col(y)] = 0

    # CDF (Φ(Xβ)) and PDF (φ(Xβ))
    Phi_XB <- pnorm(cov_sum)
    phi_XB <- dnorm(cov_sum)
    dd_F_fix = -cov_sum * phi_XB
    Phi_XB[row(Phi_XB) == col(Phi_XB)] = 0
    phi_XB[row(phi_XB) == col(phi_XB)] = 0
    Phi_XB <- pmax(Phi_XB, 1e-9)
    Phi_XB <- pmin(Phi_XB, 1 - 1e-9)

    derivative_ingredients = compute_derivatives(eta = cov_sum, y = y, X = X)

    # preparation for ingredients
    d_fix_loss =  derivative_ingredients$d_fix_loss
    d_beta_loss = derivative_ingredients$d_beta_loss
    d_fix_fix_loss = derivative_ingredients$d_fix_fix_loss# d_beta_loss = X * (y - exp(cov_sum)/( 1 + exp(cov_sum)))
    d_beta_beta_loss = derivative_ingredients$d_beta_beta_loss# d_beta_fix_loss = X * ( - exp(cov_sum)/( 1 + exp(cov_sum))^2)
    d_beta_fix_loss = derivative_ingredients$d_beta_fix_loss

    # one way to calculate the se by formula
    # Hessian matrix
    H_a_a = matrix(0, N, N)
    H_g_g = matrix(0, N, N)
    for (i in 1:N) {
      H_a_a[i, i] <- -(sum(d_fix_fix_loss[i, ]) - d_fix_fix_loss[i, i])/(N-1)  # sum of off-diagonal elements in row i
    }
    for (i in 1:N) {
      H_g_g[i, i] <- -(sum(d_fix_fix_loss[-i, i]))/(N-1) # sum of off-diagonal elements in row i
    }
    H_a_g = -d_fix_fix_loss/((N-1))

    Hessian_bar =   cbind(rbind(H_a_a, t(H_a_g)), rbind(H_a_g, H_g_g))  + c(rep(1,N), rep(-1,N)) %*% t( c(rep(1,N), rep(-1,N)) )/N
    Hessain_inverse = solve(Hessian_bar)
    Hessian_a_a = Hessain_inverse[1:(N), 1:(N)]
    Hessian_g_a = Hessain_inverse[(N+1):(N+N), 1:(N)]
    Hessian_a_g = Hessain_inverse[1:(N), (N+1):(N+N)]
    Hessian_g_g = Hessain_inverse[(N+1):(N+N), (N+1):(N+N)]

    # matrix 'Xi'
    the = matrix(0, N, N)
    for (i in 1:N) {
      for (j in 1:N) {
        temp_sum <- 0
        for (s in 1:N) {
          for (t in 1:N) {
            if (t != s) {
              temp_sum <- temp_sum + (-1/(N)) * (Hessian_a_a[i, s] + Hessian_g_a[j, s] + Hessian_a_g[i, t] + Hessian_g_g[j, t]) * d_beta_fix_loss[s, t]
            }
          }
        }
        the[i, j] <- temp_sum
      }
    }

    # W_hat based on Iva ́ n2016
    W_hat =   -(1 / ((N-1)*N)) * (sum(  d_beta_beta_loss  - d_fix_fix_loss * the * the  ) - sum(diag(d_beta_beta_loss  - d_fix_fix_loss * the * the )) ) # -(1 / N) * (  d_beta_beta_big_loss  - d_beta_fix_big_loss * solve(d_fix_2_big_loss) * d_beta_fix_big_loss  )



    #### Another way to calculate the bias corrected
    # another way to calculate the analytical corrected estimate
    H = phi_XB/(Phi_XB*(1-Phi_XB))
    small_w = H * phi_XB

    # W_hat
    X_into = matrix_to_panel_df(X)
    weight = matrix_to_panel_df(small_w)$X
    weight[which(weight==0)] = 1
    re = get_weighted_projection_fitted_exclude_t_eq_i(X_into$X, weight, X_into$id, X_into$time)
    re_matrix = matrix(re, N-1, N)
    re_matrix = shift_lower_triangle_and_add_zero_diag(re_matrix)
    tilde_X = X - re_matrix
    W_hat_another =   (1 / ((N-1)*N)) * (sum(  small_w  * tilde_X  * tilde_X   ) - sum(diag( small_w  * tilde_X  * tilde_X  )) ) # -(1 / N) * (  d_beta_beta_big_loss  - d_beta_fix_big_loss * solve(d_fix_2_big_loss) * d_beta_fix_big_loss  )


    # Omega_hat based on Iva ́ n2016
    D_beta_loss = d_beta_loss - d_fix_loss*the
    Omega_hat <- 0
    for (i in 1:N) {
      for (j in setdiff(1:N, i)) {
        for (k in setdiff(1:N, i)) {
          Omega_hat <- Omega_hat + D_beta_loss[i, j] * D_beta_loss[i, k]/(N*(N-1))
        }
      }
    }


    # the third way to calculate se, which is based on huges paper
    # matrix from huges
    d_beta_beta_big_loss = (1/(N-1))*(sum(d_beta_beta_loss) - sum(diag(d_beta_beta_loss)))
    d_beta_fix_big_loss = (1/(N-1))*(sum(d_beta_fix_loss) - sum(diag(d_beta_fix_loss)))
    d_fix_beta_big_loss = (1/(N-1))*(sum(d_beta_fix_loss) - sum(diag(d_beta_fix_loss)))
    d_fix_2_big_loss = (1/(N-1))*(sum(d_fix_fix_loss) - sum(diag(d_fix_fix_loss)))

    # W_hat huges
    W_hat_huges =   -(1 / N) * (  d_beta_beta_big_loss  - d_beta_fix_big_loss * solve(d_fix_2_big_loss) * d_beta_fix_big_loss  )
    D_beta_loss = d_beta_loss - d_fix_loss*the
    A = (D_beta_loss + t(D_beta_loss))^2

    # Omega_hat
    Omega_hat_huges =  1/(N*(N-1))* sum(A[lower.tri(A)])

    se[index_covariate] = sqrt(solve(W_hat)/(N*(N-1)))
    se_no_MLE[index_covariate] = sqrt(solve(W_hat) * Omega_hat * solve(W_hat)/(N*(N-1)))
    se_corrected[index_covariate] = sqrt(solve(W_hat_another)/(N*(N-1)))
    se_huges[index_covariate] = sqrt(solve(W_hat_huges)/(N*(N-1)))
    se_huges_no_MLE[index_covariate] = sqrt(solve(W_hat_huges) * Omega_hat_huges * solve(W_hat_huges)/(N*(N-1)))

  }

  if(is.null(beta_NULL) != 1){
    res = list(est = estimate_jack, est_MLE = cof, est_jackknife_NULL = cof_constrain_j, est_jackknife_all = cof_jack_all,
               log_likelihood_jack = log_likelihood_j,
               log_likelihood_jack_NULL = log_likelihood_constrain_j,
               Hessian_jack = Hessian_inv_1_j, Hessian_jack_NULL = Hessian_inv_2_constrain_j, X_origin = as.matrix(X_save),
               eta = eta, eta_MLE = eta_MLE,
               se_no_MLE = se_no_MLE,
               se = se_corrected, se_formula = se,
               se_huges = se_huges, se_huges_no_MLE = se_huges_no_MLE
    )
  }else{
    res = list(est = estimate_jack, est_MLE = cof, est_jackknife_all = cof_jack_all,
               log_likelihood_jack = log_likelihood_j,
               Hessian_jack = Hessian_inv_1_j, X_origin = as.matrix(X_save),
               eta = eta, eta_MLE = eta_MLE,
               se_no_MLE = se_no_MLE,
               se = se_corrected, se_formula = se,
               se_huges = se_huges, se_huges_no_MLE = se_huges_no_MLE
    )
  }
  return(res)

}


#' @export
#' @import alpaca
analytical_Amrei = function(y, X, N, index, data, link = 'probit', L = 1, beta_NULL = NULL){

  data = data.frame(y = y, X = X, data[,index[1]], data[,index[2]])
  K = dim(X)[2]
  # order the data
  if(is.null(colnames(X)) == 1){
    colnames(data)[(dim(data)[2]-1):dim(data)[2]] = index
    data = data[order(data[,index[2]], data[,index[1]]),]
    p = dim(X)[2]
    y = data$y
    if (p==1){
      X = data[,2]
    }else{
      X = data[,(2:p+1)]
    }
  }else{
    colnames(data) = c('y', colnames(X), index)
    data = data[order(data[,index[2]], data[,index[1]]),]
    p = dim(X)[2]

    y = data$y
    X = data[,colnames(X)]
  }
  mod <- alpaca::feglm(y ~ X | index.1 + index.2, data = data, family = binomial(link))
  return(est = biasCorr(mod, L = L, panel.structure = c( "network")))

}


#' @export
analytical_corrected = function(y, X, N, index, data, link = 'probit', L = 1, beta_NULL = NULL){

  ##################
  # part of estimate
  ##################

  data = data.frame(y = y, X = X, data[,index[1]], data[,index[2]])
  K = dim(X)[2]
  # order the data
  if(is.null(colnames(X)) == 1){
    colnames(data)[(dim(data)[2]-1):dim(data)[2]] = index
    data = data[order(data[,index[2]], data[,index[1]]),]
    p = dim(X)[2]
    y = data$y
    if (p==1){
      X = data[,2]
    }else{
      X = data[,(2:p+1)]
    }
  }else{
    colnames(data) = c('y', colnames(X), index)
    data = data[order(data[,index[2]], data[,index[1]]),]
    p = dim(X)[2]

    y = data$y
    X = data[,colnames(X)]
  }

  # add dummy variable
  fix_effect = matrix(0, N*N, N + N )
  for (t in seq(N)) {
    for (i in seq(N)) {
      alpha_in = rep(0, N)
      gamma_in = rep(0, N)
      alpha_in[i] = 1
      gamma_in[t] = 1
      fix_effect[i + (t - 1) * N,] = c(alpha_in - alpha_in[1], gamma_in + alpha_in[1])
    }
  }
  drop_index = NULL
  for (i in seq(N)){
    drop_index = cbind(drop_index, i + (i - 1) * N)
  }
  fix = fix_effect[-drop_index,]

  # prepare teh final data
  X_design = cbind(X, fix)
  X_design = apply(X_design, 2, as.numeric)
  data_in <- data.frame(y = y, X = X_design)

  if (sum(y) <= N*N*0.01){
    return(NULL)
  }

  # MLE
  model <-
    speedglm(y ~ . - 1, data = data_in, family = binomial(link = link))
  fit = summary(model)
  Hessian_inv = vcov(model)
  cof = unlist(as.list(fit$coefficients[, 1]))
  log_likelihood_estimate <- logLik(model)
  cof[K+1] = sum(cof[(N+K+1):(N+N+K)]) - sum(cof[(K+2):(N+K)])


  ##################
  # se
  ##################

  est = cof
  cov_sum_1 = X_design[,1] * est[1]
  cov_sum_2 = X_design[,-1] %*% est[-1]
  cov_sum = matrix(cov_sum_1 + cov_sum_2, N-1, N)
  cov_sum = shift_lower_triangle_and_add_zero_diag(cov_sum)
  X_to = as.matrix(X)
  y_to = y
  estimate_analytical = rep(0, K)
  estimate_analytical_another = rep(0, K)
  se = rep(0,K)
  se_no_MLE = rep(0, K)
  est_another = rep(0, K)
  se_another = rep(0, K)


  for (index_covariate in 1:K) {

    # X is a matrix of a single covariate, and the same for y
    X = vector_to_matrix(X_to[,index_covariate], N, ind1 = data[,index[1]], ind2 = data[,index[2]])
    y = vector_to_matrix(y_to, N, ind1 = data[,index[1]], ind2 = data[,index[2]])
    y[row(y) == col(y)] = 0


    Phi_XB <- pnorm(cov_sum)  # CDF (Φ(Xβ))
    phi_XB <- dnorm(cov_sum)  # PDF (φ(Xβ))
    dd_F_fix = -cov_sum * phi_XB
    ddd_F_fix = cov_sum^2 * phi_XB - phi_XB
    Phi_XB[row(Phi_XB) == col(Phi_XB)] = 0
    phi_XB[row(phi_XB) == col(phi_XB)] = 0
    Phi_XB <- pmax(Phi_XB, 1e-9)
    Phi_XB <- pmin(Phi_XB, 1 - 1e-9)

    derivative_ingredients = compute_derivatives(eta = cov_sum, y = y, X = X)

    # preparation for ingredients
    d_fix_loss =  derivative_ingredients$d_fix_loss
    d_beta_loss = derivative_ingredients$d_beta_loss
    d_fix_fix_loss = derivative_ingredients$d_fix_fix_loss# d_beta_loss = X * (y - exp(cov_sum)/( 1 + exp(cov_sum)))
    d_beta_beta_loss = derivative_ingredients$d_beta_beta_loss# d_beta_fix_loss = X * ( - exp(cov_sum)/( 1 + exp(cov_sum))^2)
    d_beta_fix_loss = derivative_ingredients$d_beta_fix_loss
    d_fix_fix_fix_loss = derivative_ingredients$d_fix_fix_fix_loss
    d_beta_fix_fix_loss = derivative_ingredients$d_beta_fix_fix_loss

    # one way to calculate se based on formula

    # Hessian matrix
    H_a_a = matrix(0, N, N)
    H_g_g = matrix(0, N, N)
    for (i in 1:N) {
      H_a_a[i, i] <- -(sum(d_fix_fix_loss[i, ]) - d_fix_fix_loss[i, i])/(N-1)  # sum of off-diagonal elements in row i
    }
    for (i in 1:N) {
      H_g_g[i, i] <- -(sum(d_fix_fix_loss[-i, i]))/(N-1) # sum of off-diagonal elements in row i
    }
    H_a_g = -d_fix_fix_loss/((N-1))

    Hessian_bar =   cbind(rbind(H_a_a, t(H_a_g)), rbind(H_a_g, H_g_g))  + c(rep(1,N), rep(-1,N)) %*% t( c(rep(1,N), rep(-1,N)) )/N
    Hessain_inverse = solve(Hessian_bar)
    Hessian_a_a = Hessain_inverse[1:(N), 1:(N)]
    Hessian_g_a = Hessain_inverse[(N+1):(N+N), 1:(N)]
    Hessian_a_g = Hessain_inverse[1:(N), (N+1):(N+N)]
    Hessian_g_g = Hessain_inverse[(N+1):(N+N), (N+1):(N+N)]

    # matrix 'Xi'
    the = matrix(0, N, N)
    for (i in 1:N) {
      for (j in 1:N) {
        temp_sum <- 0
        for (s in 1:N) {
          for (t in 1:N) {
            if (t != s) {
              temp_sum <- temp_sum + (-1/(N)) * (Hessian_a_a[i, s] + Hessian_g_a[j, s] + Hessian_a_g[i, t] + Hessian_g_g[j, t]) * d_beta_fix_loss[s, t]
            }
          }
        }
        the[i, j] <- temp_sum
      }
    }

    # W_hat based on Iva ́ n2016
    D_beta_loss = d_beta_loss - d_fix_loss * the
    D_beta_fix_loss = d_beta_fix_loss - d_fix_fix_loss * the
    D_beta_fix_fix_loss = d_beta_fix_fix_loss - d_fix_fix_fix_loss * the

    W_hat =   -(1 / ((N-1)*N)) * (sum(  d_beta_beta_loss  - d_fix_fix_loss * the * the  ) - sum(diag(d_beta_beta_loss  - d_fix_fix_loss * the * the )) ) # -(1 / N) * (  d_beta_beta_big_loss  - d_beta_fix_big_loss * solve(d_fix_2_big_loss) * d_beta_fix_big_loss  )
    B_hat = -(1/N) * compute_B_hat(d_fix_loss, D_beta_fix_loss, 0.5 * D_beta_fix_fix_loss, d_fix_fix_loss, L)
    D_hat = -(1 / (N-1)) * sum(colSums((d_fix_loss * D_beta_fix_loss + 0.5 * D_beta_fix_fix_loss) * (1 - diag(N))) / colSums(d_fix_fix_loss * (1 - diag(N))))
    Omega_hat <- 0
    for (i in 1:N) {
      for (j in setdiff(1:N, i)) {
        for (k in setdiff(1:N, i)) {
          Omega_hat <- Omega_hat + D_beta_loss[i, j] * D_beta_loss[i, k]/(N*(N-1))
        }
      }
    }


    #### Another way to calculate the bias corrected
    # another way to calculate the analytical corrected estimate
    H = phi_XB/(Phi_XB*(1-Phi_XB))
    small_w = H * phi_XB

    # W_hat
    X_into = matrix_to_panel_df(X)
    weight = matrix_to_panel_df(small_w)$X
    weight[which(weight==0)] = 1
    re = get_weighted_projection_fitted_exclude_t_eq_i(X_into$X, weight, X_into$id, X_into$time)
    re_matrix = matrix(re, N-1, N)
    re_matrix = shift_lower_triangle_and_add_zero_diag(re_matrix)
    tilde_X = X - re_matrix
    W_hat_another =   (1 / ((N-1)*N)) * (sum(  small_w  * tilde_X  * tilde_X   ) - sum(diag( small_w  * tilde_X  * tilde_X  )) ) # -(1 / N) * (  d_beta_beta_big_loss  - d_beta_fix_big_loss * solve(d_fix_2_big_loss) * d_beta_fix_big_loss
    D_hat_another = -(0.5 / (N-1)) * sum(colSums((H * dd_F_fix * tilde_X) * (1 - diag(N))) / colSums(small_w * (1 - diag(N))))
    B_hat_another = -(0.5/N) * compute_B_hat_another(2* H * (y - Phi_XB), small_w * tilde_X, H * dd_F_fix * tilde_X, small_w, L)


    # estimate
    estimate_analytical[index_covariate] = est[index_covariate] - solve(W_hat)*B_hat*(1/(N-1)) - solve(W_hat)*D_hat*(1/N)
    estimate_analytical_another[index_covariate] = est[index_covariate] - solve(W_hat_another)*B_hat_another*(1/(N-1)) - solve(W_hat_another)*D_hat_another*(1/N)
    se[index_covariate] = sqrt(solve(W_hat)/(N*(N-1)))
    se_no_MLE[index_covariate] = sqrt(solve(W_hat) * Omega_hat * solve(W_hat)/(N*(N-1)))
    se_another[index_covariate] = sqrt(solve(W_hat_another)/(N*(N-1)))
  }

  # calculate the estimates of fixed effects based on the bias corrected beta
  data_2 <- data_in
  formula <- as.formula( paste("y ~ -1 +", paste(colnames(data_2[,-seq(2,K+1,1)])[-1], collapse = " + "), "+ offset(offset_term)"))
  if( K != 1){
    data_2$offset_term <-  as.vector(X_design[,1:K] %*% matrix(estimate_analytical))
  }else{
    data_2$offset_term <- as.vector(estimate_analytical * X_design[,1:K])
  }
  model_j_2 <-
    speedglm(formula, data = data_2, family=binomial(link = 'probit'))
  fit_j_2 = summary(model_j_2)
  estimate_analytical = c(estimate_analytical, unlist(as.list(fit_j_2$coefficients[, 1])))
  estimate_analytical[K+1] = sum(estimate_analytical[(N+K+1):(N+N+K)]) - sum(estimate_analytical[(K+2):(N+K)])


  data_2 <- data_in
  formula <- as.formula( paste("y ~ -1 +", paste(colnames(data_2[,-seq(2,K+1,1)])[-1], collapse = " + "), "+ offset(offset_term)"))
  if( K != 1){
    data_2$offset_term <- as.vector(X_design[,1:K] %*% matrix(estimate_analytical_another))
  }else{
    data_2$offset_term <- as.vector(estimate_analytical_another * X_design[,1:K])
  }
  model_j_2 <-
    speedglm(formula = formula, data = data_2, family=binomial(link = 'probit'))
  fit_j_2 = summary(model_j_2)
  estimate_analytical_another = c(estimate_analytical_another, unlist(as.list(fit_j_2$coefficients[, 1])))
  estimate_analytical_another[K+1] = sum(estimate_analytical_another[(N+K+1):(N+N+K)]) - sum(estimate_analytical[(K+2):(N+K)])

  eta = estimate_analytical_another %*% t(X_design)
  eta_MLE = cof %*% t(X_design)

  res = list(est = estimate_analytical_another, se = se_another, se_formula = se, se_no_MLE = se_no_MLE,
             est_formula = estimate_analytical,
             eta = eta, eta_MLE = eta_MLE, est_MLE = cof)
  return(res)

}

#' @export
get_APE_bootstrap <- function(y, X, N, data, fit, model = 'probit'){
  cof_estimate = fit$est_bootstrap_all[,1]
  eta = fit$eta
  eta_MLE = matrix(fit$eta_MLE)
  cof_MLE = fit$est_MLE[1]
  APE_estimate = matrix(0,length(cof_estimate), N*(N-1))
  if (is.numeric(X) && length(unique(X)) == 2) {
    X_max = max(unique(X))
    X_min = min(unique(X))
    APE_MLE_estimate = (pnorm(cof_MLE * X_max + eta_MLE - cof_MLE * X ) - pnorm(cof_MLE * X_min + eta_MLE - cof_MLE * X) ) / (X_max - X_min)
    for (i in 1:length(cof_estimate)) {
      APE_estimate[i,] = (pnorm(cof_estimate[i]* X_max + eta[i,] - cof_estimate[i] * X ) - pnorm(cof_estimate[i]* X_min + eta[i,] - cof_estimate[i] * X )) / (X_max - X_min)
    }}else{
      APE_MLE_estimate = cof_MLE * dnorm(eta_MLE)
      for (i in 1:length(cof_estimate)) {
        APE_estimate[i,] = cof_estimate[i]* dnorm( eta[i,] )
      }
    }
  APE = rowMeans(APE_estimate)
  APE_MLE = mean(APE_MLE_estimate)
  APE_mean_bootstrap = APE_MLE - (mean(APE) - APE_MLE)
  APE_median_bootstrap = APE_MLE - (median(APE) - APE_MLE)
  se = sd(APE)
  res = list(APE_MLE = APE_MLE, APE_mean_bootstrap = APE_mean_bootstrap, APE_median_bootstrap = APE_median_bootstrap, se = se)
  return(res)
}

#' @export
get_APE_analytical<- function(y, X, N, data, index, fit, L = 1, model = 'probit'){
  cof_estimate = fit$est[1]
  cof_estimate_MLE = fit$est_MLE[1]
  eta = t(fit$eta)
  eta_MLE = matrix(fit$eta_MLE)
  if (is.numeric(X) && length(unique(X)) == 2) {
    X_max = as.numeric(max(unique(X)))
    X_min = as.numeric(min(unique(X)))
    APE_estimate = (pnorm(cof_estimate* X_max + eta - cof_estimate * X ) - pnorm(cof_estimate* X_min + eta - cof_estimate * X)) / (X_max - X_min)
    d_APE_estimate = (dnorm(cof_estimate* X_max + eta - cof_estimate * X ) - dnorm(cof_estimate* X_min + eta - cof_estimate * X)) / (X_max - X_min)
    dd_APE_estimate = (-(cof_estimate* X_max + eta - cof_estimate * X) * dnorm(cof_estimate* X_max + cof_estimate + eta - cof_estimate * X ) + (cof_estimate* X_min + eta - cof_estimate * X)  * dnorm(cof_estimate* X_min + eta - cof_estimate * X))/ (X_max - X_min)

    APE_estimate_MLE = (pnorm(cof_estimate_MLE * X_max + eta_MLE - cof_estimate_MLE * X ) - pnorm(cof_estimate_MLE* X_min + eta_MLE - cof_estimate_MLE * X)) / (X_max - X_min)
    d_APE_estimate_MLE = (dnorm(cof_estimate_MLE * X_max + eta_MLE - cof_estimate_MLE * X ) - dnorm(cof_estimate_MLE* X_min +eta_MLE - cof_estimate_MLE * X) ) / (X_max - X_min)
    dd_APE_estimate_MLE = (-(cof_estimate_MLE * X_max + eta_MLE - cof_estimate_MLE * X) * dnorm(cof_estimate_MLE* X_max + eta_MLE - cof_estimate_MLE * X ) + (cof_estimate_MLE* X_min + eta_MLE - cof_estimate_MLE * X)  * dnorm(cof_estimate_MLE* X_min + eta_MLE - cof_estimate_MLE * X)) / (X_max - X_min)

  }else{
    APE_estimate = cof_estimate * dnorm(eta)
    d_APE_estimate = cof_estimate * (-eta) * dnorm(eta)
    dd_APE_estimate = cof_estimate*(- dnorm(eta) + eta^2 * dnorm(eta))

    APE_estimate_MLE = cof_estimate_MLE * dnorm(eta_MLE)
    d_APE_estimate_MLE = cof_estimate_MLE * (-eta_MLE) * dnorm(eta_MLE)
    dd_APE_estimate_MLE = cof_estimate_MLE * (- dnorm(eta_MLE) + eta_MLE^2 * dnorm(eta_MLE))
  }
  APE_MLE = mean(APE_estimate_MLE)
  index1 = data[,index[1]]
  index2 = data[,index[2]]
  new_index1 <- match(index1, sort(unique(index1)))
  new_index2 <- match(index2, sort(unique(index2)))

  data[,index[1]] =  new_index1
  data[,index[2]] =  new_index2

  X = vector_to_matrix(X, N, ind1 = data[,index[1]], ind2 = data[,index[2]])
  y = vector_to_matrix(y, N, ind1 = data[,index[1]], ind2 = data[,index[2]])
  APE_estimate = matrix(APE_estimate, N-1,N)
  APE_estimate = shift_lower_triangle_and_add_zero_diag(APE_estimate)
  d_APE_estimate = matrix(d_APE_estimate, N-1,N)
  d_APE_estimate = shift_lower_triangle_and_add_zero_diag(d_APE_estimate)
  dd_APE_estimate = matrix(dd_APE_estimate, N-1,N)
  dd_APE_estimate = shift_lower_triangle_and_add_zero_diag(dd_APE_estimate)
  APE_estimate_MLE = matrix(APE_estimate_MLE, N-1,N)
  APE_estimate_MLE = shift_lower_triangle_and_add_zero_diag(APE_estimate_MLE)
  d_APE_estimate_MLE = matrix(d_APE_estimate_MLE, N-1,N)
  d_APE_estimate_MLE = shift_lower_triangle_and_add_zero_diag(d_APE_estimate_MLE)
  dd_APE_estimate_MLE = matrix(dd_APE_estimate_MLE, N-1,N)
  dd_APE_estimate_MLE = shift_lower_triangle_and_add_zero_diag(dd_APE_estimate_MLE)
  y[row(y) == col(y)] = 0

  cov_sum = matrix(eta_MLE, N-1,N)
  cov_sum = shift_lower_triangle_and_add_zero_diag(cov_sum)
  cov_sum[row(cov_sum) == col(cov_sum)] = 0

  Phi_XB <- pnorm(cov_sum)  # CDF (Φ(Xβ))
  phi_XB <- dnorm(cov_sum)  # PDF (φ(Xβ))

  Phi_XB[row(cov_sum) == col(cov_sum)] = 0
  phi_XB[row(cov_sum) == col(cov_sum)] = 0
  Phi_XB <- pmax(Phi_XB, 1e-9)
  Phi_XB <- pmin(Phi_XB, 1 - 1e-9)

  dd_F_fix = -cov_sum * phi_XB
  ddd_F_fix = cov_sum^2 * phi_XB - phi_XB

  # derivative_ingredients = compute_derivatives(eta = cov_sum, y = y, X = X)
  #
  # # preparation for ingredients
  # d_fix_loss =  derivative_ingredients$d_fix_loss
  # d_beta_loss = derivative_ingredients$d_beta_loss
  # d_fix_fix_loss = derivative_ingredients$d_fix_fix_loss# d_beta_loss = X * (y - exp(cov_sum)/( 1 + exp(cov_sum)))
  # d_beta_beta_loss = derivative_ingredients$d_beta_beta_loss# d_beta_fix_loss = X * ( - exp(cov_sum)/( 1 + exp(cov_sum))^2)
  # d_beta_fix_loss = derivative_ingredients$d_beta_fix_loss
  # d_fix_fix_fix_loss = derivative_ingredients$d_fix_fix_fix_loss
  # d_beta_fix_fix_loss = derivative_ingredients$d_beta_fix_fix_loss
  #
  # # one way to calculate se based on formula
  #
  # # Hessian matrix
  # H_a_a = matrix(0, N, N)
  # H_g_g = matrix(0, N, N)
  # for (i in 1:N) {
  #   H_a_a[i, i] <- -(sum(d_fix_fix_loss[i, ]) - d_fix_fix_loss[i, i])/(N-1)  # sum of off-diagonal elements in row i
  # }
  # for (i in 1:N) {
  #   H_g_g[i, i] <- -(sum(d_fix_fix_loss[-i, i]))/(N-1) # sum of off-diagonal elements in row i
  # }
  # H_a_g = -d_fix_fix_loss/((N-1))
  #
  # Hessian_bar =   cbind(rbind(H_a_a, t(H_a_g)), rbind(H_a_g, H_g_g))  + c(rep(1,N), rep(-1,N)) %*% t( c(rep(1,N), rep(-1,N)) )/N
  # Hessain_inverse = solve(Hessian_bar)
  # Hessian_a_a = Hessain_inverse[1:(N), 1:(N)]
  # Hessian_g_a = Hessain_inverse[(N+1):(N+N), 1:(N)]
  # Hessian_a_g = Hessain_inverse[1:(N), (N+1):(N+N)]
  # Hessian_g_g = Hessain_inverse[(N+1):(N+N), (N+1):(N+N)]
  #
  # # matrix 'Xi'
  # psi = matrix(0, N, N)
  # for (i in 1:N) {
  #   for (j in 1:N) {
  #     temp_sum <- 0
  #     for (s in 1:N) {
  #       for (t in 1:N) {
  #         if (t != s) {
  #           temp_sum <- temp_sum + (-1/(N)) * (Hessian_a_a[i, s] + Hessian_g_a[j, s] + Hessian_a_g[i, t] + Hessian_g_g[j, t]) * d_APE_estimate_MLE[s, t]
  #         }
  #       }
  #     }
  #     psi[i, j] <- temp_sum
  #   }
  # }
  #
  # # matrix 'Xi'
  # the = matrix(0, N, N)
  # for (i in 1:N) {
  #   for (j in 1:N) {
  #     temp_sum <- 0
  #     for (s in 1:N) {
  #       for (t in 1:N) {
  #         if (t != s) {
  #           temp_sum <- temp_sum + (-1/(N)) * (Hessian_a_a[i, s] + Hessian_g_a[j, s] + Hessian_a_g[i, t] + Hessian_g_g[j, t]) * d_beta_fix_loss[s, t]
  #         }
  #       }
  #     }
  #     the[i, j] <- temp_sum
  #   }
  # }
  #
  # # W_hat based on Iva ́ n2016
  # D_beta_loss = d_beta_loss - d_fix_loss * the
  # D_beta_fix_loss = d_beta_fix_loss - d_fix_fix_loss * the
  # D_beta_fix_fix_loss = d_beta_fix_fix_loss - d_fix_fix_fix_loss * the
  #
  # psi*d_fix_loss
  # sum(d_APE_estimate_MLE*tilde_X)/(N*(N-1))
  #
  #
  # as.numeric(sum(d_APE_estimate_MLE*X - psi*d_APE_estimate_MLE)/(N*(N-1))*solve(W_hat_another))* H * (y - Phi_XB)*tilde_X - Psi_matrix * H * (y - Phi_XB)
  #### Another way to calculate the bias corrected
  # another way to calculate the analytical corrected estimate
  H = phi_XB/(Phi_XB*(1-Phi_XB))
  small_w = H * phi_XB

  # B_hat
  X_into = matrix_to_panel_df(-d_APE_estimate_MLE/small_w)
  weight = matrix_to_panel_df(small_w)$X

  Psi = get_weighted_projection_fitted_exclude_t_eq_i(X_into$X, weight, X_into$id, X_into$time)
  Psi_matrix = matrix(Psi, N-1, N)
  Psi_matrix = shift_lower_triangle_and_add_zero_diag(Psi_matrix)
  tilde_Psi = -d_APE_estimate_MLE/small_w - Psi_matrix
  diag(tilde_Psi) = 0
  # D_hat
  D_hat_another = (0.5 / (N-1)) * sum(colSums((dd_APE_estimate_MLE -  Psi_matrix * H * dd_F_fix) * (1 - diag(N))) / colSums(small_w * (1 - diag(N))))

  # B_hat
  B_hat_another = (0.5/N) * compute_B_hat_another(2* H * (y - Phi_XB), small_w * tilde_Psi, dd_APE_estimate_MLE - H * dd_F_fix * Psi_matrix , small_w, L)

  APE_analytical = mean(APE_estimate) - B_hat_another*(1/(N-1)) - D_hat_another*(1/N)

  # residual X
  X_into = matrix_to_panel_df(X)
  weight = matrix_to_panel_df(small_w)$X
  re = get_weighted_projection_fitted_exclude_t_eq_i(X_into$X, weight, X_into$id, X_into$time)
  re_matrix = matrix(re, N-1, N)
  re_matrix = shift_lower_triangle_and_add_zero_diag(re_matrix)
  tilde_X = X - re_matrix
  W_hat_another =   (1 / ((N-1)*N)) * (sum(  small_w  * tilde_X  * tilde_X   ) - sum(diag( small_w  * tilde_X  * tilde_X  )) )
  tau = as.numeric((sum(d_APE_estimate_MLE*tilde_X)/(N*(N-1)))*solve(W_hat_another))* H * (y - Phi_XB)*tilde_X - Psi_matrix * H * (y - Phi_XB)
  APE_residual = APE_estimate_MLE - colMeans(APE_estimate_MLE) # indentical assumption
  diag(APE_residual) = 0
  se = sqrt( compute_expression(APE_residual, tau))/(N*(N-1))
  res = list(APE_MLE = APE_MLE, APE = APE_analytical, se = se)

  return(res)

}

#' @export
get_APE_jackknife<- function(y, X, N, index, data, fit, L = 1, model = 'probit'){
  cof_estimate = fit$est_jackknife_all[,1]
  cof_estimate_MLE = fit$est_MLE[1]
  eta = t(fit$eta)
  eta_MLE = t(fit$eta_MLE)
  if (is.numeric(X) && length(unique(X)) == 2) {
    X_max = as.numeric(max(unique(X)))
    X_min = as.numeric(min(unique(X)))
    APE_estimate = matrix(0,N-1,N*(N-1))
    d_APE_estimate = matrix(0,N-1,N*(N-1))
    dd_APE_estimate = matrix(0,N-1,N*(N-1))
    for (i in 1:N-1) {
      APE_estimate[i,] = (pnorm(cof_estimate[i]* X_max + eta[,i] - cof_estimate[i] * X ) - pnorm(cof_estimate[i]* X_min + eta[,i] - cof_estimate[i] * X)) / (X_max - X_min)
      d_APE_estimate[i,] = (dnorm(cof_estimate[i]* X_max + eta[,i] - cof_estimate[i] * X ) - dnorm(cof_estimate[i]* X_min + eta[,i] - cof_estimate[i] * X)) / (X_max - X_min)
      dd_APE_estimate[i,] = (-(cof_estimate[i]* X_max + eta[,i] - cof_estimate[i] * X) * dnorm(cof_estimate[i]* X_max + eta[,i] - cof_estimate[i] * X) + (cof_estimate[i]* X_min + eta[,i] - cof_estimate[i] * X)  * dnorm(cof_estimate[i]* X_min + eta[,i] - cof_estimate[i] * X))/ (X_max - X_min)
    }
    APE_estimate_MLE = (pnorm(cof_estimate_MLE * X_max + eta_MLE - cof_estimate_MLE * X ) - pnorm(cof_estimate_MLE* X_min + eta_MLE - cof_estimate_MLE * X))/(X_max - X_min)
    d_APE_estimate_MLE = (dnorm(cof_estimate_MLE * X_max + eta_MLE - cof_estimate_MLE * X ) - dnorm(cof_estimate_MLE* X_min +eta_MLE - cof_estimate_MLE * X) ) / (X_max - X_min)
    dd_APE_estimate_MLE = (-(cof_estimate_MLE * X_max + eta_MLE - cof_estimate_MLE * X) * dnorm(cof_estimate_MLE* X_max + eta_MLE - cof_estimate_MLE * X ) + (cof_estimate_MLE* X_min + eta_MLE - cof_estimate_MLE * X)  * dnorm(cof_estimate_MLE* X_min + eta_MLE - cof_estimate_MLE * X)) / (X_max - X_min)

  }else{
    APE_estimate = matrix(0,N-1,N*(N-1))
    d_APE_estimate = matrix(0,N-1,N*(N-1))
    dd_APE_estimate = matrix(0,N-1,N*(N-1))
    for (i in 1:N-1) {
      APE_estimate[i,] = cof_estimate[i] * dnorm(eta[,i])
      d_APE_estimate[i,] = cof_estimate[i] * (-eta[,i]) * dnorm(eta[,i])
      dd_APE_estimate[i,] = cof_estimate[i]*(- dnorm(eta[,i]) + eta[,i]^2 * dnorm(eta[,i]))
    }
    APE_estimate_MLE = cof_estimate_MLE * dnorm(eta_MLE)
    d_APE_estimate_MLE = cof_estimate_MLE * (-eta_MLE) * dnorm(eta_MLE)
    dd_APE_estimate_MLE = cof_estimate_MLE * (- dnorm(eta_MLE) + eta_MLE^2 * dnorm(eta_MLE))
  }
  APE_MLE = mean(APE_estimate_MLE)
  APE_jack = (N-1)*mean(APE_estimate_MLE) - (N-2)* mean( apply(APE_estimate,1,mean) )

  index1 = data[,index[1]]
  index2 = data[,index[2]]
  new_index1 <- match(index1, sort(unique(index1)))
  new_index2 <- match(index2, sort(unique(index2)))

  data[,index[1]] =  new_index1
  data[,index[2]] =  new_index2

  X = vector_to_matrix(X, N, ind1 = data[,index[1]], ind2 = data[,index[2]])
  y = vector_to_matrix(y, N, ind1 = data[,index[1]], ind2 = data[,index[2]])
  # APE_estimate = matrix(APE_estimate, N-1,N)
  # APE_estimate = shift_lower_triangle_and_add_zero_diag(APE_estimate)
  # d_APE_estimate = matrix(d_APE_estimate, N-1,N)
  # d_APE_estimate = shift_lower_triangle_and_add_zero_diag(d_APE_estimate)
  # dd_APE_estimate = matrix(dd_APE_estimate, N-1,N)
  # dd_APE_estimate = shift_lower_triangle_and_add_zero_diag(dd_APE_estimate)
  APE_estimate_MLE = matrix(APE_estimate_MLE, N-1,N)
  APE_estimate_MLE = shift_lower_triangle_and_add_zero_diag(APE_estimate_MLE)
  d_APE_estimate_MLE = matrix(d_APE_estimate_MLE, N-1,N)
  d_APE_estimate_MLE = shift_lower_triangle_and_add_zero_diag(d_APE_estimate_MLE)
  dd_APE_estimate_MLE = matrix(dd_APE_estimate_MLE, N-1,N)
  dd_APE_estimate_MLE = shift_lower_triangle_and_add_zero_diag(dd_APE_estimate_MLE)
  y[row(y) == col(y)] = 0

  cov_sum= matrix(eta_MLE, N-1,N)
  cov_sum = shift_lower_triangle_and_add_zero_diag(cov_sum)
  cov_sum[row(cov_sum) == col(cov_sum)] = 0

  Phi_XB <- pnorm(cov_sum)  # CDF (Φ(Xβ))
  phi_XB <- dnorm(cov_sum)  # PDF (φ(Xβ))

  Phi_XB[row(cov_sum) == col(cov_sum)] = 0
  phi_XB[row(cov_sum) == col(cov_sum)] = 0
  Phi_XB <- pmax(Phi_XB, 1e-9)
  Phi_XB <- pmin(Phi_XB, 1 - 1e-9)
  dd_F_fix = -cov_sum * phi_XB
  ddd_F_fix = cov_sum^2 * phi_XB - phi_XB

  #### Another way to calculate the bias corrected
  # another way to calculate the analytical corrected estimate
  H = phi_XB/(Phi_XB*(1-Phi_XB))
  small_w = H * phi_XB

  # B_hat
  X_into = matrix_to_panel_df(-d_APE_estimate_MLE/small_w)
  weight = matrix_to_panel_df(small_w)$X

  Psi = get_weighted_projection_fitted_exclude_t_eq_i(X_into$X, weight, X_into$id, X_into$time)
  Psi_matrix = matrix(Psi, N-1, N)
  Psi_matrix = shift_lower_triangle_and_add_zero_diag(Psi_matrix)
  tilde_Psi = -d_APE_estimate_MLE/small_w - Psi_matrix
  diag(tilde_Psi) = 0
  # D_hat
  D_hat_another = (0.5 / (N-1)) * sum(colSums((dd_APE_estimate_MLE -  Psi_matrix * H * dd_F_fix) * (1 - diag(N))) / colSums(small_w * (1 - diag(N))))

  # B_hat
  B_hat_another = (0.5/N) * compute_B_hat_another(2* H * (y - Phi_XB), small_w * tilde_Psi, dd_APE_estimate_MLE - H * dd_F_fix * Psi_matrix , small_w, L)


  # residual X
  X_into = matrix_to_panel_df(X)
  weight = matrix_to_panel_df(small_w)$X
  re = get_weighted_projection_fitted_exclude_t_eq_i(X_into$X, weight, X_into$id, X_into$time)
  re_matrix = matrix(re, N-1, N)
  re_matrix = shift_lower_triangle_and_add_zero_diag(re_matrix)
  tilde_X = X - re_matrix
  W_hat_another =   (1 / ((N-1)*N)) * (sum(  small_w  * tilde_X  * tilde_X   ) - sum(diag( small_w  * tilde_X  * tilde_X  )) )
  tau = as.numeric((sum(d_APE_estimate_MLE*tilde_X)/(N*(N-1)))*solve(W_hat_another))* H * (y - Phi_XB)*tilde_X - Psi_matrix * H * (y - Phi_XB)
  APE_residual = APE_estimate_MLE - colMeans(APE_estimate_MLE) # indentical assumption
  diag(APE_residual) = 0
  se = sqrt( compute_expression(APE_residual, tau))/(N*(N-1))
  res = list(APE_MLE = APE_MLE, APE = APE_jack, se = se)

  return(res)

}
