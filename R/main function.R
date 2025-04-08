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

  data = data.frame(y = y, X = X,  data[,index[1]], data[,index[2]])

  # order the data
  if(is.null(colnames(X)) == 1){
    colnames(data)[(dim(data)[2]-1):dim(data)[2]] = index
    data = data[order(data[,index[1]], data[,index[2]]),]
    p = dim(X)[2]
    y = data$y
    if (p==1){
      X = data[,2]
    }else{
      X = data[,(2:p+1)]
    }
  }else{
    colnames(data) = c('y', colnames(X), index)
    data = data[order(data[,index[1]], data[,index[2]]),]
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
  data <- data.frame(y = y, X = X_design)

  if (sum(y) <= N*N*0.01){
    return(NULL)
  }

  # MLE
  model <-
    speedglm(y ~ . - 1, data = data, family = binomial(link = link))
  fit = summary(model)
  Hessian_inv = vcov(model)
  cof = unlist(as.list(fit$coefficients[, 1]))
  log_likelihood_estimate <- logLik(model)
  cof[K+1] = sum(cof[(N+K+1):(N+N+K)]) - sum(cof[(K+2):(N+K)])

  # constrained MLE
  data_2 = data
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
      data_2$offset_term <- beta_boot * X_design[,1]
      model_B_NULL <-
        speedglm(formula = formula, data = data_2, family=binomial(link = link))
      fit_B_NULL = summary(model_B_NULL)
      cof_B_NULL = rbind(cof_B_NULL, c(beta_NULL, NA, unlist(as.list(fit_B_NULL$coefficients[, 1]))))
      log_likelihood_estimate_B_NULL <- rbind(log_likelihood_estimate_B_NULL, logLik(model_B_NULL))
    }
  }

  # get the final results
  cof_B[,K+1] = sum(cof_B[,(N+K+1):(N+N+K)]) - sum(cof_B[,(K+2):(N+K)])
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
    res = list(cof_MLE = cof, cof_mean = est_correct_mean, cof_median = est_correct_median, sd = boostrap_sd,
               cof_bootstrap_all = cof_B, cof_MLE_NULL = cof_NULL, cof_bootstrap_NULL = cof_B_NULL,
               log_likelihood_MLE = log_likelihood_estimate, log_likelihood_Bootstrap = log_likelihood_estimate_B,
               log_likelihood_MLE_NULL = log_likelihood_estimate_NULL, log_likelihood_Bootstrap_NULL = log_likelihood_estimate_B_NULL,
               Hessian_MLE = Hessian_inv, Hessian_MLE_NULL = Hessian_inv_NULL, X_origin = as.matrix(X_design)
    )
  }
  else{
    res = list(cof_MLE = cof, cof_mean = est_correct_mean, cof_median = est_correct_median, sd = boostrap_sd, cof_bootstrap_all = cof_B,
               log_likelihood_MLE = log_likelihood_estimate, log_likelihood_Bootstrap = log_likelihood_estimate_B,
               Hessian_MLE = Hessian_inv, X_origin = as.matrix(X_design)
    )
  }
  return(res)
}

#' @export
#'
APE_se = function(X, y, cov_APE, cov_APE_minus, cov_sum, APE_MLE, model = 'probit'){
  cov_sum[row(cov_sum) == col(cov_sum)] = 0
  cov_APE[row(cov_APE) == col(cov_APE)] = 0
  cov_APE_minus[row(cov_APE_minus) == col(cov_APE_minus)] = 0

  XB_pi <- cov_sum
  Phi_XB <- pnorm(cov_sum)  # CDF (Φ(Xβ))
  phi_XB <- dnorm(cov_sum)  # PDF (φ(Xβ))

  Phi_XB <- pmax(Phi_XB, 1e-9)
  Phi_XB <- pmin(Phi_XB, 1 - 1e-9)

  # preparation for ingredients
  # d_fix_loss = y - exp(cov_sum)/( 1 + exp(cov_sum))
  d_fix_loss = - (phi_XB * (y - Phi_XB) / (Phi_XB * (1 - Phi_XB)))
  #d_fix_2_loss = - exp(cov_sum)/( 1 + exp(cov_sum))^2
  d_fix_2_loss = (-cov_sum*phi_XB  *Phi_XB * ( 1 -Phi_XB ) - phi_XB*(1-2*Phi_XB)*phi_XB)*(y-Phi_XB)/((Phi_XB*(1-Phi_XB)))^2 - phi_XB^2  / ((Phi_XB*(1-Phi_XB)))

  # d_beta_loss = X * (y - exp(cov_sum)/( 1 + exp(cov_sum)))
  d_beta_loss = -X * (phi_XB * (y - Phi_XB) / (Phi_XB * (1 - Phi_XB)))

  # d_beta_beta_loss = X * (- exp(cov_sum)/( 1 + exp(cov_sum))^2) * X
  d_beta_beta_loss = (-cov_sum*phi_XB * X *Phi_XB * ( 1 -Phi_XB ) - phi_XB*(1-2*Phi_XB)*X*phi_XB)*(y-Phi_XB)/((Phi_XB*(1-Phi_XB)))^2 - X * phi_XB^2 * X / ((Phi_XB*(1-Phi_XB)))

  # d_beta_fix_loss = X * ( - exp(cov_sum)/( 1 + exp(cov_sum))^2)
  d_beta_fix_loss = (-cov_sum*phi_XB * X *Phi_XB * ( 1 -Phi_XB ) - phi_XB*(1-2*Phi_XB)*phi_XB)*(y-Phi_XB)/((Phi_XB*(1-Phi_XB)))^2 - X * phi_XB^2  / ((Phi_XB*(1-Phi_XB)))


  d_fix_APE =  (1/sqrt(2*pi) * exp( -cov_APE^2 / 2 ) - 1/sqrt(2*pi) * exp( -cov_APE_minus^2 / 2 ))/2
  d_beta_APE = (1/sqrt(2*pi) * exp( -cov_APE^2 / 2 ) + 1/sqrt(2*pi) * exp( -cov_APE_minus^2 / 2 ))/2


  # Hessian matrix
  H_a_a = -diag(rowSums(d_fix_2_loss))/sqrt(N*(N-1))
  H_g_g = -diag(colSums(d_fix_2_loss))/sqrt(N*(N-1))
  H_a_g = -d_fix_2_loss/sqrt(N*(N-1))

  Hessian_bar =   cbind(rbind(H_a_a, t(H_a_g)), rbind(H_a_g, H_g_g))  + c(rep(1,N), rep(-1,N)) %*% t( c(rep(1,N), rep(-1,N)) )/sqrt(N*(N-1))

  Hessain_inverse = solve(Hessian_bar)
  Hessian_a_a = Hessain_inverse[1:(N), 1:(N)]
  Hessian_g_a = Hessain_inverse[(N+1):(N+N), 1:(N)]
  Hessian_a_g = Hessain_inverse[1:(N), (N+1):(N+N)]
  Hessian_g_g = Hessain_inverse[(N+1):(N+N), (N+1):(N+N)]

  # the matrix
  the = matrix(0,N,N)
  for (i in 1:N) {
    for (j in 1:N) {
      the[i,j] = -(sum(Hessian_a_a[i,]*colSums(d_beta_fix_loss)) +
                     sum(Hessian_g_a[j,]*colSums(d_beta_fix_loss)) +
                     sum(Hessian_a_g[i,]*rowSums(d_beta_fix_loss)) +
                     sum(Hessian_g_g[j,]*rowSums(d_beta_fix_loss)))/sqrt(N*(N-1))
    }
  }


  # W_infity
  W_bar =  -(1 / sqrt((N)*N)) * sum(  d_beta_beta_loss  - d_fix_2_loss * the * the  )


  # phi matrix
  phi = matrix(0, N, N)
  for (i in 1:(N)) {
    for (j in 1:N) {
      phi[i,j] = -(sum(Hessian_a_a[i,]*colSums(d_fix_APE)) +
                     sum(Hessian_g_a[j,]*colSums(d_fix_APE)) +
                     sum(Hessian_a_g[i,]*rowSums(d_fix_APE)) +
                     sum(Hessian_g_g[j,]*rowSums(d_fix_APE)))/sqrt(N*(N-1))
    }
  }


  D_beta_delta = (1/((N-1)*N)) * sum( d_beta_APE - the * d_fix_APE)

  D_beta_loss = d_beta_loss - d_fix_loss*the

  D_beta_dix_loss = d_beta_fix_loss - d_fix_2_loss*the

  tau = as.numeric(D_beta_delta * solve(W_bar)) * D_beta_loss - phi * d_fix_loss

  APE_residual = APE_MLE - colMeans(APE_MLE)


  part_1 = sum(APE_residual%*%APE_residual)
  part_2 = compute_sum(APE_residual)
  part_3 = sum(tau*tau)


  return( sqrt((part_1 + part_2 + part_3)/(N^2*(N-1)^2)) )
}

#' @export
compute_sum <- function(X) {
  N <- nrow(X)
  col_sums <- colSums(X)  # Compute sum over each column

  term1 <- X * matrix(rep(col_sums, each = N), nrow = N)  # X_{it} * sum(X_{jt})
  term2 <- sum(X * X)  # Subtract diagonal terms (i = j)

  return(sum(term1) - term2)
}

#' @export
shift_lower_triangle_and_add_zero_diag <- function(A) {
  N_minus_1 <- nrow(A)
  N <- ncol(A)

  if (N_minus_1 != N - 1) {
    stop("Input matrix must be of size (N-1) x N.")
  }

  # Initialize an N x N zero matrix
  A_new <- matrix(0, nrow = N, ncol = N)

  # Copy upper triangle (including diagonal) from A
  for (i in 1:N_minus_1) {
    for (j in i:N) {
      if (i != j) {  # Exclude diagonal to enforce 0s
        A_new[i, j] <- A[i, j]
      }
    }
  }

  # Shift lower-left triangle (strictly below diagonal) down by 1 row
  for (i in 2:N_minus_1) {
    for (j in 1:(i - 1)) {
      A_new[i, j] <- A[i - 1, j]
    }
  }

  # Final row: only shifted lower-left values
  for (j in 1:(N - 1)) {
    A_new[N, j] <- A[N - 1, j]
  }

  # Diagonal explicitly set to 0 (optional, for safety)
  diag(A_new) <- 0

  return(A_new)
}

