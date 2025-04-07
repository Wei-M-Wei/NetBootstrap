#' Parametric bootstrap inference in network formation models
#'@title Parametric bootstrap inference in network formation models
#'
#' @description R package 'Netbootstrap' is dedicated to do parametric bootstrap inference in network formation models.
#'
#' @param y outcome variabe
#' @param X dependent variables
#' @param N size of the individuals
#' @param bootstrap_time number of the bootstrap times
#' @param data data which you want to apply
#' @param index name of the import side and export side
#' @param link probit or logit model
#'
#' @returns A list of fitted results is returned.
#' Within this outputted list, the following elements can be found:
#'     \item{est}{estimated coefficients from the bootstrap procedure.}
#'     \item{sd}{stand deviation form the bootstrap procedure.}
#'     \item{est_MLE}{estimated coefficients from MLE.}
#'
#' @import speedglm
#' @export
#'
#'
network_bootstrap = function(y, X, N, bootstrap_time, data, index, link = 'probit'){

  data = data.frame(y = y, X = X,  data[,index[1]], data[,index[2]])
  # order the data
  colnames(data) = c('y', colnames(X), index)
  data = data[order(data[,index[1]], data[,index[2]]),]
  K = dim(X)[2]

  y = data$y
  X = data[,colnames(X)]

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

  X_design = cbind(X, fix)
  data_input <- data.frame(y = y, X = X_design)

  # MLE
  model <-
    speedglm(y ~ . - 1, data = data_input, family = binomial(link = link))
  fit = summary(model)
  sd_MLE = vcov(model)
  cof = fit$coefficients$Estimate
  cof[K+1] = sum(cof[(N+K+1):(N+N+K)]) - sum(cof[(K+2):(N+K)])

  # bootstrap procedure
  X_design = apply(X_design, 2, as.numeric)
  cof_boost = c()
  for (k in seq(bootstrap_time)) {
    epsi_it = rnorm(n, 0, 1)
    Y = as.numeric(X_design %*%cof > epsi_it)
    data_boostrap <- data.frame(y = Y, X = X_design)
    model <-
      speedglm(y ~ . - 1 , data = data_boostrap, family = binomial(link = link))
    fit = summary(model)
    cof_boost = rbind(cof_boost, fit$coefficients$Estimate)
  }

  # results
  cof_use = cof_boost
  cof_boost[,K+1] = sum(cof_boost[,(N+K+1):(N+N+K)]) - sum(cof_boost[,(K+2):(N+K)])
  boostrap_sd = apply(cof_boost[,1:K],2,sd)
  cof_boost_mean = apply(cof_boost, 2, mean)
  est_correct = cof - (cof_boost_mean - cof)
  est_correct[K+1] = sum(est_correct[(N+K+1):(N+N+K)]) - sum(est_correct[(K+2):(N+K)])
  bias_se = (est_correct - cof)[1:K]/ boostrap_sd
  return(list(est = est_correct[1:K], sd = boostrap_sd[1:K], est_MLE = cof[1:K]))
}
