#' @export
get_APE_bootstrap <- function(y, X, N, data, fit, model = 'probit'){
  K = dim(X)[2]
  APE_MLE = rep(0 ,K)
  APE = matrix(0, length(fit$est_bootstrap_all_critical[,1]), K)
  APE_mean_bootstrap = rep(0 ,K)
  APE_median_bootstrap = rep(0 ,K)
  se = rep(0, K)
  X_to = X
  for (k in 1:K) {
    X = X_to[,k]
    cof_estimate = fit$est_bootstrap_all_critical[,k]
    eta = fit$eta_critical
    eta_MLE = matrix(fit$eta_MLE)
    cof_MLE = fit$est_MLE[k]
    APE_estimate = matrix(0,length(cof_estimate), N*(N-1))
    if (is.numeric(X) && length(unique(X)) == 2) {
      X_max = max(unique(X))
      X_min = min(unique(X))
      if (model == 'probit'){
        APE_MLE_estimate = (pnorm(cof_MLE * X_max + eta_MLE - matrix(cof_MLE * X) ) - pnorm(cof_MLE * X_min + eta_MLE - matrix(cof_MLE * X) ) )
        for (i in 1:length(cof_estimate)) {
          APE_estimate[i,] = (pnorm(cof_estimate[i]* X_max + eta[i,] - cof_estimate[i] * X ) - pnorm(cof_estimate[i]* X_min + eta[i,] - cof_estimate[i] * X ))
        }
      }else{
        APE_MLE_estimate = (plogis(cof_MLE * X_max + eta_MLE - matrix(cof_MLE * X) ) - plogis(cof_MLE * X_min + eta_MLE - matrix(cof_MLE * X) ) )
        for (i in 1:length(cof_estimate)) {
          APE_estimate[i,] = (plogis(cof_estimate[i]* X_max + eta[i,] - cof_estimate[i] * X ) - plogis(cof_estimate[i]* X_min + eta[i,] - cof_estimate[i] * X ))
        }
      }

    }else{
      if (model == 'probit'){
        APE_MLE_estimate = cof_MLE * dnorm(eta_MLE)
        for (i in 1:length(cof_estimate)) {
          APE_estimate[i,] = cof_estimate[i]* dnorm( eta[i,] )
        }
      }else{
        APE_MLE_estimate = cof_MLE * dlogis(eta_MLE)
        for (i in 1:length(cof_estimate)) {
          APE_estimate[i,] = cof_estimate[i]* dlogis( eta[i,] )
        }
      }
    }
    APE[,k] = rowMeans(APE_estimate)
    APE_MLE[k] = mean(APE_MLE_estimate)
    APE_mean_bootstrap[k] = APE_MLE[k] - (mean(APE[,k]) - APE_MLE[k])
    APE_median_bootstrap[k] = APE_MLE[k] - (median(APE[,k]) - APE_MLE[k])
    se[k] = sd(APE[,k])
  }
  res = list(APE_MLE = APE_MLE, APE_mean_bootstrap = APE_mean_bootstrap, APE_median_bootstrap = APE_median_bootstrap, APE_bootstrap_all = APE, se = se)
  return(res)
}

#' @export
get_APE_bootstrap_iter <- function(y, X, N, data, index, fit, model = 'probit'){
  original = get_APE_MLE(y = y, X = X, N = N, data = data, index  = index, fit = fit, model = model)
  K = dim(X)[2]
  APE_MLE = rep(0 ,K)
  APE = matrix(0, length(fit$est_bootstrap_all[,1]), K)
  APE_mean_bootstrap = rep(0 ,K)
  APE_median_bootstrap = rep(0 ,K)
  se = rep(0, K)
  X_to = X
  for (k in 1:K) {
    X = X_to[,k]
    cof_estimate = fit$est_bootstrap_all[,k]
    eta = fit$eta
    eta_MLE = matrix(fit$eta_critical_median)
    cof_MLE = fit$est_critical[k]

    APE_estimate = matrix(0,length(cof_estimate), N*(N-1))
    if (is.numeric(X) && length(unique(X)) == 2) {
      X_max = max(unique(X))
      X_min = min(unique(X))
      if (model == 'probit'){
        APE_MLE_estimate = (pnorm(cof_MLE * X_max + eta_MLE - matrix(cof_MLE * X) ) - pnorm(cof_MLE * X_min + eta_MLE - matrix(cof_MLE * X) ) )
        for (i in 1:length(cof_estimate)) {
          APE_estimate[i,] = (pnorm(cof_estimate[i]* X_max + eta[i,] - cof_estimate[i] * X ) - pnorm(cof_estimate[i]* X_min + eta[i,] - cof_estimate[i] * X ))
        }
      }else{
        APE_MLE_estimate = (plogis(cof_MLE * X_max + eta_MLE - matrix(cof_MLE * X) ) - plogis(cof_MLE * X_min + eta_MLE - matrix(cof_MLE * X) ) )
        for (i in 1:length(cof_estimate)) {
          APE_estimate[i,] = (plogis(cof_estimate[i]* X_max + eta[i,] - cof_estimate[i] * X ) - plogis(cof_estimate[i]* X_min + eta[i,] - cof_estimate[i] * X ))
        }
      }

    }else{
      if (model == 'probit'){
        APE_MLE_estimate = cof_MLE * dnorm(eta_MLE)
        for (i in 1:length(cof_estimate)) {
          APE_estimate[i,] = cof_estimate[i]* dnorm( eta[i,] )
        }
      }else{
        APE_MLE_estimate = cof_MLE * dlogis(eta_MLE)
        for (i in 1:length(cof_estimate)) {
          APE_estimate[i,] = cof_estimate[i]* dlogis( eta[i,] )
        }
      }
    }

    APE[,k] = rowMeans(APE_estimate)
    APE_MLE[k] = mean(APE_MLE_estimate)
    APE_mean_bootstrap[k] = original$APE[k] - (mean(APE[,k]) - APE_MLE[k])
    APE_median_bootstrap[k] = original$APE[k] - (median(APE[,k]) - APE_MLE[k])
    se[k] = sd(APE[,k])
  }
  res = list(APE_MLE = APE_MLE, APE_mean_bootstrap = APE_mean_bootstrap, APE_median_bootstrap = APE_median_bootstrap, APE_bootstrap_all = APE, se = se)
  return(res)
}


#' @export
get_APE_analytical<- function(y, X, N, data, index, fit, L = 1, model = 'probit'){
  K = dim(X)[2]
  APE_MLE = rep(0, K)
  APE_analytical = rep(0, K)
  tilde_X_list <- vector("list", K)
  tilde_Psi_list <- vector("list", K)
  APE_list <- vector("list", K)
  d_APE_list <- vector("list", K)
  APE_MLE_list <- vector("list", K)
  d_APE_MLE_list <- vector("list", K)
  Psi_list = vector("list", K)
  re_matrix_array <- array(0, dim = c(N, N, K))
  W_hat <- matrix(0, nrow = K, ncol = K)
  B_hat = matrix(0, nrow = K, ncol = 1)
  D_hat = matrix(0, nrow = K, ncol = 1)
  X_to = X
  y_to = y

  for (k in 1:K) {
    X = X_to[,k]
    cof_estimate = fit$est[k]
    cof_estimate_MLE = fit$est_MLE[k]
    eta = t(fit$eta)
    eta_MLE = t(fit$eta_MLE)
    if (is.numeric(X) && length(unique(X)) == 2) {
      X_max = as.numeric(max(unique(X)))
      X_min = as.numeric(min(unique(X)))
      #test1 = dnorm(cof_estimate* X_max + eta - cof_estimate * X )

      if (model == 'probit'){
        APE_estimate = (pnorm(cof_estimate* X_max + eta - cof_estimate * X ) - pnorm(cof_estimate* X_min + eta - cof_estimate * X))
        d_APE_estimate = (dnorm(cof_estimate* X_max + eta - cof_estimate * X ) - dnorm(cof_estimate* X_min + eta - cof_estimate * X))
        dd_APE_estimate = (-(cof_estimate* X_max + eta - cof_estimate * X) * dnorm(cof_estimate* X_max + eta - cof_estimate * X ) + (cof_estimate* X_min + eta - cof_estimate * X)  * dnorm(cof_estimate* X_min + eta - cof_estimate * X))
      }else{
        APE_estimate = (plogis(cof_estimate* X_max + eta - cof_estimate * X ) - plogis(cof_estimate* X_min + eta - cof_estimate * X))
        d_APE_estimate = (dlogis(cof_estimate* X_max + eta - cof_estimate * X ) - dlogis(cof_estimate* X_min + eta - cof_estimate * X))
        dd_APE_estimate = (dlogis(cof_estimate* X_max + eta - cof_estimate * X ) * ( 1 - 2 * plogis(cof_estimate* X_max + eta - cof_estimate * X ) )
                           - dlogis(cof_estimate* X_min + eta - cof_estimate * X) * ( 1 - 2 * plogis(cof_estimate* X_min + eta - cof_estimate * X)))
      }
    }else{
      if (model == 'probit'){
        APE_estimate = cof_estimate * dnorm(eta)
        d_APE_estimate = cof_estimate * (-eta) * dnorm(eta)
        dd_APE_estimate = cof_estimate*(- dnorm(eta) + eta^2 * dnorm(eta))
      }else{
        APE_estimate = cof_estimate * dlogis(eta)
        d_APE_estimate = cof_estimate * (1 - 2 * plogis(eta)) * dlogis(eta)
        dd_APE_estimate = cof_estimate* (- 2 * dlogis(eta)) * dlogis(eta) + cof_estimate * (1 - 2 * plogis(eta)) * dlogis(eta) * (1-2*plogis(eta))
      }

    }


    index1 = data[,index[1]]
    index2 = data[,index[2]]
    new_index1 <- match(index1, sort(unique(index1)))
    new_index2 <- match(index2, sort(unique(index2)))

    data[,index[1]] =  new_index1
    data[,index[2]] =  new_index2

    X = vector_to_matrix(X, N, ind1 = data[,index[1]], ind2 = data[,index[2]])
    y = vector_to_matrix(y_to, N, ind1 = data[,index[1]], ind2 = data[,index[2]])
    APE_estimate = vector_to_matrix( APE_estimate, N, ind1 = data[,index[1]], ind2 = data[,index[2]])
    d_APE_estimate = vector_to_matrix( d_APE_estimate, N, ind1 = data[,index[1]], ind2 = data[,index[2]])
    dd_APE_estimate = vector_to_matrix( dd_APE_estimate, N, ind1 = data[,index[1]], ind2 = data[,index[2]])

    APE_list[[k]] =  APE_estimate
    d_APE_list[[k]] =  d_APE_estimate

    # use corrected eta to calculate the martices in the formula

    cov_sum = vector_to_matrix( t(eta), N, ind1 = data[,index[1]], ind2 = data[,index[2]])
    if (model == 'probit'){
      Phi_XB <- pnorm(cov_sum)  # CDF (Φ(Xβ))
      phi_XB <- dnorm(cov_sum)  # PDF (φ(Xβ))
      Phi_XB <- pmax(Phi_XB, 1e-5)
      Phi_XB <- pmin(Phi_XB, 1 - 1e-5)
      phi_XB <- pmax(phi_XB, 1e-5)
      phi_XB <- pmin(phi_XB, 1 - 1e-5)
      dd_F_fix = -cov_sum * phi_XB
      H = phi_XB/(Phi_XB*(1-Phi_XB))
      small_w = H * phi_XB
    }else{
      Phi_XB <- plogis(cov_sum)  # CDF (Φ(Xβ))
      phi_XB <- dlogis(cov_sum)  # PDF (φ(Xβ))
      dd_F_fix = phi_XB * ( 1 - 2 * Phi_XB)
      Phi_XB[row(Phi_XB) == col(Phi_XB)] = 0
      phi_XB[row(phi_XB) == col(phi_XB)] = 0
      Phi_XB <- pmax(Phi_XB, 1e-5)
      Phi_XB <- pmin(Phi_XB, 1 - 1e-5)
      phi_XB <- pmax(phi_XB, 1e-5)
      phi_XB <- pmin(phi_XB, 1 - 1e-5)
      H = matrix(1, N, N)
      diag(H) = 0
      small_w = phi_XB
    }


    # B_hat
    X_into = matrix_to_panel_df(-d_APE_estimate/small_w)
    weight = matrix_to_panel_df(small_w)$X
    Psi = get_weighted_fe_projection(X_into$X, weight, X_into$id, X_into$time)
    Psi_list[[k]]= vector_to_matrix(Psi, N, ind1 = X_into$id, ind2 = X_into$time)
    tilde_Psi_list[[k]] = -d_APE_estimate/small_w - Psi_list[[k]]
    diag(tilde_Psi_list[[k]] ) = 0
    # D_hat
    D_hat = (0.5 / (N-1)) * sum(colSums((dd_APE_estimate +  Psi_list[[k]] * H * dd_F_fix) * (1 - diag(N))) / colSums(small_w * (1 - diag(N))))

    # B_hat
    B_hat = (0.5/N) * compute_B_hat(2* H * (y - Phi_XB), -small_w * tilde_Psi_list[[k]], dd_APE_estimate + H * dd_F_fix * Psi_list[[k]] , small_w, L)


    APE_analytical[k] = sum(APE_estimate)/(N*(N-1)) - B_hat*(1/(N-1)) - D_hat*(1/N)

    # residual X
    X_into = matrix_to_panel_df(X)
    weight = matrix_to_panel_df(small_w)$X
    re = get_weighted_fe_projection(X_into$X, weight, X_into$id, X_into$time)
    re_matrix = vector_to_matrix(re, N, ind1 = X_into$id, ind2 = X_into$time)
    tilde_X_list[[k]] = X - re_matrix
    re_matrix_array[,,k] = re_matrix
  }


  tilde_X_array <- array(unlist(tilde_X_list), dim = c(N, N, K))
  Psi_array <- array(unlist(Psi_list), dim = c(N, N, K))
  tilde_Psi_array <- array(unlist(tilde_Psi_list), dim = c(N, N, K))
  APE_array <- array(unlist(APE_list), dim = c(N, N, K))
  d_APE_array <- array(unlist(d_APE_list), dim = c(N, N, K))

  J = matrix(0, K,K)
  for (k in 1:K) {
    X = X_to[,k]
    if (is.numeric(X) && length(unique(X)) == 2){
      X_max = as.numeric(max(unique(X)))
      X_min = as.numeric(min(unique(X)))
      if (model == 'probit'){
        f1_own = dnorm(cof_estimate * X_max + eta - cof_estimate * X )
      }else{
        f1_own = dlogis(cof_estimate* X_max + eta - cof_estimate * X )
      }
      for (j in 1:K){
        diag(tilde_X_array[,,j]) = 0
        diag(d_APE_array[,,k]) = 0
        J[j, k] <- - sum(re_matrix_array[,,j] * d_APE_array[,,k]) / (N*N-1)
      }
      J[k, k] <- sum(f1_own) / (N*N-1) + J[k, k]
      J[- k, k] <- colSums(X_to[, - k, drop = FALSE] * matrix_to_panel_df(as.matrix(d_APE_array[,,k]))$X) / (N*N-1) + J[- k, k]
    }else{
      for (j in 1:K){
        diag(tilde_X_array[,,j]) = 0
        diag(d_APE_array[,,k]) = 0
        J[j, k] <- sum(tilde_X_array[,,j] * d_APE_array[,,k]) / (N*N-1)
      }
      J[k, k] <- sum(phi_XB) / (N*N-1) + J[k, k]
    }
  }

  # Step 4: Initialize sum matrix for outer products
  W_hat <- matrix(0, nrow = K, ncol = K)
  D_beta_APE <- 0
  tau = array(0, dim = c(N, N, K))
  # Step 5: For each position (i, j), extract vector and compute outer product
  for (i in 1:N) {
    for (j in 1:N) {
      if(j!=i){
        v_ij <- tilde_X_array[i, j, ]             # Vector of length p
        W_hat <- W_hat + (1 / ((N-1)*N)) * small_w[i,j] * v_ij %*% t(v_ij)  # Outer product and sum
        D_beta_APE =  D_beta_APE +  (1 / ((N-1)*N)) * d_APE_array[i,j, ] * v_ij
      }
    }
  }

  for (i in 1:N) {
    for (j in 1:N) {
      if(j!=i){
        tau[i,j,] = as.matrix(t(J) %*% solve(W_hat)) %*% as.matrix(H[i,j] * (y[i,j] - Phi_XB[i,j]) * tilde_X_array[i,j,]) - Psi_array[i,j,] * H[i,j] * (y[i,j] - Phi_XB[i,j])
      }
    }
  }


  for (k in 1:K) {
    diag(tau[,,k]) <- 0
  }
  tau_sum = 0
  for (i in 1:N) {
    for (t in 1:N) {
      v <- tau[i, t, ]            # tau_{it} as a vector
      tau_sum <- tau_sum + tcrossprod(v)  # v %*% t(v)
    }
  }

  se =  sqrt( diag( tau_sum ) )/(N*(N-1))
  res = list( APE = APE_analytical, se = se)

  return(res)

}

#' @export
get_APE_jackknife <- function(y, X, N, index, data, fit, L = 1, model = 'probit'){
  K = dim(X)[2]
  APE_jack = rep(0, K)
  tilde_X_list <- vector("list", K)
  tilde_Psi_list <- vector("list", K)
  APE_list <- vector("list", K)
  d_APE_list <- vector("list", K)
  Psi_list = vector("list", K)
  re_matrix_array = array(0, dim = c(N,N,K))
  W_hat <- matrix(0, nrow = K, ncol = K)
  X_start = X
  y_to = y
  for (k in 1:K) {
    X = X_start[,k]
    cof_estimate = fit$est_jackknife_all[,k]
    eta = t(fit$eta)
    if (is.numeric(X) && length(unique(X)) == 2) {
      X_max = as.numeric(max(unique(X)))
      X_min = as.numeric(min(unique(X)))
      APE_estimate = matrix(0,N-1,N*(N-1))
      d_APE_estimate = matrix(0,N-1,N*(N-1))
      dd_APE_estimate = matrix(0,N-1,N*(N-1))
      if (model == 'probit'){
        for (i in 1:(N-1)) {
          APE_estimate[i,] = (pnorm(cof_estimate[i]* X_max + eta[,i] - cof_estimate[i] * X ) - pnorm(cof_estimate[i]* X_min + eta[,i] - cof_estimate[i] * X))
          d_APE_estimate[i,] = (dnorm(cof_estimate[i]* X_max + eta[,i] - cof_estimate[i] * X ) - dnorm(cof_estimate[i]* X_min + eta[,i] - cof_estimate[i] * X))
          dd_APE_estimate[i,] = (-(cof_estimate[i]* X_max + eta[,i] - cof_estimate[i] * X) * dnorm(cof_estimate[i]* X_max + eta[,i] - cof_estimate[i] * X) + (cof_estimate[i]* X_min + eta[,i] - cof_estimate[i] * X)  * dnorm(cof_estimate[i]* X_min + eta[,i] - cof_estimate[i] * X))
        }}else{
          for (i in 1:(N-1)) {
            APE_estimate[i,] = (plogis(cof_estimate[i]* X_max + eta[,i] - cof_estimate[i] * X ) - plogis(cof_estimate[i]* X_min + eta[,i] - cof_estimate[i] * X))
            d_APE_estimate[i,] = (dlogis(cof_estimate[i]* X_max + eta[,i] - cof_estimate[i] * X ) - dlogis(cof_estimate[i]* X_min + eta[,i] - cof_estimate[i] * X))
            dd_APE_estimate = (dlogis(cof_estimate[i]* X_max + eta[,i] - cof_estimate[i] * X ) * ( 1 - 2*plogis(cof_estimate[i] * X_max + eta[,i] - cof_estimate[i] * X ) )
                               - dlogis(cof_estimate[i]* X_min + eta[,i] - cof_estimate[i] * X) * ( 1 - 2 * plogis(cof_estimate[i] * X_min + eta[,i] - cof_estimate[i] * X)))

          }
        }

    }else{
      APE_estimate = matrix(0,N-1,N*(N-1))
      d_APE_estimate = matrix(0,N-1,N*(N-1))
      dd_APE_estimate = matrix(0,N-1,N*(N-1))
      if (model == 'probit'){
        for (i in 1:(N-1)) {
          APE_estimate[i,] = cof_estimate[i] * dnorm(eta[,i])
          d_APE_estimate[i,] = cof_estimate[i] * (-eta[,i]) * dnorm(eta[,i])
          dd_APE_estimate[i,] = cof_estimate[i]*(- dnorm(eta[,i]) + eta[,i]^2 * dnorm(eta[,i]))
        }
      }else{
        for (i in 1:(N-1)) {
          APE_estimate[i,] = cof_estimate[i] * dlogis(eta[,i])
          d_APE_estimate[i,] = cof_estimate[i] * (1 - 2 * plogis(eta[,i])) * dlogis(eta[,i])
          dd_APE_estimate[i,] = cof_estimate[i] * (- 2 * dlogis(eta[,i])) * dlogis(eta[,i]) + cof_estimate[i] * (1 - 2 * plogis(eta[,i])) * dlogis(eta[,i]) * (1-2*plogis(eta[,i]))
        }
      }
    }

    APE_jack[k] = (N-1)*mean(APE_estimate) - (N-2)* mean( apply(APE_estimate,1,mean) )


    index1 = data[,index[1]]
    index2 = data[,index[2]]
    new_index1 <- match(index1, sort(unique(index1)))
    new_index2 <- match(index2, sort(unique(index2)))

    data[,index[1]] =  new_index1
    data[,index[2]] =  new_index2

    X = vector_to_matrix(X, N, ind1 = data[,index[1]], ind2 = data[,index[2]])
    y = vector_to_matrix(y_to, N, ind1 = data[,index[1]], ind2 = data[,index[2]])
    APE_estimate = vector_to_matrix( APE_estimate, N, ind1 = data[,index[1]], ind2 = data[,index[2]])
    d_APE_estimate = vector_to_matrix( d_APE_estimate, N, ind1 = data[,index[1]], ind2 = data[,index[2]])
    dd_APE_estimate = vector_to_matrix( dd_APE_estimate, N, ind1 = data[,index[1]], ind2 = data[,index[2]])

    APE_list[[k]] =  APE_estimate
    d_APE_list[[k]] =  d_APE_estimate

    cov_sum = vector_to_matrix( t(fit$eta_jack), N, ind1 = data[,index[1]], ind2 = data[,index[2]])
    if (model == 'probit'){
      Phi_XB <- pnorm(cov_sum)  # CDF (Φ(Xβ))
      phi_XB <- dnorm(cov_sum)  # PDF (φ(Xβ))
      Phi_XB <- pmax(Phi_XB, 1e-7)
      Phi_XB <- pmin(Phi_XB, 1 - 1e-7)
      phi_XB <- pmax(phi_XB, 1e-7)
      phi_XB <- pmin(phi_XB, 1 - 1e-7)
      dd_F_fix = -cov_sum * phi_XB

      H = phi_XB/(Phi_XB*(1-Phi_XB))
      small_w = H * phi_XB
    }else{

      Phi_XB <- plogis(cov_sum)  # CDF (Φ(Xβ))
      phi_XB <- dlogis(cov_sum)  # PDF (φ(Xβ))
      dd_F_fix = phi_XB * ( 1 - 2 * Phi_XB)
      Phi_XB[row(Phi_XB) == col(Phi_XB)] = 0
      phi_XB[row(phi_XB) == col(phi_XB)] = 0
      Phi_XB <- pmax(Phi_XB, 1e-7)
      Phi_XB <- pmin(Phi_XB, 1 - 1e-7)
      phi_XB <- pmax(phi_XB, 1e-7)
      phi_XB <- pmin(phi_XB, 1 - 1e-7)
      #### Another way to calculate the bias corrected
      # another way to calculate the analytical corrected estimate
      H = matrix(1, N, N)
      diag(H) = 0
      small_w = phi_XB
    }

    # B_hat
    X_into = matrix_to_panel_df(-d_APE_estimate/small_w)
    weight = matrix_to_panel_df(small_w)$X
    Psi = get_weighted_fe_projection(X_into$X, weight, X_into$id, X_into$time)
    Psi_list[[k]]= vector_to_matrix(Psi, N, ind1 = X_into$id, ind2 = X_into$time)
    tilde_Psi_list[[k]] = -d_APE_estimate/small_w - Psi_list[[k]]
    diag(tilde_Psi_list[[k]]) = 0

    # residual X
    X_into = matrix_to_panel_df(X)
    weight = matrix_to_panel_df(small_w)$X
    re = get_weighted_fe_projection(X_into$X, weight, X_into$id, X_into$time)
    re_matrix = vector_to_matrix(re, N, ind1 = X_into$id, ind2 = X_into$time)
    tilde_X_list[[k]] = X - re_matrix
    re_matrix_array[,,k] = re_matrix

  }

  tilde_X_array <- array(unlist(tilde_X_list), dim = c(N, N, K))
  Psi_array <- array(unlist(Psi_list), dim = c(N, N, K))
  tilde_Psi_array <- array(unlist(tilde_Psi_list), dim = c(N, N, K))
  APE_array <- array(unlist(APE_list), dim = c(N, N, K))
  d_APE_array <- array(unlist(d_APE_list), dim = c(N, N, K))

  J = matrix(0, K,K)
  for (k in 1:K) {
    X = X_start[,k]
    if (is.numeric(X) && length(unique(X)) == 2){
      X_max = as.numeric(max(unique(X)))
      X_min = as.numeric(min(unique(X)))
      if (model == 'probit'){
        f1_own = dnorm(cof_estimate * X_max + eta - cof_estimate * X )
      }else{
        f1_own = dlogis(cof_estimate * X_max + eta - cof_estimate * X )
      }
      for (j in 1:K){
        diag(tilde_X_array[,,j]) = 0
        diag(d_APE_array[,,k]) = 0
        J[j, k] <- - sum(re_matrix_array[,,j] * d_APE_array[,,k]) / (N*N-1)
      }
      J[k, k] <- sum(f1_own) / (N*N-1) + J[k, k]
      J[- k, k] <- colSums(X_start[, - k, drop = FALSE] * matrix_to_panel_df(as.matrix(d_APE_array[,,k]))$X) / (N*N-1) + J[- k, k]
    }else{
      for (j in 1:K){
        diag(tilde_X_array[,,j]) = 0
        diag(d_APE_array[,,k]) = 0
        J[j, k] <- sum(tilde_X_array[,,j] * d_APE_array[,,k]) / (N*N-1)
      }
      J[k, k] <- sum(phi_XB) / (N*N-1) + J[k, k]
    }
  }

  # Step 4: Initialize sum matrix for outer products
  W_hat <- matrix(0, nrow = K, ncol = K)
  D_beta_APE <- 0
  tau = array(0, dim = c(N, N, K))
  # Step 5: For each position (i, j), extract vector and compute outer product
  for (i in 1:N) {
    for (j in 1:N) {
      if(j!=i){
        v_ij <- tilde_X_array[i, j, ]             # Vector of length p
        W_hat <- W_hat + (1 / ((N-1)*N)) * small_w[i,j] * v_ij %*% t(v_ij)  # Outer product and sum
        D_beta_APE =  D_beta_APE +  (1 / ((N-1)*N)) * d_APE_array[i,j, ] * v_ij
      }
    }
  }

  for (i in 1:N) {
    for (j in 1:N) {
      if(j!=i){
        tau[i,j,] = as.matrix(t(J) %*% solve(W_hat)) %*% as.matrix(H[i,j] * (y[i,j] - Phi_XB[i,j]) * tilde_X_array[i,j,]) - Psi_array[i,j,] * H[i,j] * (y[i,j] - Phi_XB[i,j])
      }
    }
  }


  for (k in 1:K) {
    diag(tau[,,k]) <- 0
  }
  tau_sum = 0
  for (i in 1:N) {
    for (t in 1:N) {
      v <- tau[i, t, ]            # tau_{it} as a vector
      tau_sum <- tau_sum + tcrossprod(v)  # v %*% t(v)
    }
  }

  se =  sqrt( diag( tau_sum ) )/(N*(N-1))

  res = list(APE = APE_jack, se = se)

  return(res)

}


#' @export
get_APE_MLE <- function(y, X, N, index, data, fit,  model = 'probit'){
  K = dim(X)[2]
  APE = rep(0, K)
  tilde_X_list <- vector("list", K)
  tilde_Psi_list <- vector("list", K)
  APE_list <- vector("list", K)
  d_APE_list <- vector("list", K)
  Psi_list = vector("list", K)
  re_matrix_array = array(0, dim = c(N,N,K))
  W_hat <- matrix(0, nrow = K, ncol = K)
  X_start = X
  y_to = y
  for (k in 1:K) {
    X = X_start[,k]
    cof_estimate = fit$est_MLE[k]
    eta = t(fit$eta_MLE)
    if (is.numeric(X) && length(unique(X)) == 2) {
      X_max = as.numeric(max(unique(X)))
      X_min = as.numeric(min(unique(X)))
      if (model == 'probit'){
        APE_estimate = (pnorm(cof_estimate * X_max + eta - cof_estimate * X ) - pnorm(cof_estimate* X_min + eta - cof_estimate * X))
        d_APE_estimate = (dnorm(cof_estimate * X_max + eta - cof_estimate * X ) - dnorm(cof_estimate* X_min +eta - cof_estimate * X) )
        dd_APE_estimate = (-(cof_estimate * X_max + eta - cof_estimate * X) * dnorm(cof_estimate* X_max + eta - cof_estimate * X ) + (cof_estimate* X_min + eta - cof_estimate * X)  * dnorm(cof_estimate* X_min + eta - cof_estimate * X))
      }else{
        APE_estimate = (plogis(cof_estimate * X_max + eta - cof_estimate * X ) - plogis(cof_estimate* X_min + eta - cof_estimate * X))
        d_APE_estimate = (dlogis(cof_estimate * X_max + eta - cof_estimate * X ) - dlogis(cof_estimate* X_min +eta - cof_estimate * X) )
        dd_APE_estimate = (dlogis(cof_estimate* X_max + eta - cof_estimate * X ) * ( 1 - 2*plogis(cof_estimate* X_max + eta - cof_estimate * X ) )
                           - dlogis(cof_estimate* X_min + eta - cof_estimate * X) * ( 1 - 2 * plogis(cof_estimate* X_min + eta - cof_estimate * X)))

      }
    }else{
      if (model == 'probit'){
        APE_estimate = cof_estimate * dnorm(eta)
        d_APE_estimate = cof_estimate * (-eta) * dnorm(eta)
        dd_APE_estimate = cof_estimate * (- dnorm(eta) + eta^2 * dnorm(eta))
      }else{
        APE_estimate = cof_estimate * dlogis(eta)
        d_APE_estimate = cof_estimate * (1 - 2 * plogis(eta)) * dlogis(eta)
        dd_APE_estimate = cof_estimate* ( - 2 * dlogis(eta)) * dlogis(eta) + cof_estimate * (1 - 2 * plogis(eta)) * dlogis(eta) * (1-2*plogis(eta))
      }
    }


    APE[k] =  mean(APE_estimate)

    index1 = data[,index[1]]
    index2 = data[,index[2]]
    new_index1 <- match(index1, sort(unique(index1)))
    new_index2 <- match(index2, sort(unique(index2)))

    data[,index[1]] =  new_index1
    data[,index[2]] =  new_index2

    X = vector_to_matrix(X, N, ind1 = data[,index[1]], ind2 = data[,index[2]])
    y = vector_to_matrix(y_to, N, ind1 = data[,index[1]], ind2 = data[,index[2]])
    APE_estimate = vector_to_matrix( APE_estimate, N, ind1 = data[,index[1]], ind2 = data[,index[2]])
    d_APE_estimate = vector_to_matrix( d_APE_estimate, N, ind1 = data[,index[1]], ind2 = data[,index[2]])
    dd_APE_estimate = vector_to_matrix( dd_APE_estimate, N, ind1 = data[,index[1]], ind2 = data[,index[2]])

    APE_list[[k]] =  APE_estimate
    d_APE_list[[k]] =  d_APE_estimate

    cov_sum = vector_to_matrix( t(eta), N, ind1 = data[,index[1]], ind2 = data[,index[2]])

    if (model == 'probit'){
      Phi_XB <- pnorm(cov_sum)  # CDF (Φ(Xβ))
      phi_XB <- dnorm(cov_sum)  # PDF (φ(Xβ))
      Phi_XB <- pmax(Phi_XB, 1e-7)
      Phi_XB <- pmin(Phi_XB, 1 - 1e-7)
      phi_XB <- pmax(phi_XB, 1e-7)
      phi_XB <- pmin(phi_XB, 1 - 1e-7)
      dd_F_fix = -cov_sum * phi_XB

      H = phi_XB/(Phi_XB*(1-Phi_XB))
      small_w = H * phi_XB
    }else{

      Phi_XB <- plogis(cov_sum)  # CDF (Φ(Xβ))
      phi_XB <- dlogis(cov_sum)  # PDF (φ(Xβ))
      dd_F_fix = phi_XB * ( 1 - 2 * Phi_XB)
      Phi_XB[row(Phi_XB) == col(Phi_XB)] = 0
      phi_XB[row(phi_XB) == col(phi_XB)] = 0
      Phi_XB <- pmax(Phi_XB, 1e-7)
      Phi_XB <- pmin(Phi_XB, 1 - 1e-7)
      phi_XB <- pmax(phi_XB, 1e-7)
      phi_XB <- pmin(phi_XB, 1 - 1e-7)
      #### Another way to calculate the bias corrected
      # another way to calculate the analytical corrected estimate
      H = matrix(1, N, N)
      diag(H) = 0
      small_w = phi_XB
    }

    # B_hat
    X_into = matrix_to_panel_df(-d_APE_estimate/small_w)
    weight = matrix_to_panel_df(small_w)$X
    Psi = get_weighted_fe_projection(X_into$X, weight, X_into$id, X_into$time)
    Psi_list[[k]]= vector_to_matrix(Psi, N, ind1 = X_into$id, ind2 = X_into$time)
    tilde_Psi_list[[k]] = -d_APE_estimate/small_w - Psi_list[[k]]
    diag(tilde_Psi_list[[k]]) = 0

    # residual X
    X_into = matrix_to_panel_df(X)
    weight = matrix_to_panel_df(small_w)$X
    re = get_weighted_fe_projection(X_into$X, weight, X_into$id, X_into$time)
    re_matrix = vector_to_matrix(re, N, ind1 = X_into$id, ind2 = X_into$time)
    tilde_X_list[[k]] = X - re_matrix
    re_matrix_array[,,k] = re_matrix

  }

  tilde_X_array <- array(unlist(tilde_X_list), dim = c(N, N, K))
  Psi_array <- array(unlist(Psi_list), dim = c(N, N, K))
  tilde_Psi_array <- array(unlist(tilde_Psi_list), dim = c(N, N, K))
  APE_array <- array(unlist(APE_list), dim = c(N, N, K))
  d_APE_array <- array(unlist(d_APE_list), dim = c(N, N, K))

  J = matrix(0, K,K)
  for (k in 1:K) {
    X = X_start[,k]
    if (is.numeric(X) && length(unique(X)) == 2){
      X_max = as.numeric(max(unique(X)))
      X_min = as.numeric(min(unique(X)))
      if (model == 'probit'){
        f1_own = dnorm(cof_estimate * X_max + eta - cof_estimate * X )
      }else{
        f1_own = dlogis(cof_estimate* X_max + eta - cof_estimate * X )
      }
      for (j in 1:K){
        diag(tilde_X_array[,,j]) = 0
        diag(d_APE_array[,,k]) = 0
        J[j, k] <- - sum(re_matrix_array[,,j] * d_APE_array[,,k]) / (N*N-1)
      }
      J[k, k] <- sum(f1_own) / (N*N-1) + J[k, k]
      J[- k, k] <- colSums(X_start[, - k, drop = FALSE] * matrix_to_panel_df(as.matrix(d_APE_array[,,k]))$X) / (N*N-1) + J[- k, k]
    }else{
      for (j in 1:K){
        diag(tilde_X_array[,,j]) = 0
        diag(d_APE_array[,,k]) = 0
        J[j, k] <- sum(tilde_X_array[,,j] * d_APE_array[,,k]) / (N*N-1)
      }
      J[k, k] <- sum(phi_XB) / (N*N-1) + J[k, k]
    }
  }

  # Step 4: Initialize sum matrix for outer products
  W_hat <- matrix(0, nrow = K, ncol = K)
  D_beta_APE <- 0
  tau = array(0, dim = c(N, N, K))
  # Step 5: For each position (i, j), extract vector and compute outer product
  for (i in 1:N) {
    for (j in 1:N) {
      if(j!=i){
        v_ij <- tilde_X_array[i, j, ]             # Vector of length p
        W_hat <- W_hat + (1 / ((N-1)*N)) * small_w[i,j] * v_ij %*% t(v_ij)  # Outer product and sum
        D_beta_APE =  D_beta_APE +  (1 / ((N-1)*N)) * d_APE_array[i,j, ] * v_ij
      }
    }
  }

  for (i in 1:N) {
    for (j in 1:N) {
      if(j!=i){
        tau[i,j,] = as.matrix(t(J) %*% solve(W_hat)) %*% as.matrix(H[i,j] * (y[i,j] - Phi_XB[i,j]) * tilde_X_array[i,j,]) - Psi_array[i,j,] * H[i,j] * (y[i,j] - Phi_XB[i,j])
      }
    }
  }


  for (k in 1:K) {
    diag(tau[,,k]) <- 0
  }
  tau_sum = 0
  for (i in 1:N) {
    for (t in 1:N) {
      v <- tau[i, t, ]            # tau_{it} as a vector
      tau_sum <- tau_sum + tcrossprod(v)  # v %*% t(v)
    }
  }

  se =  sqrt( diag( tau_sum ) )/(N*(N-1))

  res = list(APE = APE, se = se)

  return(res)

}
