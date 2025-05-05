#' @export
#'
APE_se = function(fit, est, N, X, y, APE, model = 'probit'){


  cov_sum_1 = fit$X_origin[,1] * est[1]
  cov_sum_2 = fit$X_origin[,-1] %*% est[-1]
  cov_APE = matrix(est[1] + fit$X_origin[,-1] %*% est[-1], N-1, N)
  cov_APE_minus = matrix(-est[1] + fit$X_origin[,-1] %*% est[-1], N-1, N)
  cov_sum = matrix(cov_sum_1 + cov_sum_2, N-1, N)

  cov_APE = shift_lower_triangle_and_add_zero_diag(cov_APE)
  cov_APE_minus = shift_lower_triangle_and_add_zero_diag(cov_APE_minus)
  cov_sum = shift_lower_triangle_and_add_zero_diag(cov_sum)

  cov_sum[row(cov_sum) == col(cov_sum)] = 0
  cov_APE[row(cov_APE) == col(cov_APE)] = 0
  cov_APE_minus[row(cov_APE_minus) == col(cov_APE_minus)] = 0

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
  # d_fix_loss = y - exp(cov_sum)/( 1 + exp(cov_sum))
  d_fix_loss =  derivative_ingredients$d_fix_loss
  d_beta_loss = derivative_ingredients$d_beta_loss
  #d_fix_fix_loss = - exp(cov_sum)/( 1 + exp(cov_sum))^2
  d_fix_fix_loss = derivative_ingredients$d_fix_fix_loss# d_beta_loss = X * (y - exp(cov_sum)/( 1 + exp(cov_sum)))
  # d_beta_beta_loss = X * (- exp(cov_sum)/( 1 + exp(cov_sum))^2) * X
  d_beta_beta_loss = derivative_ingredients$d_beta_beta_loss# d_beta_fix_loss = X * ( - exp(cov_sum)/( 1 + exp(cov_sum))^2)
  d_beta_fix_loss = derivative_ingredients$d_beta_fix_loss


  d_fix_APE =  (1/sqrt(2*pi) * exp( -cov_APE^2 / 2 ) - 1/sqrt(2*pi) * exp( -cov_APE_minus^2 / 2 ))/2
  d_beta_APE = (1/sqrt(2*pi) * exp( -cov_APE^2 / 2 ) + 1/sqrt(2*pi) * exp( -cov_APE_minus^2 / 2 ))/2


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
  W_bar =   -(1 / ((N-1)*N)) * (sum(  d_beta_beta_loss  - d_fix_fix_loss * the * the  ) - sum(diag(d_beta_beta_loss  - d_fix_fix_loss * the * the )) )


  # phi matrix
  phi = matrix(0, N, N)
  for (i in 1:N) {
    for (j in 1:N) {
      temp_sum <- 0
      for (s in 1:N) {
        for (t in 1:N) {
          if (t != s) {
            temp_sum <- temp_sum + (-1/(N)) * (Hessian_a_a[i, s] + Hessian_g_a[j, s] + Hessian_a_g[i, t] + Hessian_g_g[j, t]) * d_fix_APE[s, t]
          }
        }
      }
      phi[i, j] <- temp_sum
    }
  }


  D_beta_delta = (1/((N-1)*N)) * (sum(d_beta_APE - the * d_fix_APE) - sum(diag(d_beta_APE - the * d_fix_APE)))

  D_beta_loss = d_beta_loss - d_fix_loss*the

  tau = as.numeric(solve(W_bar) * D_beta_delta)  * D_beta_loss - phi * d_fix_loss

  APE_residual = APE - colMeans(APE) # indentical assumption


  part_1 = sum(APE_residual%*%APE_residual) - sum(diag(APE_residual%*%APE_residual))
  part_2 = compute_sum(APE_residual)
  part_3 = sum(tau*tau) - sum(diag(tau*tau))


  return( sqrt((part_1 + part_2 + part_3)/(N^2*(N-1)^2)) )
}

#' @export
se_formula = function(fit, est, N, X, y, model = 'probit'){

  cov_sum_1 = fit$X_origin[,1] * est[1]
  cov_sum_2 = fit$X_origin[,-1] %*% est[-1]
  cov_sum = matrix(cov_sum_1 + cov_sum_2, N-1, N)
  cov_sum = shift_lower_triangle_and_add_zero_diag(cov_sum)
  cov_sum[row(cov_sum) == col(cov_sum)] = 0

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
  # d_fix_loss = y - exp(cov_sum)/( 1 + exp(cov_sum))
  d_fix_loss =  derivative_ingredients$d_fix_loss
  d_beta_loss = derivative_ingredients$d_beta_loss
  #d_fix_fix_loss = - exp(cov_sum)/( 1 + exp(cov_sum))^2
  d_fix_fix_loss = derivative_ingredients$d_fix_fix_loss# d_beta_loss = X * (y - exp(cov_sum)/( 1 + exp(cov_sum)))
  # d_beta_beta_loss = X * (- exp(cov_sum)/( 1 + exp(cov_sum))^2) * X
  d_beta_beta_loss = derivative_ingredients$d_beta_beta_loss# d_beta_fix_loss = X * ( - exp(cov_sum)/( 1 + exp(cov_sum))^2)
  d_beta_fix_loss = derivative_ingredients$d_beta_fix_loss


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

  se_version1 = sqrt(solve(W_hat) * Omega_hat * solve(W_hat)/N^2)
  res = list(se = se_version1)

  return( res )

}

#' @export
se_formula_corrected = function(y, X, N, data, index, est, model = 'probit'){
  # order the data to match the index and the summy variable of the fixed effects
  data = data.frame(y = y, X = X, data[,index[1]], data[,index[2]])
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


  # prepare the dummy variable for fixed effects
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


  # prepare the design matrix with dummy variables
  X_design = cbind(X, fix)
  X_design = apply(X_design, 2, as.numeric)


  # calculate the X'beta + pi, formed as a matrix
  cov_sum_1 = X_design[,1] * est[1]
  cov_sum_2 = X_design[,-1] %*% est[-1]
  cov_sum = matrix(cov_sum_1 + cov_sum_2, N-1, N)
  cov_sum = shift_lower_triangle_and_add_zero_diag(cov_sum)
  cov_sum[row(cov_sum) == col(cov_sum)] = 0


  # X is a matrix of a single covariate, and the same for y
  X = vector_to_matrix(X, N, ind1 = data[,index[1]], ind2 = data[,index[2]])
  y = vector_to_matrix(y, N, ind1 = data[,index[1]], ind2 = data[,index[2]])
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
  # d_fix_loss = y - exp(cov_sum)/( 1 + exp(cov_sum))
  d_fix_loss =  derivative_ingredients$d_fix_loss
  d_beta_loss = derivative_ingredients$d_beta_loss
  #d_fix_fix_loss = - exp(cov_sum)/( 1 + exp(cov_sum))^2
  d_fix_fix_loss = derivative_ingredients$d_fix_fix_loss# d_beta_loss = X * (y - exp(cov_sum)/( 1 + exp(cov_sum)))
  # d_beta_beta_loss = X * (- exp(cov_sum)/( 1 + exp(cov_sum))^2) * X
  d_beta_beta_loss = derivative_ingredients$d_beta_beta_loss# d_beta_fix_loss = X * ( - exp(cov_sum)/( 1 + exp(cov_sum))^2)
  d_beta_fix_loss = derivative_ingredients$d_beta_fix_loss


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
  x1 = matrix_to_panel_df(d_beta_fix_loss)
  x2 = matrix_to_panel_df(d_fix_fix_loss)
  to_in = x1$X/x2$X
  to_in[is.nan(to_in)] <- 0
  weight = -x2$X
  weight[which(weight==0)] = 1
  re = get_weighted_projection_fitted_exclude_t_eq_i(to_in, weight, x1$id, x1$time)
  re_matrix = matrix(re, N-1, N)
  re_matrix = shift_lower_triangle_and_add_zero_diag(re_matrix)
  tilde_X = X - re_matrix
  W_hat_another =   (1 / ((N-1)*N)) * (sum(  small_w  * tilde_X  * tilde_X   ) - sum(diag( small_w  * tilde_X  * tilde_X  )) )

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

  # matrix from huges
  d_beta_beta_big_loss = (1/(N-1))*(sum(d_beta_beta_loss) - sum(diag(d_beta_beta_loss)))
  d_beta_fix_big_loss = (1/(N-1))*(sum(d_beta_fix_loss) - sum(diag(d_beta_fix_loss)))
  d_fix_beta_big_loss = (1/(N-1))*(sum(d_beta_fix_loss) - sum(diag(d_beta_fix_loss)))
  d_fix_2_big_loss = (1/(N-1))*(sum(d_fix_fix_loss) - sum(diag(d_fix_fix_loss)))

  # W_hat
  W_hat =   -(1 / N) * (  d_beta_beta_big_loss  - d_beta_fix_big_loss * solve(d_fix_2_big_loss) * d_beta_fix_big_loss  )
  D_beta_loss = d_beta_loss - d_fix_loss*the
  A = (D_beta_loss + t(D_beta_loss))^2

  # Omega_hat
  Omega_hat =  1/(N*(N-1))* sum(A[lower.tri(A)])

  res = list(se_weidner = sqrt(solve(W_hat)/(N*(N-1))), se_no_MLE =sqrt(solve(W_hat) * Omega_hat * solve(W_hat)/(N*(N-1))),
             se_huges = sqrt(solve(W_hat)/(N*(N-1))), se_huges_no_MLE = sqrt(solve(W_hat) * Omega_hat * solve(W_hat)/(N*(N-1))))


  return(res)

}

# #' @export
# se_formula_corrected = function(y, X, N, data, index, est, model = 'probit'){
#   # order the data to match the index and the summy variable of the fixed effects
#   data = data.frame(y = y, X = X, data[,index[1]], data[,index[2]])
#   if(is.null(colnames(X)) == 1){
#     colnames(data)[(dim(data)[2]-1):dim(data)[2]] = index
#     data = data[order(data[,index[2]], data[,index[1]]),]
#     p = dim(X)[2]
#     y = data$y
#     if (p==1){
#       X = data[,2]
#     }else{
#       X = data[,(2:p+1)]
#     }
#   }else{
#     colnames(data) = c('y', colnames(X), index)
#     data = data[order(data[,index[2]], data[,index[1]]),]
#     p = dim(X)[2]

#     y = data$y
#     X = data[,colnames(X)]
#   }


#   # prepare the dummy variable for fixed effects
#   fix_effect = matrix(0, N*N, N + N )
#   for (t in seq(N)) {
#     for (i in seq(N)) {
#       alpha_in = rep(0, N)
#       gamma_in = rep(0, N)
#       alpha_in[i] = 1
#       gamma_in[t] = 1
#       fix_effect[i + (t - 1) * N,] = c(alpha_in - alpha_in[1], gamma_in + alpha_in[1])
#     }
#   }
#   drop_index = NULL
#   for (i in seq(N)){
#     drop_index = cbind(drop_index, i + (i - 1) * N)
#   }
#   fix = fix_effect[-drop_index,]


#   # prepare the design matrix with dummy variables
#   X_design = cbind(X, fix)
#   X_design = apply(X_design, 2, as.numeric)


#   # calculate the X'beta + pi, formed as a matrix
#   cov_sum_1 = X_design[,1] * est[1]
#   cov_sum_2 = X_design[,-1] %*% est[-1]
#   cov_sum = matrix(cov_sum_1 + cov_sum_2, N-1, N)
#   cov_sum = shift_lower_triangle_and_add_zero_diag(cov_sum)
#   cov_sum[row(cov_sum) == col(cov_sum)] = 0


#   # X is a matrix of a single covariate, and the same for y
#   X = vector_to_matrix(X, N, ind1 = data[,index[1]], ind2 = data[,index[2]])
#   y = vector_to_matrix(y, N, ind1 = data[,index[1]], ind2 = data[,index[2]])
#   y[row(y) == col(y)] = 0


#   # CDF (Φ(Xβ)) and PDF (φ(Xβ))
#   Phi_XB <- pnorm(cov_sum)
#   phi_XB <- dnorm(cov_sum)
#   dd_F_fix = -cov_sum * phi_XB
#   Phi_XB[row(Phi_XB) == col(Phi_XB)] = 0
#   phi_XB[row(phi_XB) == col(phi_XB)] = 0
#   Phi_XB <- pmax(Phi_XB, 1e-9)
#   Phi_XB <- pmin(Phi_XB, 1 - 1e-9)


#   # preparation for ingredients
#   # d_fix_loss = y - exp(cov_sum)/( 1 + exp(cov_sum))
#   d_fix_loss =  y * (phi_XB) / Phi_XB  +  (1-y)*(-phi_XB)/(1-Phi_XB) # (phi_XB * (y - Phi_XB) / (Phi_XB * (1 - Phi_XB)))
#   #d_fix_2_loss = - exp(cov_sum)/( 1 + exp(cov_sum))^2
#   d_fix_2_loss = (dd_F_fix * Phi_XB * ( 1 -Phi_XB ) - phi_XB*(1-2*Phi_XB)*phi_XB )*(y-Phi_XB)/((Phi_XB*(1-Phi_XB)))^2 - phi_XB^2  / ((Phi_XB*(1-Phi_XB)))
#   # d_beta_loss = X * (y - exp(cov_sum)/( 1 + exp(cov_sum)))
#   d_beta_loss = X * (phi_XB * (y - Phi_XB) / (Phi_XB * (1 - Phi_XB)))
#   # d_beta_beta_loss = X * (- exp(cov_sum)/( 1 + exp(cov_sum))^2) * X
#   d_beta_beta_loss = (dd_F_fix * Phi_XB * ( 1 -Phi_XB ) - phi_XB*(1-2*Phi_XB)*phi_XB )* (X^2 *(y-Phi_XB)) /((Phi_XB*(1-Phi_XB)))^2 - X^2 * phi_XB^2  / ((Phi_XB*(1-Phi_XB)))
#   # d_beta_fix_loss = X * ( - exp(cov_sum)/( 1 + exp(cov_sum))^2)
#   d_beta_fix_loss = (dd_F_fix * Phi_XB * ( 1 -Phi_XB ) - phi_XB*(1-2*Phi_XB)*phi_XB )* (X *(y-Phi_XB)) /((Phi_XB*(1-Phi_XB)))^2 - X * phi_XB^2  / ((Phi_XB*(1-Phi_XB)))
#   d_beta_beta_big_loss = (1/(N-1))*(sum(d_beta_beta_loss) - sum(diag(d_beta_beta_loss)))
#   d_beta_fix_big_loss = (1/(N-1))*(sum(d_beta_fix_loss) - sum(diag(d_beta_fix_loss)))
#   d_fix_beta_big_loss = (1/(N-1))*(sum(d_beta_fix_loss) - sum(diag(d_beta_fix_loss)))
#   d_fix_2_big_loss = (1/(N-1))*(sum(d_fix_2_loss) - sum(diag(d_fix_2_loss)))


#   # Hessian matrix, based on huges's paper
#   #H_a_a = diag(-rowSums(d_fix_2_loss))/sqrt(N*(N-1))
#   #H_g_g = diag(-colSums(d_fix_2_loss))/sqrt(N*(N-1))

#   H_a_a = matrix(0, N, N)
#   H_g_g = matrix(0, N, N)
#   for (i in 1:nrow(d_fix_2_loss)) {
#     H_a_a[i, i] <- -(sum(d_fix_2_loss[i, ]) - d_fix_2_loss[i, i])/(N-1)  # sum of off-diagonal elements in row i
#   }
#   for (i in 1:nrow(d_fix_2_loss)) {
#     H_g_g[i, i] <- -(sum(d_fix_2_loss[-i, i]))/(N-1) # sum of off-diagonal elements in row i
#   }

#   H_a_g = -d_fix_2_loss/((N-1))
#   Hessian_bar =   cbind(rbind(H_a_a, t(H_a_g)), rbind(H_a_g, H_g_g))  + c(rep(1,N), rep(-1,N)) %*% t( c(rep(1,N), rep(-1,N)) )/N
#   Hessain_inverse = solve(Hessian_bar)
#   Hessian_a_a = Hessain_inverse[1:(N), 1:(N)]
#   Hessian_g_a = Hessain_inverse[(N+1):(N+N), 1:(N)]
#   Hessian_a_g = Hessain_inverse[1:(N), (N+1):(N+N)]
#   Hessian_g_g = Hessain_inverse[(N+1):(N+N), (N+1):(N+N)]

#   # matrix 'The', \based on huges's paper
#   the = matrix(0, N, N)
#   for (i in 1:N) {
#     for (j in 1:N) {
#       temp_sum <- 0
#       for (s in 1:N) {
#         for (t in 1:N) {
#           if (t != s) {
#             temp_sum <- temp_sum + (-1/(N)) * (Hessian_a_a[i, s] + Hessian_g_a[j, s] + Hessian_a_g[i, t] + Hessian_g_g[j, t]) * d_beta_fix_loss[s, t]
#           }
#         }
#       }
#       the[i, j] <- temp_sum
#     }
#   }

#   # W_hat based on Iva ́ n2016
#   W_hat =   -(1 / ((N-1)*N)) * (sum(  d_beta_beta_loss  - d_fix_2_loss * the * the  ) - sum(diag(d_beta_beta_loss  - d_fix_2_loss * the * the )) ) # -(1 / N) * (  d_beta_beta_big_loss  - d_beta_fix_big_loss * solve(d_fix_2_big_loss) * d_beta_fix_big_loss  )

#   # Omega_hat based on Iva ́ n2016
#   D_beta_loss = d_beta_loss - d_fix_loss*the
#   Omega_hat <- 0
#   for (i in 1:N) {
#     for (j in setdiff(1:N, i)) {
#       for (k in setdiff(1:N, i)) {
#         Omega_hat <- Omega_hat + D_beta_loss[i, j] * D_beta_loss[i, k]/(N*(N-1))
#       }
#     }
#   }

#   se_version1 = sqrt(solve(W_hat) * Omega_hat * solve(W_hat)/N^2)

#   # matrix from huges
#   # W_hat
#   W_hat =   -(1 / N) * (  d_beta_beta_big_loss  - d_beta_fix_big_loss * solve(d_fix_2_big_loss) * d_beta_fix_big_loss  )

#   D_beta_loss = d_beta_loss - d_fix_loss*the

#   A = (D_beta_loss + t(D_beta_loss))^2

#   Omega_hat =  1/(N*(N-1))* sum(A[lower.tri(A)])

#   se_version2 = sqrt(solve(W_hat) * Omega_hat * solve(W_hat)/N^2)

#   res = list(se_1 = se_version1, se_2 = se_version2)
#   return( res)

# }


vector_to_matrix <- function(v, N, ind1, ind2) {
  # Determine matrix size
  nrow <- N
  ncol <- N

  # Initialize zero matrix
  M <- matrix(0, nrow = nrow, ncol = ncol)

  # Assign values
  for (i in seq_along(v)) {
    M[ind1[i], ind2[i]] <- v[i]
  }

  return(M)
}

compute_sum <- function(A, N) {
  N = dim(A)[1]
  # Ensure that a and b are matrices of size N x N
  result <- sum(sapply(1:N, function(i) {
    sum(sapply(setdiff(1:N, i), function(t) {
      sum(sapply(setdiff(1:N, c(i, t)), function(j) {
        A[i, t] * A[j, t]
      }))
    }))
  }))
  return(result)
}


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

compute_derivatives <- function(eta, y, X, model = 'probit') {

  z <- eta
  phi_z <- dnorm(z)
  Phi_z <- pnorm(z)
  Phi_z[row(Phi_z) == col(Phi_z)] = 0
  phi_z[row(phi_z) == col(phi_z)] = 0

  # Avoid division by 0
  Phi_z <- pmax(Phi_z, 1e-10)
  Phi_z <- pmin(Phi_z, 1 - 1e-9)

  phi_z <- pmax(phi_z, 1e-10)
  phi_z <- pmin(phi_z, 1 - 1e-9)

  one_minus_Phi_z <- pmax(1 - Phi_z, 1e-10)
  dd_F_fix = - eta * phi_z
  ddd_F_fix = eta^2 * phi_z - phi_z


  score1 <- (y / Phi_z - (1 - y) / one_minus_Phi_z) * phi_z

  # Second derivative w.r.t. z
  term1 <- y * ((-z * phi_z) / Phi_z - (phi_z^2) / (Phi_z^2))
  term2 <- (1 - y) * ((-z * phi_z) / one_minus_Phi_z + (phi_z^2) / (one_minus_Phi_z^2))
  score2 <- term1 - term2

  # Third derivative w.r.t. z
  term1 = y * ( (ddd_F_fix*Phi_z - phi_z*dd_F_fix) / Phi_z^2  - (2*phi_z * dd_F_fix * Phi_z^2 - 2 * Phi_z * phi_z^3) / Phi_z^4   )
  term2 = (1-y) * ( (ddd_F_fix*(1-Phi_z) + phi_z*dd_F_fix) / (1-Phi_z)^2  + (2*phi_z * dd_F_fix * (1-Phi_z)^2 + 2 * (1-Phi_z) * phi_z^3) / (1-Phi_z)^4   )
  score3 <- term1 - term2

  if (model == 'probit'){

    d_beta_loss = X * score1
    d_fix_loss = score1
    # Second-order derivatives
    d_beta_beta_loss <- score2 * X^2      # d²ℓ / dβ²
    d_fix_fix_loss  <- score2            # d²ℓ / dπ²
    d_beta_fix_loss <- score2 * X       # d²ℓ / dβ dπ

    d_beta_beta_beta_loss <- score3 * X^3      # d²ℓ / dβ²
    d_fix_fix_fix_loss  <- score3           # d²ℓ / dπ²
    d_beta_fix_fix_loss <- score3 * X       # d²ℓ / dβ dπ

  }

  diag(d_fix_loss) = 0
  diag(d_beta_loss) = 0
  diag(d_beta_beta_loss) = 0
  diag(d_fix_fix_loss) = 0
  diag(d_beta_fix_loss) = 0
  diag(d_beta_beta_loss) = 0
  diag(d_beta_beta_beta_loss) = 0
  diag(d_beta_fix_fix_loss) = 0
  diag(d_fix_fix_fix_loss) = 0

  res = list(d_fix_loss = d_fix_loss, d_fix_fix_loss = d_fix_fix_loss,
             d_beta_loss = d_beta_loss, d_beta_beta_loss = d_beta_beta_loss, d_beta_fix_loss = d_beta_fix_loss,
             d_beta_beta_beta_loss = d_beta_beta_beta_loss, d_beta_fix_fix_loss = d_beta_fix_fix_loss, d_fix_fix_fix_loss = d_fix_fix_fix_loss)
  return(res)
}

#' @importFrom fixest feols
get_weighted_projection_fitted <- function(X, weight, id, time) {
  # Input validation: check if all vectors have the same length
  if (length(X) != length(weight) || length(X) != length(id) || length(X) != length(time)) {
    stop("All input vectors (X, weight, id, time) must have the same length.")
  }

  # Replace NaN in X with 0
  X[is.nan(X)] <- 0

  # Create data frame from input vectors
  df <- data.frame(Xit = X, weight = weight, id = id, time = time)

  # Run the weighted two-way fixed effects regression
  model <- feols(Xit ~ 1 | id + time, data = df, weights = ~weight)

  # Return fitted (predicted) values from the model
  return(fitted(model))
}

#' @importFrom fixest feols
get_weighted_projection_fitted_exclude_t_eq_i <- function(X, weight, id, time) {
  # Input validation: check if all vectors have the same length
  if (length(X) != length(weight) || length(X) != length(id) || length(X) != length(time)) {
    stop("All input vectors (X, weight, id, time) must have the same length.")
  }

  # Replace NaN in X with 0
  X[is.nan(X)] <- 0

  # Create data frame from input vectors
  df <- data.frame(Xit = X, weight = weight, id = id, time = time)

  # Exclude rows where time == id (t != i)
  df <- df[df$id != df$time, ]

  # Run the weighted two-way fixed effects regression
  model <- feols(Xit ~ 1 | id + time, data = df, weights = ~weight)

  # Return fitted (predicted) values from the model
  return(fitted(model))
}

matrix_to_panel_df <- function(X_mat) {
  if (!is.matrix(X_mat)) {
    stop("Input must be a matrix.")
  }

  # Flatten matrix into vector
  X_vec <- as.vector(X_mat)

  # Get dimensions
  n_rows <- nrow(X_mat)
  n_cols <- ncol(X_mat)

  # Create row (id) and column (time) indices
  id <- rep(1:n_rows, times = n_cols)
  time <- rep(1:n_cols, each = n_rows)

  # Return data frame
  return(data.frame(id = id, time = time, X = X_vec))
}
