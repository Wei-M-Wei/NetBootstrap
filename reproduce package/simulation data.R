data_generation = function(N, beta, design = 1, model = 'probit'){
  # sparsity of the fixed effects
  K = length(beta)
  alpha = rep(0, N)
  gamma = rep(0, N)
  if (design == 1){
    dens_N_up = log(log(N))
    dens_N_low = -log(log(N))
    dens_T_up = log(log(N))
    dens_T_low = -log(log(N))
  }else if( design == 2){
    dens_N_up = 0
    dens_N_low = -log(log(N))
    dens_T_up = 0
    dens_T_low = -log(log(N))
  }else if( design == 3){
    dens_N_up = 0
    dens_N_low = -log(N)^(0.5)
    dens_T_up = 0
    dens_T_low = -log(N)^(0.5)
  }else{
    dens_N_up = 0
    dens_N_low = -0.5 * log(N)
    dens_T_up = 0
    dens_T_low = -0.5 * log(N)
  }
  for (i in seq(N)){
    alpha[i] = dens_N_low + ((i - 1) / ( N - 1)) * (dens_N_up - dens_N_low)
  }
  for (i in seq(N)){
    gamma[i] = dens_T_low + ((i - 1) / ( N - 1)) * (dens_T_up - dens_T_low)
  }
  Z = array(0, dim = c(N, N, K))
  Y = matrix(0, nrow = N, ncol = N)
  for (t in seq(N)) {
    for (i in seq(N)) {
      for (k in seq(K)){
      epsi_it = rnorm(1, 0, 1)
      if ( t %% 2 == 0 & i %% 2 == 0){
        Z[i, t, k] = (1 - 2 * 1*(0)) * (1 - 2 * 1*(0))
      } else if(t %% 2 == 0 & i %% 2 == 1){
        Z[i, t, k] = (1 - 2 * 1*(1)) * (1 - 2 * 1*(0))
      } else if(t %% 2 == 1 & i %% 2 == 0){
        Z[i, t, k] = (1 - 2 * 1*(0)) * (1 - 2 * 1*(1))
      } else{
        Z[i, t, k] = (1 - 2 * 1*(1)) * (1 - 2 * 1*(1))
      }
      }

      Y[i, t] = (beta %*% Z[i, t, ] + alpha[i] + gamma[t] > epsi_it)

    }
  }

  # remove the index c(i,j) where i = j
  y = c(Y[!row(Y) == col(Y)])
  X_in = matrix(0, N * N, K )
  index_matrix = matrix(0, N * N, 2)
  for (t in seq(N)) {
    for (i in seq(N)) {
      X_in[i + (t - 1) * N,] = c(Z[i,t,])
      index_matrix[i + (t - 1) * N,] = c(i,t)
    }
  }
  drop_index = NULL
  for (i in seq(N)){
    drop_index = cbind(drop_index, i + (i - 1) * N)
  }
  X_design = as.matrix(X_in[-drop_index,])
  index_in = index_matrix[-drop_index,]

  # final data
  data_in = data.frame(y = y, X = X_design, index = index_in)

  return(list(data = data_in, y = y, X = X_design, Y_original = Y, X_original = Z, alpha = alpha, gamma = gamma))
}

# data_generation = function(N, beta, design = 1, model = 'probit'){
#   # sparsity of the fixed effects
#   K = length(beta)
#   alpha = rep(0, N)
#   gamma = rep(0, N)
#   if (design == 1){
#     dens_N_up = log(log(N))
#     dens_N_low = -log(log(N))
#     dens_T_up = log(log(N))
#     dens_T_low = -log(log(N))
#   }else if( design == 2){
#     dens_N_up = 0
#     dens_N_low = -log(log(N))
#     dens_T_up = 0
#     dens_T_low = -log(log(N))
#   }else if( design == 3){
#     dens_N_up = 0
#     dens_N_low = -log(N)^(0.5)
#     dens_T_up = 0
#     dens_T_low = -log(N)^(0.5)
#   }else{
#     dens_N_up = 0
#     dens_N_low = -log(N)
#     dens_T_up = 0
#     dens_T_low = -log(N)
#   }
#   for (i in seq(N)){
#     alpha[i] = dens_N_low + ((i - 1) / ( N - 1)) * (dens_N_up - dens_N_low)
#   }
#   for (i in seq(N)){
#     gamma[i] = dens_T_low + ((i - 1) / ( N - 1)) * (dens_T_up - dens_T_low)
#   }
#   Z = array(0, dim = c(N, N, K))
#   Y = matrix(0, nrow = N, ncol = N)
#   for (k in seq(K)){
#     x = rnorm(N, 0, 1)
#     for (t in seq(N)){
#       for (i in seq(N)) {
#         Z[i, t, k] = x[i]*x[t]
#         epsi_it = rnorm(1, 0, 1)
#         Y[i, t] = (beta %*% Z[i, t, ] + alpha[i] + gamma[t] > epsi_it)
#       }
#     }
#   }
#
#   # remove the index c(i,j) where i = j
#   y = c(Y[!row(Y) == col(Y)])
#   X_in = matrix(0, N * N, K )
#   index_matrix = matrix(0, N * N, 2)
#   for (t in seq(N)) {
#     for (i in seq(N)) {
#       X_in[i + (t - 1) * N,] = c(Z[i,t,])
#       index_matrix[i + (t - 1) * N,] = c(i,t)
#     }
#   }
#   drop_index = NULL
#   for (i in seq(N)){
#     drop_index = cbind(drop_index, i + (i - 1) * N)
#   }
#   X_design = as.matrix(X_in[-drop_index,])
#   index_in = index_matrix[-drop_index,]
#
#   # final data
#   data_in = data.frame(y = y, X = X_design, index = index_in)
#
#   return(list(data = data_in, y = y, X = X_design, Y_original = Y, X_original = Z, alpha = alpha, gamma = gamma))
# }

