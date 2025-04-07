## NetBootstrap

This is the first version of the R package 'NetBootstrap'. 

To install and use the package, we recommend installing the package by
```{r }
install.packages("devtools")  # If not already installed
library(devtools)
install_github("Wei-M-Wei/NetBootstrap")
```
Once installed, load the package with
```{r }
library(NetBootstrap)
```

## Features
- **Main functionality**: Perform parametric bootstrap inference in the netwrok formulation models, see the paper [^1].
- **Validation example**: An example 'test_example.R' is included. 'network_bootstrap(y, X, N, bootstrap_time, data, index, link = 'probit')' is the main function.
- ```{r }
  help(network_bootstrap) # check an example provided
  ```

## Example
```{r }
rm(list = ls())
library(NetBootstrap)

#generate the simulated data
N = 30
K = 1
beta = 1
alpha = rep(0, N)
gamma = rep(0, N)
dens_N_up = log(log(N))
dens_N_low = -log(log(N))
dens_T_up = log(log(N))
dens_T_low = -log(log(N))
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
    v_it = rnorm(1, 0, 1 / sqrt(2))
    epsi_it = rnorm(1, 0, 1)
    if ( t %% 2 == 0 & i %% 2 == 0){
      Z[i, t, ] = (1 - 2 * 1*(0)) * (1 - 2 * 1*(0))
    } else if(t %% 2 == 0 & i %% 2 == 1){
      Z[i, t, ] = (1 - 2 * 1*(1)) * (1 - 2 * 1*(0))
    } else if(t %% 2 == 1 & i %% 2 == 0){
      Z[i, t, ] = (1 - 2 * 1*(0)) * (1 - 2 * 1*(1))
    } else{
      Z[i, t, ] = (1 - 2 * 1*(1)) * (1 - 2 * 1*(1))
    }
    Y[i, t] = (beta * Z[i, t, ] + alpha[i] + gamma[t] > epsi_it)

  }
}


# final data
y = c(Y[!row(Y) == col(Y)])
X_in = matrix(0, N * N, K )
index_matrix = matrix(0, N * N, 2)
for (t in seq(N)) {
  for (i in seq(N)) {
    alpha_in = rep(0, N)
    gamma_in = rep(0, N)
    alpha_in[i] = 1
    gamma_in[t] = 1
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

# estimation procedure
data_in = data.frame(y = y, X = X_design, index = index_in)
index_name = colnames(data_in)[(K+2):(K+3)]
fit = network_bootstrap(y, X = X_design, N, bootstrap_time = 599, index = c('index.1', 'index.2'), data = data_in, link = 'probit', beta_NULL = NULL)

# get the results
est_MLE = fit$cof_MLE[1]
est_corrected = fit$cof[1]
sd = fit$sd

```

A CRAN release is coming soon.

## Reference
[^1]: Haoyaun, X., Wei, M., Dhaene, G., Behyum, J. Parametric bootstrap inference in network formation models. 
