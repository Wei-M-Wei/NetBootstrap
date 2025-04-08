library(NetBootstrap)
library(readstata13)


# Read the Stata file
data <- read.dta13("Log of Gravity.dta")

# initial settings and the data
bootstrap_time = 599
N = 136
y = as.numeric(data$trade != 0)
X = data[, c("ldist", "border", "comlang", "colony","comfrt")]

# estimation and the index names can be found in the data
est = network_bootstrap(y, X, N, bootstrap_time = bootstrap_time, data = data, index= c('s1_im', 's2_ex'))

# results
cof_MLE = est$cof_MLE[1:5]
cof_mean_bootstrap = est$cof_mean[1:5]
sd_MLE = sqrt(diag(est$Hessian_MLE))[1:5]
sd_bootstrap = est$sd
