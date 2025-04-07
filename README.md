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

A CRAN release is coming soon.

## Reference
[^1]: Haoyaun, X., Wei, M., Dhaene, G., Behyum, J. Parametric bootstrap inference in network formation models. 
