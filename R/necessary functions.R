vector_to_matrix <- function(v, N, ind1, ind2) {


  ind1 <- match(ind1, sort(unique(ind1)))
  ind2 <- match(ind2, sort(unique(ind2)))

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
  Phi_z <- pmin(Phi_z, 1 - 1e-10)

  phi_z <- pmax(phi_z, 1e-10)
  phi_z <- pmin(phi_z, 1 - 1e-10)

  one_minus_Phi_z <- 1 - Phi_z

  dd_F_fix = - eta * phi_z
  ddd_F_fix = eta^2 * phi_z - phi_z


  score1 <- (y / Phi_z - (1 - y) / one_minus_Phi_z) * phi_z

  # Second derivative w.r.t. z
  term1 <- y * ((-z * phi_z) / Phi_z - (phi_z^2) / (Phi_z^2))
  term2 <- (1 - y) * ((-z * phi_z) / one_minus_Phi_z + (phi_z^2) / (one_minus_Phi_z^2))
  score2 <- term1 - term2
  score2 <- pmin(score2, 0)

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

  # Run the weighted two-way fixed effects regression
  model <- feols(Xit ~ -1 | id + time, data = df, weights = ~weight)



  # Return fitted (predicted) values from the model
  return(fitted(model))
}

library(MASS) # for ginv()

get_weighted_fe_projection <- function(X, weight, id, time) {
  n <- length(X)
  stopifnot(length(weight) == n, length(id) == n, length(time) == n)

  id_f <- factor(id)
  time_f <- factor(time)

  # Design matrix with full dummy variables (no intercept)
  D_id <- model.matrix(~ id_f - 1)
  D_time <- model.matrix(~ time_f - 1)
  D <- cbind(D_id, D_time)  # n x (N + T)

  # Weight matrix as vector sqrt for WLS transformation
  w_sqrt <- sqrt(weight)

  # Weighted X and D
  X_w <- X * w_sqrt
  D_w <- D * w_sqrt

  # Compute pseudo-inverse (to handle rank-deficiency)
  DtWD_inv <- ginv(crossprod(D_w))

  # Compute projection matrix P = D (DᵀWD)^- Dᵀ W
  # But better to apply projection as vector operation, to save memory
  fitted_vals <- D %*% (DtWD_inv %*% crossprod(D_w, X_w))

  # Center fitted values weighted by weight to match feols normalization
  # fitted_centered <- fitted_vals - weighted.mean(fitted_vals, weight)

  return(as.vector(fitted_vals))
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

  df = data.frame(id = id, time = time, X = X_vec)
  df <- df[df$id != df$time, ]
  # Return data frame
  return(df)
}


compute_B_hat <- function(D, E, B, C, L) {
  N <- nrow(D)
  result <- 0
  for (i in 1:N) {
    sum_l_term <- 0
    for (l in 1:L) {
      coeff <- N / (N - l)
      E_trans <- E[i, -i]
      D_trans <- D[i, -i]
      for (j in (l + 1):(N-1)) {
        d_index <- j - l
        #print(c(j, d_index, E_trans[ j], D_trans[d_index] ))
        if (d_index >= 1 && d_index <= N) {
          sum_l_term <- sum_l_term + coeff * D_trans[ d_index] * E_trans[j]

        }
      }
    }

    sum_B <- sum(B[i, -i])       # exclude j == i
    sum_C <- sum(C[i, -i])       # exclude j == i
    result <- result + (sum_l_term + sum_B) / sum_C
  }

  return(result)
}


compute_expression_array <- function(A, B) {
  N <- dim(A)[1]
  T <- dim(A)[2]
  d <- dim(A)[3]

  total_sum <- matrix(0, nrow = d, ncol = d)

  for (i in 1:N) {
    t_idx <- setdiff(1:T, i)

    # Term 1: sum_{t ≠ i} sum_{τ ≠ i} A[i,t,] %*% t(A[i,τ,])
    for (t in t_idx) {
      A_it <- A[i, t, ]
      for (tau in t_idx) {
        A_itau <- A[i, tau, ]
        total_sum <- total_sum + A_it %*% t(A_itau)
      }
    }

    # Term 2: sum_{t ≠ i} sum_{j ≠ i, j ≠ t} A[i,t,] %*% t(A[j,t,])
    for (t in t_idx) {
      A_it <- A[i, t, ]
      j_idx <- setdiff(1:N, c(i, t))
      for (j in j_idx) {
        A_jt <- A[j, t, ]
        total_sum <- total_sum + A_it %*% t(A_jt)
      }
    }

    # Term 3: sum_{t ≠ i} B[i,t,] %*% t(B[i,t,])
    for (t in t_idx) {
      B_it <- B[i, t, ]
      total_sum <- total_sum + B_it %*% t(B_it)
    }
  }

  return(total_sum)
}


compute_ape_se <- function(model, variables = NULL, K, beta_hat) {
  stopifnot(inherits(model, "glm"))

  link <- model$family$link
  if (!link %in% c("probit", "logit")) stop("Only probit and logit models are supported.")

  X <- model.matrix(model)
  vcov_beta <- vcov(model)
  vcov_beta[is.na(vcov_beta)] <- 0
  eta <- X %*% beta_hat
  N <- nrow(X)

  # Set up link-specific functions
  pdf_fun <- switch(link,
                    "probit" = dnorm,
                    "logit"  = dlogis)
  cdf_fun <- switch(link,
                    "probit" = pnorm,
                    "logit"  = plogis)

  if (is.null(variables)) {
    variables <- colnames(X)
  }

  results <- data.frame(variable = character(), APE = numeric(), SE = numeric(), stringsAsFactors = FALSE)

  for (var in variables) {
    var_index <- which(colnames(X) == var)
    if (length(var_index) == 0) next

    if (is.numeric(X[, var_index]) && length(unique(X[, var_index])) == 2) {
      is_binary = TRUE
    }else{
      is_binary = FALSE
    }

    if (is_binary) {
      # Binary variable: APE via discrete difference
      X1 <- X; X1[, var_index] <- max(X1[, var_index])
      X0 <- X; X0[, var_index] <- min(X0[, var_index])
      eta1 <- X1 %*% beta_hat
      eta0 <- X0 %*% beta_hat
      p1 <- cdf_fun(eta1)
      p0 <- cdf_fun(eta0)
      ape <- mean(p1 - p0)/(max(X1[, var_index]) - min(X0[, var_index]))

      pdf_eta1 <- pdf_fun(eta1)
      pdf_eta0 <- pdf_fun(eta0)
      pdf_eta1_rep <- matrix(pdf_eta1, nrow = N, ncol = ncol(X), byrow = FALSE)
      pdf_eta0_rep <- matrix(pdf_eta0, nrow = N, ncol = ncol(X), byrow = FALSE)

      grad_i <- (pdf_eta1_rep * X1 - pdf_eta0_rep * X0)
      grad <- colMeans(grad_i)

    } else {
      # Continuous variable: APE via average derivative
      pdf_eta <- pdf_fun(eta)
      ape <- mean(pdf_eta * beta_hat[var_index])
      pdf_eta_rep <- matrix(pdf_eta, nrow = N, ncol = ncol(X), byrow = FALSE)
      variable_indicator <- matrix(rep(colnames(X) == var, each = N), nrow = N)
      grad <- colMeans(pdf_eta_rep * X) * beta_hat[var_index] + colMeans(pdf_eta_rep * variable_indicator)
    }

    se <- sqrt(t(grad) %*% vcov_beta %*% grad)
    results <- rbind(results, data.frame(variable = var, APE = ape, SE = se))
  }

  return(results)
}

compute_TN <- function(Y, p) {
  N <- nrow(Y)
  TN <- 0

  for (i in 1:N) {
    for (j in setdiff(1:N, i)) {
      for (k in setdiff(1:N, c(i, j))) {
        # Observed product
        Y_term <- Y[i, j] * Y[i, k] * Y[k, j]
        # Expected product
        p_term <- p[i,j] * p[i,k] * p[k,j]

        # Add to sum
        TN <- TN + (Y_term - p_term)
      }
    }
  }

  # Normalize
  TN <- TN / (N * (N - 1) * (N - 2))

  return(TN)
}

get_mode <- function(x, plot = FALSE) {
  # Remove NA values
  x <- x[!is.na(x)]

  # Kernel density estimate
  d <- density(x)

  # Find mode (peak of the density)
  mode_val <- d$x[which.max(d$y)]

  # Optionally plot
  if (plot) {
    plot(d, main = "Density Estimate with Mode")
    abline(v = mode_val, col = "red", lwd = 2)
    text(mode_val, max(d$y), labels = paste("Mode ≈", round(mode_val, 2)), pos = 4)
  }

  return(mode_val)
}
