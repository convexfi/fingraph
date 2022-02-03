library(spectralGraphTopology)

#' @export
#' @import spectralGraphTopology
learn_regular_kcomp_kyfan_graph <- function(S, k = 1, w0 = "naive", d = 1,
                                            rho = 100, update_rho = FALSE, beta = 1e-8, update_beta = TRUE,
                                            maxiter = 10000, reltol = 1e-5, early_stopping = FALSE, verbose = TRUE) {
  # number of nodes
  p <- nrow(S)
  # w-initialization
  w <- spectralGraphTopology:::w_init(w0, MASS::ginv(S))
  A0 <- A(w)
  A0 <- A0 / rowSums(A0)
  w <- spectralGraphTopology:::Ainv(A0)
  # Theta-initilization
  Lw <- L(w)
  Theta <- Lw
  LstarS <- Lstar(S)
  U <- eigen(Lw, symmetric = TRUE)$vectors[, (p - k + 1):p]
  Y <- matrix(0, p, p)
  y <- rep(0, p)
  # ADMM constants
  mu <- 2
  tau <- 2
  # residual vectors
  primal_lap_residual <- c()
  primal_deg_residual <- c()
  dual_residual <- c()
  # augmented lagrangian vector
  lagrangian <- c()
  if (verbose)
    pb <- progress::progress_bar$new(format = "<:bar> :current/:total  eta: :eta",
                                     total = maxiter, clear = FALSE, width = 80)
  elapsed_time <- c()
  start_time <- proc.time()[3]
  for (i in 1:maxiter) {
    # update w
    LstarLw <- Lstar(Lw)
    DstarDw <- Dstar(diag(Lw))
    grad <- LstarS + Lstar(beta * crossprod(t(U)) - Y - rho * Theta) + Dstar(y - rho * d) + rho * (LstarLw + DstarDw)
    eta <- 1 / (2*rho * (2*p - 1))
    wi <- w - eta * grad
    wi[wi < 0] <- 0
    Lwi <- L(wi)
    U <- eigen(Lwi, symmetric = TRUE)$vectors[, (p - k + 1):p]
    # update Theta
    eig <- eigen(rho * Lwi - Y, symmetric = TRUE)
    V <- eig$vectors[,1:(p-k)]
    gamma <- eig$values[1:(p-k)]
    Thetai <- V %*% diag((gamma + sqrt(gamma^2 + 4 * rho)) / (2 * rho)) %*% t(V)
    # update Y
    R1 <- Thetai - Lwi
    Y <- Y + rho * R1
    # update y
    R2 <- diag(Lwi) - d
    y <- y + rho * R2
    primal_lap_residual <- c(primal_lap_residual, norm(R1, "F"))
    primal_deg_residual <- c(primal_deg_residual, norm(R2, "2"))
    dual_residual <- c(dual_residual, rho*norm(Lstar(Theta - Thetai), "2"))
    lagrangian <- c(lagrangian, compute_augmented_lagrangian_kcomp(wi, S, Thetai, U, Y, y, d, rho, beta, p, k))
    # update rho
    if (update_rho) {
      s <- rho * norm(Lstar(Theta - Thetai), "2")
      r <- norm(R1, "F") #+ norm(R2, "2")
      if (r > mu * s)
        rho <- rho * tau
      else if (s > mu * r)
        rho <- rho / tau
    }
    if (update_beta) {
      eig_vals <- spectralGraphTopology:::eigval_sym(L(wi))
      n_zero_eigenvalues <- sum(eig_vals < 1e-9)
      if (k < n_zero_eigenvalues)
        beta <- .5 * beta
      else if (k > n_zero_eigenvalues)
        beta <- 2 * beta
      else {
        if (early_stopping) {
          has_converged <- TRUE
          break
        }
      }
    }
    if (verbose)
      pb$tick()
    has_converged <- (norm(Lwi - Lw, 'F') / norm(Lw, 'F') < reltol) && (i > 1)
    elapsed_time <- c(elapsed_time, proc.time()[3] - start_time)
    if (has_converged)
      break
    w <- wi
    Lw <- Lwi
    Theta <- Thetai
  }
  results <- list(laplacian = L(wi),
                  adjacency = A(wi),
                  maxiter = i,
                  convergence = has_converged,
                  primal_lap_residual = primal_lap_residual,
                  primal_deg_residual = primal_deg_residual,
                  dual_residual = dual_residual,
                  lagrangian = lagrangian,
                  elapsed_time = elapsed_time)
  return(results)
}

compute_augmented_lagrangian_kcomp <- function(w, S, Theta, U, Y, y, d, rho, beta, p, k) {
    eig <- eigen(Theta, symmetric = TRUE, only.values = TRUE)$values[1:(p-k)]
    Lw <- L(w)
    Dw <- diag(Lw)
    return(sum(w * Lstar(S + beta * crossprod(t(U))))
           - sum(log(eig)) + sum(y * (Dw - d)) + sum(diag(Y %*% (Theta - Lw)))
           + .5 * rho * (norm(Dw - d, "2")^2 + norm(Lw - Theta, "F")^2))
}
