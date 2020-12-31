library(spectralGraphTopology)

#' @export
learn_regular_heavytail_graph <- function(X,
                                          heavy_type = "gaussian", nu = NULL,
                                          w0 = "naive", d = 1,
                                          rho = 1, update_rho = TRUE, maxiter = 10000, reltol = 1e-5, verbose = TRUE,
                                          record_objective = FALSE) {
  X <- as.matrix(X)
  # number of nodes
  p <- ncol(X)
  # number of observations
  n <- nrow(X)
  LstarSq <- vector(mode = "list", length = n)
  for (i in 1:n)
    LstarSq[[i]] <- Lstar(X[i, ] %*% t(X[i, ])) / (n-1)
  # w-initialization
  w <- spectralGraphTopology:::w_init(w0, MASS::ginv(cor(X)))
  A0 <- A(w)
  A0 <- A0 / rowSums(A0)
  w <- spectralGraphTopology:::Ainv(A0)
  J <- matrix(1, p, p) / p
  # Theta-initilization
  Lw <- L(w)
  Theta <- Lw
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
    LstarSweighted <- rep(0, .5*p*(p-1))
    if (heavy_type == "student") {
      for (q in 1:n)
        LstarSweighted <- LstarSweighted + LstarSq[[q]] * compute_student_weights(w, LstarSq[[q]], p, nu)
    } else if(heavy_type == "gaussian") {
      for (q in 1:n)
        LstarSweighted <- LstarSweighted + LstarSq[[q]]
    }
    grad <- LstarSweighted - Lstar(rho * Theta + Y) + Dstar(y - rho * d) + rho * (LstarLw + DstarDw)
    eta <- 1 / (2*rho * (2*p - 1))
    wi <- w - eta * grad
    wi[wi < 0] <- 0
    Lwi <- L(wi)
    # update Theta
    eig <- eigen(rho * (Lwi + J) - Y, symmetric = TRUE)
    V <- eig$vectors
    gamma <- eig$values
    Thetai <- V %*% diag((gamma + sqrt(gamma^2 + 4 * rho)) / (2 * rho)) %*% t(V) - J
    # update Y
    R1 <- Thetai - Lwi
    Y <- Y + rho * R1
    # update y
    R2 <- diag(Lwi) - d
    y <- y + rho * R2
    # compute primal, dual residuals, & lagrangian
    primal_lap_residual <- c(primal_lap_residual, norm(R1, "F"))
    primal_deg_residual <- c(primal_deg_residual, norm(R2, "2"))
    dual_residual <- c(dual_residual, rho*norm(Lstar(Theta - Thetai), "2"))
    lagrangian <- c(lagrangian, compute_augmented_lagrangian_ht(wi, LstarSq, Thetai, J, Y, y, d, heavy_type, n, p, rho, nu))
    # update rho
    if (update_rho) {
      s <- rho * norm(Lstar(Theta - Thetai), "2")
      r <- norm(R1, "F")
      if (r > mu * s)
        rho <- rho * tau
      else if (s > mu * r)
        rho <- rho / tau
    }
    if (verbose)
      pb$tick()
    has_converged <- (norm(Lw - Lwi, 'F') / norm(Lw, 'F') < reltol) && (i > 1)
    elapsed_time <- c(elapsed_time, proc.time()[3] - start_time)
    if (has_converged)
      break
    w <- wi
    Lw <- Lwi
    Theta <- Thetai
  }
  results <- list(laplacian = L(wi),
                  adjacency = A(wi),
                  theta = Thetai,
                  maxiter = i,
                  convergence = has_converged,
                  primal_lap_residual = primal_lap_residual,
                  primal_deg_residual = primal_deg_residual,
                  dual_residual = dual_residual,
                  lagrangian = lagrangian,
                  elapsed_time = elapsed_time)
  return(results)
}

compute_student_weights <- function(w, LstarSq, p, nu) {
  return((p + nu) / (sum(w * LstarSq) + nu))
}

compute_augmented_lagrangian_ht <- function(w, LstarSq, Theta, J, Y, y, d, heavy_type, n, p, rho, nu) {
  eig <- eigen(Theta + J, symmetric = TRUE, only.values = TRUE)$values
  Lw <- L(w)
  Dw <- diag(Lw)
  u_func <- 0
  if (heavy_type == "student") {
    for (q in 1:n)
      u_func <- u_func + (p + nu) * log(1 + n * sum(w * LstarSq[[q]]) / nu)
  } else if (heavy_type == "gaussian"){
    for (q in 1:n)
      u_func <- u_func + sum(n * w * LstarSq[[q]])
  }
  u_func <- u_func / n
  return(u_func - sum(log(eig)) + sum(y * (Dw - d)) + sum(diag(Y %*% (Theta - Lw)))
         + .5 * rho * (norm(Dw - d, "2")^2 + norm(Lw - Theta, "F")^2))
}
