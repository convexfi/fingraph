library(spectralGraphTopology)

#' @export
learn_kcomp_graph <- function(S, k = 1, w0 = "naive", alpha = 0, eps = 0, d = 1,
                              rho = 100, update_rho = FALSE, maxiter = 10000,
                              reltol = 1e-5, verbose = TRUE) {
  # number of nodes
  p <- nrow(S)
  # w-initialization
  w <- spectralGraphTopology:::w_init(w0, MASS::ginv(S))
  A0 <- A(w)
  A0 <- A0 / rowSums(A0)
  w <- spectralGraphTopology:::Ainv(A0)
  J <- matrix(1, p, p) / p
  H <- alpha * (diag(p) - p * J)
  K <- S
  # Theta-initilization
  Lw <- L(w)
  Theta <- Lw
  Y <- matrix(0, p, p)
  y <- rep(0, p)
  # ADMM constants
  mu <- 2
  tau <- 2
  if (verbose)
    pb <- progress::progress_bar$new(format = "<:bar> :current/:total  eta: :eta",
                                     total = maxiter, clear = FALSE, width = 80)
  for (i in 1:maxiter) {
    # update w
    LstarLw <- Lstar(Lw)
    DstarDw <- Dstar(diag(Lw))
    grad <- Lstar(K - Y - rho * Theta) + Dstar(y - rho * d) + rho * (LstarLw + DstarDw)
    eta <- 1 / (2*rho * (2*p - 1))
    wi <- w - eta * grad
    wi[wi < 0] <- 0
    Lwi <- L(wi)
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
    has_converged <- (norm(Lwi - Lw, 'F') / norm(Lw, 'F') < reltol) && (i > 1)
    if (has_converged)
      break
    w <- wi
    Lw <- Lwi
    Theta <- Thetai
    if (eps > 0 & alpha > 0)
      K <- S + H / (-Lw + eps)
  }
  results <- list(laplacian = L(wi), adjacency = A(wi), maxiter = i,
                  convergence = has_converged)
  return(results)
}
