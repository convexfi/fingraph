library(spectralGraphTopology)

#' @title Laplacian matrix of a connected graph with Gaussian data
#'
#' Computes the Laplacian matrix of a graph on the basis of an observed data matrix,
#' where we assume the data to be Gaussian distributed.
#'
#' @param S a p x p covariance matrix, where p is the number of nodes in the graph
#' @param w0 initial vector of graph weights. Either a vector of length p(p-1)/2 or
#'        a string indicating the method to compute an initial value.
#' @param d the nodes' degrees. Either a vector or a single value.
#' @param rho ADMM hyperparameter.
#' @param maxiter maximum number of iterations.
#' @param reltol relative tolerance as a convergence criteria.
#' @param verbose whether or not to show a progress bar during the iterations.
#' @export
#' @import spectralGraphTopology
 learn_connected_graph <- function(S,
                                   w0 = "naive",
                                   d = 1,
                                   rho = 1,
                                   maxiter = 10000,
                                   reltol = 1e-5,
                                   verbose = TRUE) {
  # number of nodes
  p <- nrow(S)
  # w-initialization
  w <- spectralGraphTopology:::w_init(w0, MASS::ginv(S))
  A0 <- A(w)
  A0 <- A0 / rowSums(A0)
  w <- spectralGraphTopology:::Ainv(A0)
  LstarS <- Lstar(S)
  J <- matrix(1, p, p) / p
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
    grad <- LstarS - Lstar(Y + rho * Theta) + Dstar(y - rho * d) + rho * (LstarLw + DstarDw)
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
    # update rho
    s <- rho * norm(Lstar(Theta - Thetai), "2")
    r <- norm(R1, "F")
    if (r > mu * s)
      rho <- rho * tau
    else if (s > mu * r)
      rho <- rho / tau
    if (verbose)
      pb$tick()
    has_converged <- (norm(Lwi - Lw, 'F') / norm(Lw, 'F') < reltol) && (i > 1)
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
                  convergence = has_converged)
  return(results)
}
