library(MBASIC)
library(testthat)

myOptim <- function(nstates, scalar, p, q) {
  S <- length(nstates) / 2
  fn <- function(w) {
    sum(nstates[seq(S)] * (1 - w) ^ p) +
      scalar * sum(nstates[seq(S) + S] * (1 - w) ^ q) +
        sum(nstates[seq(S)]) * sum(w ^ p) -
          sum(nstates[seq(S)] * w ^ p) +
            scalar * sum(nstates[seq(S) + S]) * sum(w ^ q) -
              scalar * sum(nstates[seq(S) + S] * w ^ q)
  }
  f <- function(w) {
    fn(c(w, 1 - sum(w)))
  }
  ui <- matrix(0, nrow = 2 * S, ncol = S - 1)
  for(s in seq(S - 1)) {
    ui[seq(2) + 2 * (s - 1), s] <- c(1, -1)
  }
  ui[2 * S - 1, ] <- 1
  ui[2 * S, ] <- -1
  ci <- rep(c(0, -1), S)
  theta <- apply(matrix(nstates, ncol = 2) + 1, 1, sum)
  theta <- theta / sum(theta)
  theta <- theta[seq(length(theta) - 1)]
  constrOptim(theta = theta, f = f, grad = NULL, ui = ui, ci = ci)
}


test_that("Error: solvant not correct", {
  p <- 2
  q <- 0.5

  nstates <- c(200, 300, 500, 200)
  S <- length(nstates) / 2
  scalar <- S ^ (q-p)
  .Call("RSolveW", nstates, scalar, p, q)
  myOptim(nstates, scalar, p, q)$par

  nstates <- c(1, 200, 300, 10, 500, 200)
  S <- length(nstates) / 2
  scalar <- S ^ (q-p)
  .Call("RSolveW", nstates, scalar, p, q)
  myOptim(nstates, scalar, p, q)$par
  
  nstates <- c(0 ,200, 300, 20, 200, 500, 200, 0)
  S <- length(nstates) / 2
  scalar <- S ^ (q-p)
  .Call("RSolveW", nstates, scalar, p, q)
  myOptim(nstates, scalar, p, q)$par
  
})
