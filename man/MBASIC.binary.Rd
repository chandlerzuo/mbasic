\name{MBASIC.binary}
\alias{MBASIC.binary}
\title{Bayesian clustering model for binary state matrix with prior estimated background means.}
\usage{
MBASIC.binary(Y, Mu0, fac, J = NULL, zeta = 0.2, maxitr = 100,
  burnin = 20, outfile = NULL, out = NULL, init.mod = NULL,
  struct = NULL, family = "lognormal", tol = 1e-04, nsig = 1,
  min.count = 5)
}
\arguments{
  \item{Y}{An N by I matrix containing the data from N
  experiments across I observation units (loci).}

  \item{Mu0}{An N by I matrix for the prior estimated mean
  for the background state, for N experiments across the I
  observation units (loci).}

  \item{fac}{A vector of length N denoting the experimental
  condition for each replicate.}

  \item{struct}{A matrix indicating the levels of the
  signal matrix.}

  \item{J}{The number of clusters to be identified.}

  \item{family}{The distribution of family to be used.
  Either "lognormal" or "negbin". See details for more
  information.}

  \item{burnin}{An integer value for the number of
  iterations in initialization. Default: 20.}

  \item{maxitr}{The maximum number of iterations in the E-M
  algorithm. Default: 100.}

  \item{tol}{Tolerance for error in checking the E-M
  algorithm's convergence. Default: 1e-04.}

  \item{nsig}{The number of mixture components for the
  distribution of the signal state.}

  \item{zeta}{The initial value of the proportion of
  unclustered units. Default: 0.2.}

  \item{out}{The file directory for writing fitting
  information in each E-M iteration. Default: NULL ( no
  information is outputted ).}

  \item{outfile}{The file directory for writing the
  intermediate results every 10 E-M iterations. This can be
  useful when the running time until final convergence is
  long. Default: NULL ( no intermediate result is saved ).}

  \item{init.mod}{A 'MBASICFit' class object taken as the
  initial model to initialize parameters. This object can
  be an intermediate result from a not yet converged E-M
  algorithm, or a fitted model with smaller number of
  clusters. The user must be cautious in providing this
  initial model, since it must be fitted using the same
  data as the argument of this function.}
}
\value{
A 'MBASICFit' class object.
}
\description{
Bayesian clustering model for binary state matrix with
prior estimated background means.
}
\details{
Function MBASIC.binary currently supports two different
distributional families: log-normal and negative binomial.
This should be specified by the 'family' argument.\cr For
the log-normal distributions, log(Y+1) is modeled as normal
distributions. For experiment n, if locus i is unenriched,
distribution for log(Y[n,i]+1) is N( e[n] * Mu0[n,i],
sigma0[n] ). If locus i is enriched and the enrichment
state is s, the distribution of log(Y[n,i]+1) is N(
mu1[n,s], sigma1[n,s] ).\cr For the negative binomial
distributions, the meanings of mu1, sigma1, sigma0 changes.
For experiment n, if locus i is unenriched, distribution of
Y[n,i] is NB( Mu0[n,i] * e[n], sigma0[n] ). Otherwise, if
locus i is enriched and the enrichment component is s, the
distribution of Y[n]-1 is NB( mu1[n,s], sigma[n,s] ). In
this package, NB( mu, a ) denotes the negative-binomial
distribution with mean mu and size a (i.e. the variance is
mu*(1+mu/a) ).
}
\examples{
\dontrun{
## Simulate a data
dat <- MBASIC.sim.binary( I = 100, fac = rep( 1:5, each = 2 ), J = 3, f = 5, family = "lognormal")
## Fit the model
fit <- MBASIC.binary( t( dat$Y ), Mu0 = t( dat$Mu0 ), fac = rep( 1:5, each = 2 ), J=3, struct = NULL, family="lognormal" )
}
}
\author{
Chandler Zuo \email{zuo@stat.wisc.edu}
}

