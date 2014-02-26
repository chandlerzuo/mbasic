\name{MBASIC.sim}
\alias{MBASIC.sim}
\title{Simulate data for the general MBASIC model.}
\usage{
MBASIC.sim(xi, family = "lognormal", struct = NULL, I, fac, J, S = 2,
  f = 5, delta = NULL, delta.non = NULL, zeta = 0.1)
}
\arguments{
  \item{xi}{Parameter for the magnitude of each
  observations. See details for more information.}

  \item{family}{A parameter for the family of distribution
  to be used. Either "lognormal" or "negbin". Default:
  "lognormal".}

  \item{f}{A numerical value that determine the difference
  of the means between different states. See details for
  more information. Default: 5.}

  \item{I}{An integer for the total number of units.}

  \item{fac}{A vector of length N denoting the experimental
  condition for each replicate.}

  \item{J}{An integer for the number of clusters.}

  \item{S}{An integer for the number of states. Default:
  2.}

  \item{struct}{An K by J integer matrix. The j-th column
  denotes the levels for the cluster level parameter. See
  details for more information. Default: NULL.}

  \item{delta}{A vector of length S, or NULL. This is the
  dirichlet prior parameter to simulate the probability
  across the S states for each CLUSTERED unit and each
  experiment. If NULL, rep(0.1,S) is used.}

  \item{delta.non}{A vector of length S, or NULL. This is
  the dirichlet prior parameter to simulate the probability
  across the S states for each UNCLUSTERED unit and each
  experiment. If NULL, rep(0.1,S) is used.}

  \item{zeta}{The probability that each unit does not
  belong to any cluster. Default: 0.1.}
}
\value{
A list containing: \tabular{ll}{ Y \tab A N by I matrix.
The (n,i)-th entry is the observed value at the i-th unit
for the n-th experiment. \cr Theta \tab A K by I matrix.
The (k,i)-th element is the indicator of whether the i-th
unit is binding for condition k.\cr W \tab A K by (J*S)
matrix. The (k,J*(s-1)+j)-th entry is the probability of
units in the j-th cluster have state s under the k-th
experimental condition. \cr Z \tab An I by J matrix. Each
column is the indicator for an individual loci set.\cr
delta \tab Same as the input argument 'delta'.\cr zeta \tab
Same as the input argumnent 'zeta'.\cr non.id \tab A vector
of length I indicating the number of each unit not to
belong to any cluster.\cr prior.mean \tab A vector of
length S, the hyper means for the S states for each
experiment.\cr prior.sd \tab A vector of length S, the
hyper parameters for the dispersions for the S states for
each experiment.\cr Mu \tab A N by S matrix. The (n,s)-th
entry is the mean values of the response for the s-th state
in the n-th experiment.\cr stdev \tab A vector of length S.
The dispersion parameter for the S states for the observed
data.\cr }
}
\description{
Simulate data for the general MBASIC model.
}
\details{
MBASIC.sim allows two types of distributions. For the
"lognormal" family, entries in the matrix Y follows
distribution: log( Y[n,i] + 1 ) | Theta[n,i]=s ~ N(
Mu[n,s], stdev[s] ). For the "negbin" family, entris in the
matrix Y follows distribution: Y[n,i] | Theta[n,i]=s ~ NB(
Mu[n,s], stdev[s] ). In this package, NB(mu,size) denotes a
Negative binomial distribution with mean mu and variance
mu(1+mu/size). For both distributions,
Mu[n,s]~N(prior.mean[s],prior.sd[s]). Hyper paramters
prior.mean and prior.sd are set differently under the two
distributional families. For the "lognormal" family, where
prior.mean[s] = xi+log((s-1)(f-1)+1), and
prior.sd=log(f)/30. For the "negbin" family,
prior.mean[s]=xi*((s-1)(f-1)+1), and prior.sd=(f-1)*xi/6.
In general, xi is the mean for the state S=1, and f is
roughly the ratio between the means from state S=2 and S=1.
}
\examples{
dat.sim <- MBASIC.sim( xi = 2, I = 100, fac = rep( 1:5, each = 2 ), J = 3 )
}
\author{
Chandler Zuo \email{zuo@stat.wisc.edu}
}
\seealso{
\code{\link{MBASIC.sim.state}}
}

