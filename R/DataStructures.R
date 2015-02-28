#'An S-4 class containing the model fit information for MBASIC model.
#'
#'\describe{
#'\item{Theta}{A K by I matrix. The (k,i)-th element is the posterior probability for the i-th locus to be UNENRICHED under the k-th experimental condition.}
#'\item{W}{A K by J matrix. The (k,j)-th element is the posterior probability that loci in the j-th group are UNENRICHED under the k-th experimental condition.}
#'\item{V}{A K by M matrix. The (k,m)-th element is the weight of the m-th component for the k-th experimental condition. For each row, entries corresponding to the same state should sum to 1.}
#'\item{Z}{An I by J matrix. The (i,j)-th element is the posterior probability that the i-th unit (locus) belongs to the j-th cluster.}
#' \item{clustProb}{An I by (J + 1) matrix. The(i, j+1)-th element is the posterior probability that the i-th unit belongs to the j-th cluster, and the (i, 1)-th element is the posterior probability that the i-th unit belongs to the singleton cluster.}
#'\item{b}{A vector of length I. The i-th element is the posterior probability that the i-th unit (locus) does not belong to any cluster.}
#'\item{aic}{The AIC value of the fitted model.}
#'\item{bic}{The BIC value of the fitted model.}
#'\item{aicc}{The AICC value of the fitted model.}
#'\item{lik}{The log-likelihood after the final iteration.}
#'\item{alllik}{The vector for the log-likelihood after each E-M iteration.}
#'\item{zeta}{The hyper probability for each unit (locus) to belong to some cluster.}
#'\item{Mu}{This is an N by M matrix, where the (n,m)-th entry is the mean parameter for the m-th component.}
#'\item{Sigma}{An N by M matrix, where the (n,m)-th entry is the dispersion parameter for the m-th component.}
#'\item{sigma0}{(To be discarded) NULL except if this object is fitted by the 'MBASIC.binary' function. In that case, this is a vector of length N. The n-th entry is the dispersion parameter for the unenrichment component of the n-th experiment.}
#'\item{e}{(To be discarded) NULL except if this object is fitted by the 'MBASIC.binary' function. In that case, this is a vector of length N. The n-th entry is the normalization value for the background of the n-th experiment.}
#'\item{probz}{A vector of length J. The j-th entry is the hyper probability of any locus to belong to the j-th cluster.}
#'\item{P}{A matrix of length I by S. The (i,s)-th entry is the probability for the i-th locus to have state s, conditional on that this locus does not belong to any cluster.}
#'\item{converged}{Whether the final model is converged.}
#'\item{Theta.err}{Mean squared error in the Theta matrix.}
#'\item{Mu.err}{Mean squared error in the Mu parameter.}
#'\item{ARI}{Adjusted Rand Index between the fitted clustering and the true clustering.}
#'\item{W.err}{Mean squared error in the W matrix.}
#'\item{MisClassRate}{Mis-classification rate by comparing the fitted clustering and the true clustering.}
#'\item{Iter}{Number of iterations taken.}
#'\item{Loss}{A list object for different terms of loss function values.}
#'\item{Struct}{A matrix for the structure constraints for each cluster.}
#'}
#'
#'@examples showClass("MBASICFit")
#'@name MBASICFit-class
#'@rdname MBASICFit-class
#'@exportClass MBASICFit
setClass("MBASICFit",
         representation=representation(
           Theta = "matrix",
           W = "matrix",
           V = "matrix",
           Z = "matrix",
           clustProb = "matrix",
           b = "numeric",
           aic = "numeric",
           bic = "numeric",
           aicc = "numeric",
           lik = "numeric",
           alllik = "numeric",
           zeta = "numeric",
           Mu = "matrix",
           Sigma = "matrix",
           sigma0 = "numeric",
           e = "numeric",
           probz = "numeric",
           P = "matrix",
           converged = "logical",
           Theta.err = "numeric",
           ARI = "numeric",
           W.err = "numeric",
           MisClassRate = "numeric",
           Mu.err = "numeric",
           Iter = "numeric",
           AssociationMatrix = "matrix",
           Loss = "list",
           Struct = "matrix"
         )
       )
