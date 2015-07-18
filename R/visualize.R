#' @name plot
#' @title Plot the state probability for a fitted MBASIC model.
#' @param x An \linkS4class{MBASICFit} object.
#' @param y Missing.
#' @param ... Other arguments passed to the function \code{\link[gplots]{heatmap.2}}.
#' @return A heatmap showing the probability for a given state.
#' @author Chander Zuo \email{zuo@@stat.wisc.edu}
#' @aliases plot, MBASICFit-methods
#' @docType methods
#' @rdname plot-methods
if (!isGeneric("plot"))
      setGeneric("plot", function(x, y, ...) standardGeneric("plot"))

#' @rdname plot-methods
#' @param state The probability for this state is plotted.
#' @param xlab The x-axis label.
#' @param ylab The y-axis label.
#' @param main The title.
#' @param slot The slot to be plotted, "Theta" (default), "W", "Mu", or "Sigma".
#' @param Rowv,Colv,trace,dendrogram,... Other arguments passed to \code{\link[gplots]{heatmap.2}}.
#' @details This function provides four differerent kinds of plots.\cr
#' If \code{slot='Theta'}, a heatmap for the probability for each unit to have state \code{state} under each condition is drawn.\cr
#' If \code{slot='W'}, a heatmap for the probability for each cluster to have state \code{state} under each condition is drawn.\cr
#' If \code{slot=}'Mu' or 'Sigma', the range of the fitted 'Mu' or 'Sigma' parameters among different replicates for the same condition is drawn.\cr
#' @aliases plot, MBASICFit-methods
#' @importFrom gplots heatmap.2
#' @import ggplot2
#' @export
setMethod("plot",
          signature = c(x = "MBASICFit", y = "missing"),
          definition = function(x, y = numeric(0), state = 1, slot = "Theta", xlab = NULL, ylab = NULL, main = NULL, Rowv= FALSE, Colv = FALSE, trace = "none", dendrogram = "none", ...) {
            if(!slot %in% c("Theta", "W", "Mu", "Sigma")) {
              stop("Plot for slot ", slot, " is not currently supported.")
            }
            if(is.null(rownames(x@Theta))) {
              stop("Row names for the slot Theta is NULL. There is no way to decide the number of states.")
            }
            K <- length(unique(rownames(x@Theta)))
            S <- nrow(x@Theta) / K
            if(S == 1) {
              S <- max(x@Theta)
            }
            ids <- NULL
            J <- ncol(x@clustProb) - 1
            clusterLabels <- apply(x@clustProb, 1, which.max) - 1
            clusterSizes <- rep(0, J)
            for(j in seq(J)) {
              ids <- c(ids, which(clusterLabels == j))
              clusterSizes[j] <- sum(clusterLabels == j)
            }
            if(slot == "Theta") {
              if(is.null(xlab)) {
                xlab <- "Unit"
              }
              if(is.null(ylab)) {
                ylab <- "Condition"
              }
              if(is.null(main)) {
                main <- paste("State", state)
              }
              if(dim(x@Theta) > K) {
                Theta <- x@Theta[seq(K) + K * (state - 1), ]
              } else {
                Theta <- matrix(as.numeric(x@Theta == state), nrow = nrow(x@Theta))
              }
              colnames(Theta) <- seq(ncol(Theta))
              Theta <- Theta[, ids]
              xline.pos <- cumsum(clusterSizes) + 0.5
              xline.pos <<- xline.pos[-length(xline.pos)]
              heatmap.2(Theta, xlab = xlab, ylab = ylab, add.expr = abline(v = xline.pos), Rowv = Rowv, Colv = Colv, trace = trace, dendrogram = dendrogram, main = main, ...)
            } else if(slot == "W") {
              if(is.null(xlab)) {
                xlab <- "Cluster"
              }
              if(is.null(ylab)) {
                ylab <- "Condition"
              }
              if(is.null(main)) {
                main <- paste("State", state)
              }
              W <- x@W[seq(K) + K * (state - 1), ]
              colnames(W) <- paste("Cluster ", seq(J), ", size ", clusterSizes, sep = "")
              heatmap.2(W, xlab = xlab, ylab = ylab, Rowv = Rowv, Colv = Colv, trace = trace, dendrogram = dendrogram, main = main, ...)
            } else if(slot %in% c("Mu", "Sigma")) {
              if(slot == "Mu") {
                dat <- x@Mu
              } else {
                dat <- x@Sigma
              }
              allconds <- rownames(dat)
              uniqconds <- unique(allconds)
              minvals <- maxvals <- rep(0, length(uniqconds))
              names(minvals) <- names(maxvals) <- uniqconds
              for(cond in uniqconds) {
                minvals[cond] <- min(dat[allconds == cond, state])
                maxvals[cond] <- max(dat[allconds == cond, state])
              }
              plotdat <- data.frame(Condition = uniqconds,
                                    min = minvals,
                                    max = maxvals)
              if(is.null(ylab)) {
                ylab <- slot
              }
              if(is.null(xlab)) {
                xlab <- "Condition"
              }
              if(is.null(main)) {
                main <- paste("Component", state)
              }
              ggplot(data = plotdat, aes(x = Condition)) + geom_errorbar(aes(ymax = max, ymin = min)) + xlab(xlab) + ylab(ylab) + ggtitle(main)
            }
          }
          )
