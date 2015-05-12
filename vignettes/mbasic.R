## ----eval=TRUE, echo=FALSE,results="asis"-----------------------------------------------
BiocStyle::latex()
knitr::opts_chunk$set(eval = TRUE, message = FALSE, cache = FALSE, echo = TRUE, results = "markup", tidy = TRUE, warning = FALSE)

## ----echo=FALSE,results="hide"------------------
rm(list = ls())
mydir <- getwd()
options(width = 50)
setwd("~/mbasic_git/test/")

## -----------------------------------------------
library(MBASIC)

## ----results="hide"-----------------------------
## Step 2': calculate the mappability and GC-content scores for each locus
target <- averageMGC(target = target, m.prefix = "syntheticData/mgc/", m.suffix = "_M.txt", gc.prefix = "syntheticData/mgc/", gc.suffix = "_GC.txt")
## Step 2": compute the normalized input counts
dat$input1 <- bkng_mean(inputdat = dat$input, target = target, family = "negbin")

## ----eval=FALSE,echo=FALSE,cache=FALSE----------
#  save(dat, conds, file = "dat.Rda")

## ----eval=FALSE,echo=FALSE,cache=FALSE----------
#  setwd("~/mbasic_git/test/")
#  library(MBASIC)
#  load("dat.Rda")

## ----results="hide"-----------------------------
## Step 3: Fit an MBASIC model
fit <- MBASIC(Y = t(dat$chip), Gamma = t(dat$input), S = 2, fac = conds, J=3, maxitr = 10, family="negbin")

## ----results="hide"-----------------------------
## Step 3: Fit an MBASIC model
allfits <- MBASIC.full(Y = t(dat$chip), Gamma = t(dat$input), S = 2, fac = conds, J=3:10, maxitr = 10, family="negbin", ncores = 10)

## ----eval=FALSE---------------------------------
#  showClass("MBASICFit")

## -----------------------------------------------
dim(fit@Theta)
rownames(fit@Theta)
head(fit@Theta[1, ])

## ----fig.align="center",dpi=600,fig.width=6,fig.height=6,cache=FALSE----
plot(fit, slot = "Theta", xlab = "Locus", state = 2, cexRow = 0.6, cexCol = 0.4)

## -----------------------------------------------
dim(fit@clustProb)
round(head(fit@clustProb),3)
clusterLabels <- apply(fit@clustProb, 1, which.max) - 1
table(clusterLabels)

## -----------------------------------------------
rownames(fit@W)
dim(fit@W)
round(head(fit@W[seq(10) + 10, ]),3)

## ----fig.align="center",dpi=600,fig.width=6,fig.height=6,cache=FALSE----
plot(fit, slot = "W", state = 2, cexRow = 0.6, cexCol = 1, srtCol = 0, adjCol = c(0.5, 1))

## -----------------------------------------------
dim(fit@Mu)
dim(fit@Sigma)

## ----fig.align="center",dpi=600,fig.width=4,fig.height=4,cache=FALSE----
plot(fit, slot = "Mu", state = 2) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

## ----fig.align="center",dpi=600,fig.width=4,fig.height=5,cache=FALSE----
## which replicate to plot
repid <- 1
chip.counts <- dat$chip[, repid]
input.counts <- dat$input[, repid]
pred.states <- as.character(apply(fit@Theta[rownames(fit@Theta) == conds[repid], ], 2, which.max))
ggplot() + geom_point(aes(x = input.counts, y = chip.counts, color = pred.states)) +
  theme(legend.position = "top")

## ----results="hide"-----------------------------
fit.mix <- MBASIC(Y = t(dat$chip), Gamma = t(dat$input), S = 2, fac = conds, J = 3, maxitr = 10, family = "negbin", statemap = c(1, 2, 2))

## -----------------------------------------------
mincount.thresholds <- apply(dat$input, 2, function(x) quantile(0.25)) * apply(dat$chip, 2, mean) / apply(dat$input, 2, mean)
mincount.thresholds[mincount.thresholds < 5] <- 5
fit.threshold1 <- MBASIC(Y = t(dat$chip), Gamma = t(dat$input), S = 2, fac = conds, J = 3, maxitr = 10, family = "negbin", statemap = c(1, 2, 2), min.count = mincount.thresholds)

## ----results="hide"-----------------------------
fit.threshold2 <- MBASIC(Y = t(dat$chip), Gamma = t(dat$input), S = 2, fac = conds, J = 3, maxitr = 10, family = "negbin", statemap = c(1, 2, 2), min.count = 5)

## -----------------------------------------------
fit.threshold3 <- MBASIC(Y = t(dat$chip), Gamma = t(dat$input), S = 2, fac = conds, J = 3, maxitr = 10, family = "negbin", statemap = c(1, 2, 2), min.count = c(0, 5, 10))

## -----------------------------------------------
mincount.mat <- cbind(0, mincount.thresholds, 2 * mincount.thresholds)
fit.threshold4 <- MBASIC(Y = t(dat$chip), Gamma = t(dat$input), S = 2, fac = conds, J = 3, maxitr = 10, family = "negbin", statemap = c(1, 2, 2), min.count = mincount.mat)

## ----results="hide"-----------------------------
fit.update <- MBASIC(Y = t(dat$chip), Gamma = t(dat$input), S = 2, fac = conds, J = 3,
                     maxitr = 10, family="negbin", initial = fit)

## -----------------------------------------------
dat.asb <- MBASIC.sim(xi = 10, family = "binom", I = 1000, fac = rep(seq(10), each = 2), J = 3, S = 3, zeta = 0.1)
## The counts from either the maternal or the paternal allele
dim(dat.asb$Y)
## The total number of counts from both alleles
dim(dat.asb$X)

## -----------------------------------------------
## Using binomial distributions
fit.asb.bin <- MBASIC(Y = dat.asb$Y, Gamma = dat.asb$X, S = 3, fac = dat.asb$fac, J = 3, maxitr = 5, para = dat.asb, family = "binom")
## Using gamma-binomial distributions
fit.asb.gb <- MBASIC(Y = dat.asb$Y, Gamma = dat.asb$X, S = 3, fac = dat.asb$fac, J = 3, maxitr = 5, para = dat.asb, family = "gamma-binom")

## -----------------------------------------------
## Simulate data across I=1000 units with J=3 clusters
## There are S=3 states
dat.sim <- MBASIC.sim(xi = 2, family = "lognormal", I = 1000, fac = rep(1:10, each = 2), J = 3, S = 3, zeta = 0.1)

## -----------------------------------------------
names(dat.sim)
dim(dat.sim$Y)
dim(dat.sim$W)
dim(dat.sim$Theta)

## ----results="hide"-----------------------------
dat.sim.fit <- MBASIC(Y = dat.sim$Y, S = 3, fac = dat.sim$fac, J = 3, maxitr = 3, para = dat.sim, family = "lognormal")

## -----------------------------------------------
dat.sim.fit@ARI
dat.sim.fit@W.err
dat.sim.fit@Theta.err
dat.sim.fit@MisClassRate

## ----results="hide"-----------------------------
state.sim <- MBASIC.sim.state(I = 1000, K = 10, J = 4, S = 3, zeta = 0.1)

## ----results="hide"-----------------------------
state.sim.fit <- MBASIC.state(Theta = state.sim$Theta, J = 4, zeta = 0.1)

## ----echo=FALSE---------------------------------
print(sessionInfo())

## ----echo=FALSE,results="hide"------------------
file.remove("data.Rda")
file.remove("syntheticData")
setwd(mydir)

## ----eval=FALSE, echo=TRUE, results="hide"------
#  ## Fit a MBASIC model with 4 clusters
#  MBASIC.binary(Y = t(dat$chip),  Mu0 = t(Mu0),
#                fac = conds,  J=4,  zeta=0.2,
#                maxitr = 100, burnin = 20,
#                init.mod = fit,
#                struct = NULL, family="negbin",
#                tol = 1e-4,  nsig = 2)
#  
#  ## Fit a MBASIC model with more iterations
#  MBASIC.binary(Y = t(dat$chip),  Mu0 = t(Mu0),
#                fac = conds,  J=3,  zeta=0.2,
#                maxitr = 200, burnin = 20,
#                init.mod = fit,
#                struct = NULL, family="negbin",
#                tol = 1e-4,  nsig = 2)

