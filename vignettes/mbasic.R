## ----style-knitr, eval=TRUE, echo=FALSE, results="asis"---------------------------------
  BiocStyle::latex()

## ----eval=TRUE, echo=FALSE, results="hide"------
rm(list = ls())
mydir <- getwd()
options(width = 50)
setwd("~/mbasic_git/test/")

## ----include=TRUE,eval=FALSE,echo=TRUE,results="markup"----
#  library(MBASIC)

## ----eval=TRUE,echo=FALSE,results="hide"--------
library(MBASIC)

## ----eval=TRUE,echo=TRUE,results="markup"-------
target <- generateSyntheticData(dir = "syntheticData")
head(target)
system("ls syntheticData/*/*", intern = TRUE)[1:5]

## ----eval=TRUE,echo=TRUE,results="markup",tidy=TRUE----
tbl <- ChIPInputMatch(dir = paste("syntheticData/", c("chip", "input"), sep = ""), celltypes = c("Cell1", "Cell2"), suffices = c(".bam", ".bed"), depth = 5)
head(tbl)

## ----eval=TRUE,echo=TRUE,results="markup",cache=FALSE----
conds <- paste(tbl$cell, tbl$factor, sep = ".")

## ----eval=TRUE,echo=FALSE,results="hide",tidy=TRUE----
## remove file 'data.Rda' since it will be generated
## try(file.remove("data.Rda"))
fit <- MBASIC.pipeline(chipfile = tbl$chipfile, inputfile = tbl$inputfile, input.suffix = ".bed", target = target, chipformat = tbl$chipformat, inputformat = tbl$inputformat, fragLen = 150, pairedEnd = FALSE, unique = TRUE, fac = conds, J = 3, S = 2, family = "negbin", datafile = "data.Rda")
class(fit)

## ----eval=FALSE,echo=TRUE,results="hide",tidy=TRUE----
#  ## remove file 'data.Rda' since it will be generated
#  ## try(file.remove("data.Rda"))
#  fit <- MBASIC.pipeline(chipfile = tbl$chipfile, inputfile = tbl$inputfile, input.suffix = ".bed", target = target, chipformat = tbl$chipformat, inputformat = tbl$inputformat, fragLen = 150, pairedEnd = FALSE, unique = TRUE, fac = conds, J = 3, S = 2, family = "negbin", datafile = "data.Rda")
#  class(fit)

## ----eval=TRUE,echo=FALSE,results="hide",tidy=TRUE----
allfits <- MBASIC.pipeline(chipfile = tbl$chipfile, inputfile = tbl$inputfile, input.suffix = ".bed", target = target, chipformat = tbl$chipformat, inputformat = tbl$inputformat, fragLen = 150, pairedEnd = FALSE, unique = TRUE, fac = conds, J = 3:10, S = 2, family = "negbin", datafile = "data.Rda")
names(allfits)
class(allfits$BestFit)

## ----eval=FALSE,echo=TRUE,results="hide",tidy=TRUE----
#  allfits <- MBASIC.pipeline(chipfile = tbl$chipfile, inputfile = tbl$inputfile, input.suffix = ".bed", target = target, chipformat = tbl$chipformat, inputformat = tbl$inputformat, fragLen = 150, pairedEnd = FALSE, unique = TRUE, fac = conds, J = 3:10, S = 2, family = "negbin", datafile = "data.Rda")
#  names(allfits)
#  class(allfits$BestFit)

## ----eval=FALSE, echo=TRUE, results="hide",tidy=TRUE----
#  ## Step 2: Generate mapped count matrices
#  dat <- generateReadMatrices(chipfile = tbl$chipfile, inputfile = tbl$inputfile, input.suffix = ".bed", target = target, chipformat = tbl$chipformat, inputformat = tbl$inputformat, fragLen = 150, pairedEnd = FALSE, unique = TRUE)

## ----eval=FALSE, echo=TRUE, results="hide"------
#  conds <- paste(tbl$cell, tbl$factor, sep = ".")

## ----eval=TRUE, echo=FALSE, results="hide"------
## Step 2: Generate mapped count matrices
dat <- generateReadMatrices(chipfile = tbl$chipfile, inputfile = tbl$inputfile, input.suffix = ".bed", target = target, chipformat = tbl$chipformat, inputformat = tbl$inputformat, fragLen = 150, pairedEnd = FALSE, unique = TRUE)
conds <- paste(tbl$cell, tbl$factor, sep = ".")
target <- averageMGC(target = target, m.prefix = "syntheticData/mgc/", m.suffix = "_M.txt", gc.prefix = "syntheticData/mgc/", gc.suffix = "_GC.txt")
dat$input1 <- bkng_mean(inputdat = dat$input, target = target, family = "negbin")

## ----eval=FALSE,echo=FALSE----------------------
#  save(dat, conds, file = "dat.Rda")

## ----eval=FALSE,echo=FALSE----------------------
#  setwd("~/mbasic_git/test/")
#  library(MBASIC)
#  load("dat.Rda")

## ----eval=TRUE, echo=TRUE, results="hide", tidy=TRUE----
## Step 3: Fit an MBASIC model
fit <- MBASIC(Y = t(dat$chip), Gamma = t(dat$input), S = 2, fac = conds, J=3, maxitr = 10, family="negbin")

## ----eval=TRUE, echo=TRUE, results="hide", tidy=TRUE----
## Step 3: Fit an MBASIC model
allfits <- MBASIC.full(Y = t(dat$chip), Gamma = t(dat$input), S = 2, fac = conds, J=3:10, maxitr = 10, family="negbin", ncores = 10)

## ----eval=FALSE,echo=TRUE,results="markup"------
#  showClass("MBASICFit")

## ----eval=TRUE,echo=TRUE,results="markup"-------
dim(fit@Theta)
rownames(fit@Theta)
head(fit@Theta[1, ])

## ----include=TRUE,eval=TRUE,echo=TRUE,fig.align="center",dpi=600,fig.width=6,fig.height=6,tidy=TRUE----
library(gplots)
heatmap.2(fit@Theta[seq(10) + 10, ], Rowv = FALSE, Colv= FALSE, dendrogram = "none", trace = "none")

## ----eval=TRUE, echo=TRUE, results="markup"-----
dim(fit@clustProb)
round(head(fit@clustProb),3)

## ----eval=TRUE, echo=TRUE, results="markup"-----
rownames(fit@W)
dim(fit@W)
round(head(fit@W[seq(10) + 10, ]),3)

## ----eval=TRUE, echo=TRUE, results="hide", tidy=TRUE----
fit.update <- MBASIC(Y = t(dat$chip), Gamma = t(dat$input), S = 2, fac = conds, J = 3, maxitr = 10, family="negbin", Mu.init = fit@Mu, Sigma.init = fit@Sigma, V.init = fit@V, ProbMat.init = fit@Theta, W.init = fit@W, Z.init = fit@Z, b.init = fit@b, P.init = fit@P)

## ----eval=TRUE, echo=TRUE, results="hide", tidy=TRUE----
fit.mix <- MBASIC(Y = t(dat$chip), Gamma = t(dat$input), S = 2, fac = conds, J = 3, maxitr = 10, family = "negbin", statemap = c(1, 2, 2))

## ----eval=TRUE, echo=TRUE, results="markup", tidy=TRUE----
## Simulate data across I=1000 units with J=3 clusters
## There are S=3 states
dat.sim <- MBASIC.sim(xi = 2, family = "lognormal", I = 1000, fac = rep(1:10, each = 2), J = 3, S = 3, zeta = 0.1)

## ----eval=TRUE, echo=TRUE, results="markup"-----
names(dat.sim)
dim(dat.sim$Y)
dim(dat.sim$W)
dim(dat.sim$Theta)

## ----eval=TRUE, echo=TRUE, results="hide", tidy=TRUE----
dat.sim.fit <- MBASIC(Y = dat.sim$Y, S = 3, fac = dat.sim$fac, J = 3, maxitr = 3, para = dat.sim, family = "lognormal")

## ----eval=TRUE,echo=TRUE,results="markup"-------
dat.sim.fit@ARI
dat.sim.fit@W.err
dat.sim.fit@Theta.err
dat.sim.fit@MisClassRate

## ----eval=FALSE,echo=TRUE,results="hide",tidy=TRUE----
#  state.sim <- MBASIC.sim.state(I = 1000, K = 10, J = 4, S = 3, zeta = 0.1)

## ----eval=FALSE,echo=TRUE,results="hide",tidy=TRUE----
#  state.sim.fit <- MBASIC.state(Theta = state.sim$Theta, J = 4, zeta = 0.1)

## ----12,eval=TRUE,echo=FALSE,results="markup",cache=FALSE----
print(sessionInfo())

## ----final,eval=TRUE,echo=FALSE,results="hide"----
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

