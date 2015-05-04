## ----eval=TRUE, echo=FALSE,results="asis"-----------------------------------------------
BiocStyle::latex()
knitr::opts_chunk$set(eval = TRUE, message = FALSE, cache = TRUE, echo = TRUE, results = "markup", tidy = TRUE, warning = FALSE)

## ----cache=FALSE------------------------------------------------------------------------
library(MBASIC)

## ----eval=FALSE,echo=FALSE,cache=FALSE--------------------------------------------------
#  save(dat, conds, file = "dat.Rda")

## ----eval=FALSE,echo=FALSE,cache=FALSE--------------------------------------------------
#  setwd("~/mbasic_git/test/")
#  library(MBASIC)
#  load("dat.Rda")

## ----fig.align="center",dpi=600,fig.width=6,fig.height=6,cache=FALSE--------------------
plot(fit, slot = "Theta", xlab = "Locus", state = 2, cexRow = 0.6, cexCol = 0.4)

## ----fig.align="center",dpi=600,fig.width=6,fig.height=6,cache=FALSE--------------------
plot(fit, slot = "W", state = 2, cexRow = 0.6, cexCol = 1, srtCol = 0, adjCol = c(0.5, 1))

## ----fig.align="center",dpi=600,fig.width=4,fig.height=4,cache=FALSE--------------------
plot(fit, slot = "Mu", state = 2) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

## ----fig.align="center",dpi=600,fig.width=4,fig.height=5,cache=FALSE--------------------
## which replicate to plot
repid <- 1
chip.counts <- dat$chip[, repid]
input.counts <- dat$input[, repid]
pred.states <- as.character(apply(fit@Theta[rownames(fit@Theta) == conds[repid], ], 2, which.max))
ggplot() + geom_point(aes(x = input.counts, y = chip.counts, color = pred.states)) +
  theme(legend.position = "top")

