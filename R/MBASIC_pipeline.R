#' @name MBASIC.pipeline
#' @title The pipeline for fitting a MBASIC model for sequencing data.
#' @param chipfile A string vector for the ChIP files.
#' @param inputfile A string vector for the matching input files. The length must be the same as 'chipfile'.
#' @param input.suffix A string for the suffix of input files. If NULL, 'inputfile' will be treated as the full names of the input files. Otherwise, all inputfiles with the initial 'inputfile' and this suffix will be merged.
#' @param target A GenomicRanges object for the target intervals where the reads are mapped.
#' @param chipformat A string specifying the type of the ChIP file. Currently two file types are allowed: "BAM" or "BED". Default: "BAM".
#' @param inputformat A string specifying the type of the input files. Currently two file types are allowed: "BAM" or "BED". Default: "BAM".
#' @param fragLen Either a single value or a 2-column matrix of the fragment lengths for the chip and input files.  Default: 150.
#' @param pairedEnd Either a boolean value or a 2-column boolean matrix for whether each file is a paired-end data set. Currently this function only allows "BAM" files for paired-end data. Default: FALSE.
#' @param unique A boolean value for whether only reads with distinct genomic coordinates or strands are mapped. Default: TRUE.
#' @param m.prefix A string for the prefix of the mappability files.
#' @param m.suffix A string for the suffix of the mappability files. See details for more information. Default: NULL.
#' @param gc.prefix A string for the prefix of the GC files.
#' @param gc.suffix A string for the suffix of the GC files. See details for more information. Default: NULL.
#' @param fac A vector for the experimental conditions corresponding to the ChIP files.
#' @param struct A matrix indicating the levels of the signal matrix.
#' @param J The number of clusters to be identified.
#' @param family The distribution of family to be used. Either "lognormal" or "negbin". See details for more information.
#' @param burnin An integer value for the number of iterations in initialization. Default: 20.
#' @param maxitr The maximum number of iterations in the E-M algorithm. Default: 100.
#' @param tol Tolerance for error in checking the E-M algorithm's convergence. Default: 1e-04.
#' @param nsig The number of mixture components for the distribution of the signal state.
#' @param datafile The file location to save the data matrix.
#' @return A 'MBASICFit' class object.
#' @author Chandler Zuo \email{zuo@@stat.wisc.edu}
#' @examples
#' \dontrun{
#' ## This is the example in our vignette
#' target <- generateSyntheticData(dir = "syntheticData")
#' tbl <- ChIPInputMatch(dir = paste("syntheticData/", c("chip", "input"), sep = ""),suffix = ".bed", depth = 5)
#' conds <- paste(tbl$cell, tbl$factor, sep = ".")
#' MBASIC.fit <- MBASIC.pipeline(chipfile = tbl$chipfile, inputfile = tbl$inputfile, input.suffix = ".bed", target = target, format = "BED", fragLen = 150, pairedEnd = FALSE, unique = TRUE, m.prefix = "syntheticData/mgc/", m.suffix = "_M.txt", gc.prefix = "syntheticData/mgc/", gc.suffix = "_GC.txt", fac = conds, struct = NULL, J = 3, family = "negbin", burnin = 20, maxitr = 100, tol = 1e-4, nsig = 2, datafile = NULL)
#'}
#' @export
MBASIC.pipeline <- function(chipfile, inputfile, input.suffix, target, chipformat, inputformat, fragLen, pairedEnd, unique, m.prefix, m.suffix, gc.prefix, gc.suffix, fac, struct, J, family, burnin = 20, maxitr = 100, tol = 1e-4, nsig = 2, datafile) {

  target <- averageMGC(target = target, m.prefix = m.prefix, m.suffix = m.suffix, gc.prefix = gc.prefix, gc.suffix = gc.suffix)

  dat <- generateReadMatrices(chipfile = chipfile, inputfile = inputfile, input.suffix = input.suffix, target = target, chipformat = chipformat, inputformat = inputformat, fragLen = fragLen, pairedEnd = pairedEnd, unique = unique)
  
  Mu0 <- bkng_mean(inputdat = dat$input, target = target, family = family)

  if(!is.null(datafile))
    save(dat, Mu0, file = datafile)
  
  return(MBASIC.binary(t(dat$chip), t(Mu0), fac, J=J, zeta=0.2, maxitr = maxitr, burnin = burnin, outfile=NULL, out=NULL, init.mod = NULL, struct = struct, family=family, tol = tol, nsig = nsig))
  
}
