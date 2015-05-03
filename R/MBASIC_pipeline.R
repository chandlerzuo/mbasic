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
#' @param datafile The file location to save or load the data matrix. See details.
#' @param ... Parameters for function MBASIC.
#' @details
#' This function executes three steps:\cr
#' The first step uses the "generateReadMatrices" function to get the ChIP and Input counts for each locus.\cr
#' The second step is to compute the covariate matrix. If any of 'm.prefix', 'm.suffix', 'gc.prefix', 'gc.suffix' is NULL, then the input count matrix is directly used as the covariate matrix for MBASIC. Alternatively, it will use the 'bkng_mean' to normalize the input count data according to the mappability and GC scores to produce the covariate matrix.\cr
#' The final step is to call the MBASIC function for model fitting.\cr
#' Because the first two steps are time consuming, we recommend in specifying a file location for 'datafile'. Then, when this function executes, it first checks whether 'datafile' exists. If it exists, it will be loaded and the function will jump to the final step. If it does not exist, after the function executes the first two steps, the ChIP data matrix and the covariate matrix will be saved to this file, so that when you rerun this function you do not need to repeat the first two steps.\cr
#' @return A 'MBASICFit' class object.
#' @author Chandler Zuo \email{zuo@@stat.wisc.edu}
#' @examples
#' \dontrun{
#' ## This is the example in our vignette
#' target <- generateSyntheticData(dir = "syntheticData")
#' tbl <- ChIPInputMatch(dir = paste("syntheticData/", c("chip", "input"), sep = ""),suffix = ".bed", depth = 5)
#' conds <- paste(tbl$cell, tbl$factor, sep = ".")
#' MBASIC.fit <- MBASIC.pipeline(chipfile = tbl$chipfile, inputfile = tbl$inputfile, input.suffix = ".bed", target = target, format = "BED", fragLen = 150, pairedEnd = FALSE, unique = TRUE, fac = conds, struct = NULL, S = 2, J = 3, family = "lognormal", maxitr = 10, statemap = NULL)
#'}
#' @export
MBASIC.pipeline <- function(chipfile, inputfile, input.suffix, target, chipformat, inputformat, fragLen, pairedEnd, unique, m.prefix = NULL, m.suffix = NULL, gc.prefix = NULL, gc.suffix = NULL, datafile = NULL, ncores = 10, J, ...) {

  if(is.null(datafile) | (!is.null(datafile) & !file.exists(datafile))) {
    dat <- generateReadMatrices(chipfile = chipfile, inputfile = inputfile, input.suffix = input.suffix, target = target, chipformat = chipformat, inputformat = inputformat, fragLen = fragLen, pairedEnd = pairedEnd, unique = unique)
    if(!is.null(m.prefix) & !is.null(m.suffix) & !is.null(gc.prefix) & !is.null(gc.suffix)) {
      target <- averageMGC(target = target, m.prefix = m.prefix, m.suffix = m.suffix, gc.prefix = gc.prefix, gc.suffix = gc.suffix)
      Gamma <- bkng_mean(inputdat = dat$input, target = target, family = family)
    } else {
      Gamma <- t(dat$input)
    }
    if(!is.null(datafile)) {
      save(dat, Gamma, file = datafile)
    }
  } else {
    load(datafile)
  }

  if(length(J) == 1) {
    return(MBASIC(Y = t(dat$chip), Gamma = Gamma, J = J, ...))
  } else {
    return(MBASIC.full(Y = t(dat$chip), Gamma = Gamma, J = J, ...))
  }
}
