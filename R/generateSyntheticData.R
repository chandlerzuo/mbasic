#' @name generateSyntheticData
#' @title Generate synthetic BED data, M and GC files.
#' @param dir The directory for the generated datasets. If this directory already exists, it will be removed before being recreated.
#' @param nchr The number of chromosomes. Default: 5.
#' @param K The number of different transcription factors. Default: 5.
#' @param I The number of loci. Default: 100.
#' @param J The number of clusters among the loci. Default: 3.
#' @description
#' This function generates three set of files. ChIP BED files are stored in 'chip/', input BED files are stored in 'input/', mappability and GC files are stored in 'mgc/'. All BED files follow the name convention for the ENCODE consortium.
#' Each chromosome has a size of 10K, and each locus has 20 bp.
#' ChIP data for two celllines, each with K transcription factors, are generated. All ChIP experiments with the same celltype have the same matching input. Each ChIP experiment has randomly 1-3 replicates. Each input has 3 replicates.
#' @return A RangedData object.
#' @author Chandler Zuo \email{zuo@@stat.wisc.edu}
#' @examples \dontrun{generateSyntheticData("tmpData/", nchr = 5, K = 5, I = 100, J = 3) }
#' @export
generateSyntheticData <- function(dir, nchr = 5, K = 5, I = 100, J = 3) {

  if(file.exists(dir))
    unlink(dir, recursive = TRUE)
  dir.create(dir)
  
  dir.create(chipdir <- paste(dir, "/chip/", sep = ""))
  dir.create(inputdir <- paste(dir, "/input/", sep = ""))
  dir.create(mgcdir <- paste(dir, "/mgc/", sep = ""))
  
  n <- sample(1:3, 2 * K, replace = TRUE)
  dat.sim <-  MBASIC.sim(2, family = "negbin", struct = NULL, I = 100, fac = rep(seq_along(n), n), J = J, S = 2, f = 5, zeta = 0.2)

  cellid <- rep(1:2, each = K)
  tfid <- rep(seq_len(K), 2)

  ## simulate the target intervals
  ir.start <- sample(seq_len(I * nchr * 20 / 40), I) * 40
  ir.end <- ir.start + 20
  ir.start[ ir.end > I * 20 * nchr - 10 ] <- I * 20 * nchr - 30
  ir.start[ ir.start < 30 ] <- 30
  ir.end <- ir.start + 20
    
  target <- RangedData(chromosome = sample(paste("chr", seq_len(nchr), sep = ""), length(ir.start), replace = TRUE), IRanges(start = ir.start, end = ir.end))

  Y <- dat.sim$Y

  idrow <- 0
  idfac <- 0
  for(idcell in unique(cellid)) {
    for(idtf in seq_len(K)) {
      idfac <- idfac + 1
      for(idrep in seq_len(n[ idfac ])) {
        idrow <- idrow + 1
        midpos <- rep(start(target), Y[ idrow, ]) + sample(1:20, sum(Y[ idrow, ]), replace = TRUE)
        chr <- rep(target$chromosome, Y[ idrow, ])
        strand <- sample(c("+", "-"), sum(Y[ idrow, ]), replace = TRUE)
        dat <- data.frame(chr = chr, start = midpos - 10, end = midpos + 9, unknown = NA, unknown2 = NA, strand = strand)
        filename <- paste(chipdir, "wgEncodeLabExpCell", idcell, "Fac", idtf, "CtrlAlnRep", idrep, ".bed", sep = "")
        write.table(dat, file = filename, col.names = FALSE, row.names = FALSE, quote = FALSE)
      }
    }
  }

  for(idcell in unique(cellid)) {
    for(idrep in seq_len(3)) {
      inputcount <- sample(1:3, I, replace = TRUE)
      midpos <- rep(start(target), inputcount) + sample(1:20, sum(inputcount), replace = TRUE)
      chr <- rep(target$chromosome, inputcount)
      strand <- sample(c("+", "-"), sum(inputcount), replace = TRUE)
      dat <- data.frame(chr = chr, start = midpos - 10, end = midpos + 9, unknown = NA, unknown2 = NA, strand = strand)
      filename <- paste(inputdir, "wgEncodeLabExpCell", idcell, "InputCtrlAlnRep", idrep, ".bed", sep = "")
      write.table(dat, file = filename, col.names = FALSE, row.names = FALSE, quote = FALSE)
    }
  }

  ## generate M and GC files

  for(i in seq_len(nchr)) {
    write.table(cbind((0 : 99) * 100, sqrt(runif(100))), col.names = FALSE, row.names = FALSE, file = paste(mgcdir, "chr", i, "_M.txt", sep = ""), quote = FALSE)
    write.table(cbind((0 : 99) * 100, sqrt(runif(100))), col.names = FALSE, row.names = FALSE, file = paste(mgcdir, "chr", i, "_GC.txt", sep = ""), quote = FALSE)
  }
  gc()

  return(target)
}
