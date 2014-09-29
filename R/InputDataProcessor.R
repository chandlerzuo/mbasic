#' @name ChIPInputMatch
#' @title Match the ChIP files with input files.
#' @param dir A string vector of length 2. The first entry is the ChIP file directory and the second is the input file directory. It is required that ChIP and input files must be in two different directories.
#' @param depth The maximum depth under the dir paths that the files are stored. Default: 5.
#' @param suffices A vector of strings for the suffices of each file name.
#' @details
#' This function matches the ChIP files with the control files according to their filenames. The file names must follow the ENCODE consortium convention, i.e. the prefix of the file name is wgEncode<lab><experiment_type><cellline><factor><control>Aln[replicate]. Each string in <.> must has one single initial capital letter followed by lower case letters or digits. \cr
#' @return A data.frame object with the following fields:
#' \tabular{ll}{
#' chipfile \tab The name of the ChIP file.\cr
#' inputfile \tab The prefix of the matching input files (i.e. without replicate identifiers).\cr
#' lab \tab The string for the lab information.\cr
#' experiment \tab The string for the experiment type.\cr
#' cellline \tab The string for the cellline.\cr
#' factor \tab The string for the factor.\cr
#' control \tab The string for the control experiment.\cr
#' chiptype \tab The file format of the ChIP data.\cr
#' inputtype \tab The file format of the input data.\cr
#'}
#' @examples \dontrun{ ChIPInputMatch(c("ChIPDataDir/", "inputDataDir/"), suffices = ".tagAlign", depth = 2) }
#' @author Chandler Zuo \email{zuo@@stat.wisc.edu}
#' @export
ChIPInputMatch <- function(dir, suffices, depth = 5) {
  if(length(dir) != 2) {
    stop("Error: the length of dir must be 2.")
  }
  if(dir[1] == dir[2]) {
    stop("Error: the ChIP and input files must be in two different directories.")
  }
  chipdir <- dir[1]
  inputdir <- dir[2]

  chiplab <- chipexp <- chipcell <- chipfac <- chipctrl <-
    inputlab <- inputexp <- inputcell <- inputfac <- inputctrl <- chiptype <- inputtype <-
      chiplist <- inputlist <- NULL
  for(suffix in suffices) {
    chipfilelist <- system(paste("find ", chipdir, " -maxdepth ", depth, " -name *", suffix, sep = ""), intern = TRUE)
    inputfilelist <- system(paste("find ", inputdir, " -maxdepth ", depth, " -name *", suffix, sep = ""), intern = TRUE)
    
    chippath <- extract(chipfilelist, "", "wgEncode*")
    
    exp1 <- "[A-Z][a-z|0-9]*"
    chiplab <- c(chiplab, extract(chipfilelist, "*wgEncode", exp1))
    chipexp <- c(chipexp, extract(chipfilelist, paste("wgEncode", paste(rep(exp1, 1), collapse = ""), sep = ""), exp1))
    chipcell <- c(chipcell, extract(chipfilelist, paste("wgEncode", paste(rep(exp1, 2), collapse = ""), sep = ""), exp1))
    chipfac <- c(chipfac, extract(chipfilelist, paste("wgEncode", paste(rep(exp1, 3), collapse = ""), sep = ""), exp1))
    chipctrl <- c(chipctrl, extract(chipfilelist, paste("wgEncode", paste(rep(exp1, 4), collapse = ""), sep = ""), exp1))
    chiptype <- c(chiptype, rep(suffix, length(chipfilelist)))
    
    inputlab <- c(inputlab, extract(inputfilelist, "*wgEncode", exp1))
    inputexp <- c(inputexp, extract(inputfilelist, paste("wgEncode", paste(rep(exp1, 1), collapse = ""), sep = ""), exp1))
    inputcell <- c(inputcell, extract(inputfilelist, paste("wgEncode", paste(rep(exp1, 2), collapse = ""), sep = ""), exp1))
    inputfac <- c(inputfac, extract(inputfilelist, paste("wgEncode", paste(rep(exp1, 3), collapse = ""), sep = ""), exp1))
    inputctrl <- c(inputctrl, extract(inputfilelist, paste("wgEncode", paste(rep(exp1, 4), collapse = ""), sep = ""), exp1))
    inputtype <- c(inputtype, rep(suffix, length(inputfilelist)))

    chiplist <- c(chiplist, chipfilelist)
    inputlist <- c(inputlist, inputfilelist)
  }

  ## to extract directories
  strRev <- function(x)
    sapply(strsplit(x, split = ""), function(str) { paste(rev(str), collapse = "") })

  chiptype <- toupper(extract(chiptype, "\\.", "[A-Z|a-z]*"))
  inputtype <- toupper(extract(inputtype, "\\.", "[A-Z|a-z]*"))
  
  chipprefix <- strRev(extract(strRev(chiplist), "w*", "/.*"))
  inputprefix <- strRev(extract(strRev(inputlist), "w*", "/.*"))

  chipconds <- paste(chiplab, chipexp, chipcell, chipctrl, sep = ".")
  inputconds <- paste(inputlab, inputexp, inputcell, inputctrl, sep = ".")
  uniqueconds <- sort(unique(chipconds))
  
  chipfile <- inputfile <- lab <- exper <- cell <- fac <- ctrl <- ctype <- itype <- NULL
  for(cond in uniqueconds) {
    id.chip <- which(chipconds == cond)
    id.input <- which(inputconds == cond)
    if(length(id.input) > 1)
      message("Warning:", cond, "has multiple matching input files")
    inputfilehead <- unique(paste(inputprefix, "wgEncode", inputlab, inputexp, inputcell, inputfac, inputctrl, sep = "")[ id.input])
    if(length(id.chip) > 0) {
      for(i in id.chip) {
        if(length(inputfilehead) > 0) {
          for(i.input in seq_along(id.input)) {
            chipfile <- c(chipfile, chiplist[ i ])
            inputfile <- c(inputfile, inputfilehead[i.input])
            lab <- c(lab, chiplab[ i ])
            exper <- c(exper, chipexp[ i ])
            cell <- c(cell, chipcell[ i ])
            fac <- c(fac, chipfac[ i ])
            ctrl <- c(ctrl, chipctrl[ i ])
            ctype <- c(ctype, chiptype[i])
            itype <- c(itype, inputtype[id.input[i.input]])
          }
        } else {
          chipfile <- c(chipfile, chiplist[ i ])
          inputfile <- c(inputfile, NA)
          lab <- c(lab, chiplab[ i ])
          exper <- c(exper, chipexp[ i ])
          cell <- c(cell, chipcell[ i ])
          fac <- c(fac, chipfac[ i ])
          ctrl <- c(ctrl, chipctrl[ i ])
          ctype <- c(ctype, chiptype[i])
          itype <- c(itype, NA)
        }
      }
    }
  }
    
  return(data.frame(
                     chipfile = chipfile,
                     inputfile = inputfile,
                     lab = lab,
                     experiment = exper,
                     cell = cell,
                     factor = fac,
                     control = ctrl,
                     chipformat = ctype,
                     inputformat = itype
                  )
      )
  
}

#' @name generateReadMatrices
#' @title Map the reads for each pair of matched ChIP and input files to a specified set of genomic intervals.
#' @param chipfile A string vector for the ChIP files.
#' @param inputfile A string vector for the matching input files. The length must be the same as 'chipfile'.
#' @param input.suffix A string for the suffix of input files.
#' @param target A GenomicRanges object for the target intervals where the reads are mapped.
#' @param chipformat A vector of strings specifying the type of the ChIP files. Can be a single value if all ChIP files have the same format. Currently two file types are allowed: "BAM" or "BED". Default: "BAM".
#' @param inputformat A vector of strings specifying the type of the input file. Can be a single string if all input files have the same format. Currently two file types are allowed: "BAM" or "BED". Default: "BAM".
#' @param fragLen Either a single value or a 2-column matrix of the fragment lengths for the chip and input files.  Default: 150.
#' @param pairedEnd Either a boolean value or a 2-column boolean matrix for whether each file is a paired-end data set. Currently this function only allows "BAM" files for paired-end data. Default: FALSE.
#' @param unique A boolean value for whether only reads with distinct genomic coordinates or strands are mapped. Default: TRUE.
#' @param ncores The number of cores. Default: 1.
#' @details
#' This function uses the readGAlignments and readGAlignmentsPaired from the \link{GenomicRanges} package to read BAM files. It uses \link{scan} function to read the "BED" formatted data, assuming that chr, start, end, strand information are in column 1, 2, 3, 6.\cr
#' For the input files, read counts from all files with the same prefix are added.
#' @return A list of two matrices:
#' \tabular{ll}{
#' chip \tab A matrix for the number of mapped reads at each target interval (row) from each ChIP file (column).\cr
#' input \tab A matrix for the number of mapped reads at each target interval(row) from matching input files for the ChIP file (column).\cr
#' target \tab A GRanges object with sorted elements.\cr
#'}
#' @author Chandler Zuo \email{zuo@@stat.wisc.edu}
#' @export
generateReadMatrices <- function(chipfile, inputfile, input.suffix, target, chipformat = "BAM", inputformat = "BAM", fragLen = 150, pairedEnd = FALSE, unique = TRUE, ncores = 1) {
  ## Check the arguments
  require(doMC)
  nfiles <- length(chipfile)
  if(length(inputfile) != nfiles)
    stop("Error: number of matching input files must be the same as the number of ChIP files!")
  if(length(chipformat) == 1)
    chipformat <- rep(chipformat, nfiles)
  if(length(inputformat) == 1)
    inputformat <- rep(inputformat, nfiles)
  if(length(chipformat) != nfiles)
    stop("Error: number of ChIP file formats must match the number of ChIP files!")
  if(length(inputformat) != nfiles)
    stop("Error: number of input file formats must match the number of input files!")
  ## Convert strings to upper cases
  chipformat <- toupper(chipformat)
  inputformat <- toupper(inputformat)
  
  if(prod(c(inputformat, chipformat) %in% c("TAGALIGN", "BAM", "BED", NA)) == 0)
    stop("Error: only BAM, BED and TAGALIGN files are allowed.")
  chipformat[chipformat == "TAGALIGN"] <- "BED"
  inputformat[inputformat == "TAGALIGN"] <- "BED"
  checkMatrixDim <- function(x) {
    if(!is.matrix(x))
      if(length(x) == 1)
        return(TRUE)
    if(is.matrix(x))
      if(nrow(x) == nfiles & ncol(x) == 2)
        return(TRUE)
    return(FALSE)
  }
  if(! checkMatrixDim(fragLen))
    stop(paste("Error: fragLen must be either a single value, or a matrix with", nfiles, "rows and 2 columns."))
  if(! checkMatrixDim(pairedEnd))
    stop(paste("Error: pairedEnd must be either a single value, or a matrix with", nfiles, "rows and 2 columns."))

  if(!is.matrix(fragLen))
    fragLen <- matrix(fragLen, nrow = nfiles, ncol = 2)
  if(!is.matrix(pairedEnd))
    pairedEnd <- matrix(pairedEnd, nrow = nfiles, ncol = 2)

  if(sum(pairedEnd[,1] > 0 & chipformat != "BAM") > 0)
    stop("Error: for paired-end data only BAM format is allowed.")
  if(sum(pairedEnd[,2] > 0 & inputformat != "BAM") > 0)
    stop("Error: for paired-end data only BAM format is allowed.")
  
  target <- sort(target)

  inputfile <- as.character(inputfile)
  chipfile <- as.character(chipfile)

  registerDoMC(ncores)
  
  ## process all input files
  uniqueInputCounts <- foreach(file = na.omit(unique(inputfile))) %dopar% {
      ## For input file, must read all replicates
      if(!is.null(input.suffix))
          listinputstr <- paste("ls ", file, "*", input.suffix, sep = "")
      else
          listinputstr <- paste("ls ", file, sep = "")
      uniqueInputCount <- 0
      for(ifile in system(listinputstr, intern = TRUE)) {
          message(paste("processing input file", ifile))
          rds <- readReads(ifile, extended = TRUE, fragLen = fragLen[ which(inputfile == file)[1], 1 ], pairedEnd = pairedEnd[ which(inputfile == file)[1], 1 ], format = inputformat[which(inputfile == file)[1]])
          if(unique)
              rds <- unique(rds)
          uniqueInputCount <- uniqueInputCount + countOverlaps(target, rds)
      }
      gc()
      uniqueInputCount
  }

  uniqueInputCountsWithoutNA <- NULL
  if(length(uniqueInputCounts) > 0) {
      uniqueInputCountsWithoutNA <- matrix(unlist(uniqueInputCounts), ncol = length(uniqueInputCounts))
  }

  ## process all chip files
  uniqueChIPCounts <- foreach(file = na.omit(unique(chipfile))) %dopar% {
      message(paste("processing chip file", file))
      rds <- readReads(file, extended = TRUE, fragLen = fragLen[ which(chipfile == file)[1], 1 ], pairedEnd = pairedEnd[ which(chipfile == file)[1], 1 ], format = chipformat[which(chipfile == file)[1]])
      if(unique)
          rds <- unique(rds)
      gc()
      countOverlaps(target, rds)
  }

  uniqueChIPCountsWithoutNA <- NULL
  if(length(uniqueChIPCounts) > 0) {
      uniqueChIPCountsWithoutNA <- matrix(unlist(uniqueChIPCounts), ncol = length(uniqueChIPCounts))
  }

  ##message(dim(uniqueChIPCountsWithoutNA))
  ##message(dim(uniqueInputCountsWithoutNA))
  
  uniquechipcounts <- matrix(1, nrow = length(target), ncol = length(unique(chipfile)))
  uniqueinputcounts <- matrix(1, nrow = length(target), ncol = length(unique(inputfile)))
  inputfile[is.na(inputfile)] <- "NA"
  chipfile[is.na(chipfile)] <- "NA"
  colnames(uniqueinputcounts) <- unique(inputfile)
  colnames(uniquechipcounts) <- unique(chipfile)
  if(!is.null(uniqueChIPCountsWithoutNA)) {
      uniquechipcounts[, colnames(uniquechipcounts) != "NA"] <- uniqueChIPCountsWithoutNA
  }
  if(!is.null(uniqueInputCountsWithoutNA)) {
      uniqueinputcounts[, colnames(uniqueinputcounts) != "NA"] <- uniqueInputCountsWithoutNA
  }
  
  allchipcounts <- allinputcounts <- matrix(0, nrow = length(target), ncol = nfiles)

  for(i in seq_len(nfiles)) {
    allchipcounts[ , i ] <- uniquechipcounts[ , chipfile[i] ]
    allinputcounts[ , i ] <- uniqueinputcounts[ , inputfile[i] ]
  }

  gc()
  return(list(chip = allchipcounts, input = allinputcounts, target = target))
  
}

#' @name averageMGC
#' @title Compute the average mappability and GC scores over a set of genomic intervals.
#' @param target A \link{GRanges} object for the target intervals.
#' @param m.prefix A string for the prefix of the mappability files.
#' @param m.suffix A string for the suffix of the mappability files. See details for more information. Default: NULL.
#' @param gc.prefix A string for the prefix of the GC files.
#' @param gc.suffix A string for the suffix of the GC files. See details for more information. Default: NULL.
#' @details
#' If 'm.suffix' is NULL, then a single file with name 'm.prefix' should include mappability scores of all chromosomes, and this file is read. Alternatively, all mappability files with 'm.prefix'<chrXX>'m.suffix' are read. The 'gc.suffix' argument has similar effects.\cr
#' 'target' has to be a sorted GRanges object. If it is not sorted, then the elements are reordered. Users have to make sure that other data sources must follow the same ordering in the elements.
#' @return A GRanges object with two extra fields from the input argument 'target':
#' \tabular{ll}{
#' mappability \tab The average mappability scores for the sorted elements of 'target'.\cr
#' GC \tab The average GC scores for the sorted elements of 'target'.\cr
#' }
#' @author Chandler Zuo \email{zuo@@stat.wisc.edu}
#' @useDynLib MBASIC
#' @export
averageMGC <- function(target, m.prefix, m.suffix = NULL, gc.prefix, gc.suffix = NULL) {
  if(class(target) != "GRanges")
    stop("Error: target must be a GRanges object.")
  if(is.null(m.prefix) & is.null(m.suffix))
    stop("Error: either m.prefix or m.suffix must not be null.")
  if(is.null(gc.prefix) & is.null(gc.suffix))
    stop("Error: either gc.prefix or gc.suffix must not be null.")

  if(is.null(m.suffix) & !file.exists(m.prefix))
    stop(paste("Error: the file", m.prefix, "does not exist."))
  if(is.null(gc.suffix) & !file.exists(gc.prefix))
    stop(paste("Error: the file", gc.prefix, "does not exist."))
  
  target.sorted <- sort(target)
  if(prod(target == target.sorted) == 0)
    message("Warning: target is not sorted. The elements of target will be reordered.")
  target <- target.sorted

  allmap <- allgc <- NULL
  
  if(is.null(m.suffix)) {
    allmap <- read.table(m.prefix)
    names(allmap) <- c("chr", "pos", "score")
    allmap$chr <- as.character(allmap$chr)
  }
  if(is.null(gc.suffix)) {
    allgc <- read.table(gc.prefix)
    names(allgc) <- c("chr", "pos", "score")
    allgc$chr <- as.character(allgc$chr)
  }
  
  allmgc <- NULL
  for(chr in as.character(unique(seqnames(target)))) {
    message(paste("processing", chr))
    if(!is.null(m.suffix))
      map <- as.matrix(read.table(paste(m.prefix, chr, m.suffix, sep = "")))
    else
      map <- as.matrix(allmap[ allmap$chr == chr, -1 ])
    message("read mappability file")
    
    if(!is.null(gc.suffix))
      gc <- as.matrix(read.table(paste(gc.prefix, chr, gc.suffix, sep = "")))
    else
      gc <- as.matrix(allgc[ allgc$chr == chr, -1 ])
    message("read gc file")
    
    genecoord <- cbind(start(ranges(target[ seqnames(target) == chr ])),
                   end(ranges(target[ seqnames(target) == chr ])))
    avgmgc <- .Call("avg_score", genecoord, map[ , 2 ], gc[ , 2 ], diff(map[ 1:2, 1 ]), diff(gc[ 1:2, 1 ]), package = "MBASIC")
    allmgc <- rbind(allmgc, avgmgc)
  }

  elementMetadata(target)[[ "mappability" ]] <- allmgc[ , 1 ]
  elementMetadata(target)[[ "GC" ]] <- allmgc[ , 2 ]
  return(target)
  
}

#' @name bkng_mean
#' @title Compute the background means.
#' @param inputdat A matrix for the input counts at each locus (column) for each experiment (row).
#' @param target A \link{GRanges} object with two fiels named "mappability" and "GC".The length of target must be the same as the column of inputdat.
#' @param nknots A integer for the number of knots for the spline function. Default: 2.
#' @param family A string for the distributional family. Must be either "lognormal" (default) or "negbin".
#' @return A numeric matrix for the background counts at each locus (column) for each experiment (row).
#' @author Chandler Zuo \email{zuo@@stat.wisc.edu}
#' @export
bkng_mean <- function(inputdat, target, nknots = 2, family = "lognormal") {
  options(warn = -1)
  if(class(target) != "GRanges")
    stop("Error: target must be a GRanges object.")
  if(length(target) != nrow(inputdat))
    stop("Error: length of target must be the same as number of columns in inputdat.")
  if(! family %in% c("lognormal", "negbin"))
    stop("Error: family must be either lognormal or negbin.")

  gc_bs <- splines::bs(target$GC, knots = quantile(target$GC, prob = seq(0, 1, length = nknots + 2)[ 2 : (nknots + 1) ]))
  map_bs <- splines::bs(target$mappability, knots = quantile(target$mappability, prob = seq(0, 1, length = nknots + 2)[ 2 : (nknots + 1) ]), degree = 1)

  Mu0 <- inputdat
  for(i in seq_len(ncol(inputdat))) {
    if(var(inputdat[, i]) == 0) {
      Mu0[, i] <- inputdat[1, 1]
      next
    }
    if(family == "lognormal") {
      fit <- lm(log(inputdat[ , i ] + 1) ~ gc_bs + map_bs)
      Mu0[ , i ] <- predict(fit)
    }  else {
      fit <- glm.nb(inputdat[ , i ] ~ gc_bs + map_bs)
      Mu0[ , i ] <- exp(predict.glm(fit))
    }
  }

  return(Mu0)
}

 
#' @name readReads
#' @title Read sequencing reads from either a BAM or a BED file.
#' @param reads The sequencing file.
#' @param extended A boolean value for whether each read will be extended to the fragment length.
#' @param fragLen A numeric value for the fragment length.
#' @param pairedEnd A boolean value for whether the sequencing file is paired end.
#' @param use.names A boolean value to be passed to \link{GenomicRanges} functions.
#' @param format A string of file format. Must be either "BAM" or "BED".
#' @return A \link{GRanges} object.
#' @author Samuel Younkin \email{syounkin@@stat.wisc.edu}
#' @export
readReads <- function(reads, extended, fragLen = 200, pairedEnd = FALSE, use.names = FALSE, format) {
  if(pairedEnd) { # Paired-end reads
    if(format == "BAM") {
      ga.pairs <- readGAlignmentPairs(file = reads, format = "BAM", use.names = use.names)
      reads.gr <- granges(ga.pairs)
      rm(ga.pairs)
      gc()
      return(reads.gr)
    } else {
      stop("Paired-end format must be BAM.")
    }
  } else { # Single-end reads
    if(format == "BED") {
      reads <- scan(file = reads, what = list("",1L,1L,"",1L,""))
      names(reads) <- c("chr","start","end", "unknown", "unknown2", "strand")
      reads.gr <- with(reads, GRanges(seqnames = chr, ranges = IRanges(start=start,end=end), strand = strand))
      rm(reads)
      gc()
    }else if (format == "BAM") {
      ga.single <- readGAlignments(file = reads, format = format, use.names = use.names)
      reads.gr <- granges(ga.single)
      rm(ga.single)
      gc()
    } else {
      stop("Single end format must be BED or BAM.")
    }
    if(extended) {
      reads.forward <- reads.gr[strand(reads.gr)=="+"]
      reads.reverse <- reads.gr[strand(reads.gr)=="-"]
      reads.forward.extended <- GRanges(seqnames=seqnames(reads.forward), ranges = IRanges(start = start(reads.forward), width = fragLen), strand = "+")
      reads.reverse.extended <- GRanges(seqnames=seqnames(reads.reverse), ranges = IRanges(end = end(reads.reverse), width = fragLen), strand = "-")
      reads.gr <- c(reads.forward.extended, reads.reverse.extended)
      names(reads.gr) <- names(c(reads.forward,reads.reverse))
    }
    return(reads.gr)
  }
}
