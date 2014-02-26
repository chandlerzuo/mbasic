\name{generateReadMatrices}
\alias{generateReadMatrices}
\title{Map the reads for each pair of matched ChIP and input files to a specified set of genomic intervals.}
\usage{
generateReadMatrices(chipfile, inputfile, input.suffix, target,
  format = "BAM", fragLen = 150, pairedEnd = FALSE, unique = TRUE)
}
\arguments{
  \item{chipfile}{A string vector for the ChIP files.}

  \item{inputfile}{A string vector for the matching input
  files. The length must be the same as 'chipfile'.}

  \item{input.suffix}{A string for the suffix of input
  files.}

  \item{target}{A GenomicRanges object for the target
  intervals where the reads are mapped.}

  \item{format}{A string specifying the type of the file.
  Currently two file types are allowed: "BAM" or "BED".
  Default: "BAM".}

  \item{fragLen}{Either a single value or a 2-column matrix
  of the fragment lengths for the chip and input files.
  Default: 150.}

  \item{pairedEnd}{Either a boolean value or a 2-column
  boolean matrix for whether each file is a paired-end data
  set. Currently this function only allows "BAM" files for
  paired-end data. Default: FALSE.}

  \item{unique}{A boolean value for whether only reads with
  distinct genomic coordinates or strands are mapped.
  Default: TRUE.}
}
\value{
A list of two matrices: \tabular{ll}{ chip \tab A matrix
for the number of mapped reads at each target interval
(row) from each ChIP file (column).\cr input \tab A matrix
for the number of mapped reads at each target interval(
row) from matching input files for the ChIP file
(column).\cr target \tab A GRanges object with sorted
elements.\cr }
}
\description{
Map the reads for each pair of matched ChIP and input files
to a specified set of genomic intervals.
}
\details{
This function uses the readGAlignments and
readGAlignmentsPaired from the \link{GenomicRanges} package
to read BAM files. It uses \link{scan} function to read the
"BED" formatted data, assuming that chr, start, end, strand
information are in column 1, 2, 3, 6.\cr For the input
files, read counts from all files with the same prefix are
added.
}
\author{
Chandler Zuo \email{zuo@stat.wisc.edu}
}

