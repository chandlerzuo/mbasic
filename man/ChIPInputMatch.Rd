\name{ChIPInputMatch}
\alias{ChIPInputMatch}
\title{Match the ChIP files with input files.}
\usage{
ChIPInputMatch(dir, suffix, depth = 5)
}
\arguments{
  \item{dir}{A string vector of length 2. The first entry
  is the ChIP file directory and the second is the input
  file directory. It is required that ChIP and input files
  must be in two different directories.}

  \item{depth}{The maximum depth under the dir paths that
  the files are stored. Default: 5.}

  \item{suffix}{A string for the suffix of each file name.}
}
\value{
A data.frame object with the following fields:
\tabular{ll}{ chipfile \tab The name of the ChIP file.\cr
inputfile \tab The prefix of the matching input files (i.e.
without replicate identifiers).\cr lab \tab The string for
the lab information.\cr experiment \tab The string for the
experiment type.\cr cellline \tab The string for the
cellline.\cr factor \tab The string for the factor.\cr
control \tab The string for the control experiment.\cr }
}
\description{
Match the ChIP files with input files.
}
\details{
This function matches the ChIP files with the control files
according to their filenames. The file names must follow
the ENCODE consortium convention, i.e. the prefix of the
file name is
wgEncode<lab><experiment_type><cellline><factor><control>Aln[replicate].
Each string in <.> must has one single initial capital
letter followed by lower case letters or digits. \cr
}
\examples{
\dontrun{ ChIPInputMatch( c( "ChIPDataDir/", "inputDataDir/" ), suffix = ".tagAlign", depth = 2 ) }
}
\author{
Chandler Zuo \email{zuo@stat.wisc.edu}
}

