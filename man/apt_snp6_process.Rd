\name{apt.snp6.process}
\alias{apt.snp6.process}
\title{Normalize GenomeWideSNP.6 arrays using Affymetrix Power Tools.}
\description{
  This function performs the normalization of a CEL file for Affymetrix GenomeWide_SNP.6 microarrays using the binary apt-copynumber-workflow (v1.20.0 2015-08-28) from Affymetrix Power Tools.
}
\usage{
  apt.snp6.process(CEL = NULL, samplename = NULL, out.dir = getwd(),
  temp.files.keep = FALSE, force.OS = NULL, apt.build = "na35.r1")
}
\arguments{
  \item{CEL}{A CEL file.}
  \item{samplename}{Name of the analyzed sample.}
  \item{out.dir}{The path to output results.}
  \item{temp.files.keep}{Keep the (numerous) intermediate files.}
  \item{force.OS}{Specifiy the OS type of the execution machine (see Details).}
  \item{apt.build}{Affymetrix genome build version.}
}
\details{
  \code{force.OS} is \code{NULL} by default, but can take any of these values :
  \itemize{
    \item{linux}
    \item{windows}
    \item{osx}
  }
}
\author{
Bastien Job
}

