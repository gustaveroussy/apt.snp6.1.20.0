\name{apt.snp6.process.batch}
\alias{apt.snp6.process.batch}
\title{Run apt.snp6.process() with multithreading}
\description{
  Same as apt.snp6.process() in batch mode, with multithreading.
}
\usage{
  apt.snp6.process.batch(CEL.list.file = NULL, nthread = 1,
  cluster.type = "PSOCK", ...)
}
\arguments{
  \item{CEL.list.file}{A tab-separated file giving the names of the CEL files and corresponding sample name (see Details section).}
  \item{nthread}{Number of threads for multithreading (see Details).}
  \item{cluster.type}{Type of multithreading cluster (see \code{?parallel::makeCluster()}).}
  \item{...}{Any parameter for the function \code{apt.snp6.process()}, except \code{CEL} and \code{samplename}.}
}
\details{
  \code{CEL.list.file} is a tab-separated text file containing 2 columns :
  \itemize{
    \item{cel_files : Name (and path) of the CEL file(s)}
    \item{SampleName : The output sample name}
  }
}
\author{
  Bastien Job
}

