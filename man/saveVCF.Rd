% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/saveVCF.R
\name{saveObject,VCF-method}
\alias{saveObject,VCF-method}
\alias{stageObject,VCF-method}
\alias{stageObject,VCFHeader-method}
\title{Save a VCF object to disk}
\usage{
\S4method{saveObject}{VCF}(x, path, ...)
}
\arguments{
\item{x}{Any instance of a \linkS4class{VCF} class or one of its subclasses.}

\item{path}{String containing the path to a directory in which to save \code{x}.}

\item{...}{Further arguments to pass to specific methods.}
}
\value{
\code{x} is saved to file inside \code{path}, and \code{NULL} is returned.
}
\description{
Save a \linkS4class{VCF} object to its on-disk representation, namely a VCF file with the same contents.
}
\examples{
fl <- system.file("extdata", "structural.vcf", package="VariantAnnotation")
vcf <- readVcf(fl)

tmp <- tempfile()
saveObject(vcf, tmp)

}
\seealso{
\code{\link{readVCF}}, to read a VCF object back to the R session.
}
\author{
Aaron Lun
}
