#' Load a VCF header
#'
#' Load the headers of a VCF file into a \linkS4class{VCFHeader} object.
#'
#' @param info Named list of metadata for a VCF file.
#' @param project Any argument accepted by the acquisition functions, see \code{?\link{acquireFile}}.
#' By default, this should be a string containing the path to a staging directory.
#'
#' @return A \linkS4class{VCFHeader} object.
#'
#' @details
#' As the name suggests, this only loads the headers of the VCF file.
#' To load all contents into memory, use \code{\link{scanVcf}} instead.
#'
#' Users can override this function in \code{\link{loadVCF}} by calling \code{.altLoadVCFHeader}.
#' This should be set to a function that accepts the same arguments as \code{loadVCFHeader} and returns a VCFHeader object.
#' 
#' @author Aaron Lun
#'
#' @examples
#' fl <- system.file("extdata", "structural.vcf", package="VariantAnnotation")
#' hdr <- scanVcfHeader(fl)
#'
#' tmp <- tempfile()
#' dir.create(tmp)
#' info <- stageObject(hdr, dir=tmp, path="header")
#' loadVCFHeader(info, tmp)
#' 
#' @export
#' @aliases .altLoadVCFHeader
#' @importFrom VariantAnnotation scanVcfHeader
loadVCFHeader <- function(info, project) {
    header.path <- acquireFile(project, info$path)
    scanVcfHeader(header.path)
}

header.env <- new.env()
header.env$fun <- loadVCFHeader

#' @export
.altLoadVCFHeader <- function(fun) {
    prev <- header.env$fun
    if (missing(fun)) {
        return(prev)
    } else {
        header.env$fun <- fun
        return(invisible(prev))
    }
}
