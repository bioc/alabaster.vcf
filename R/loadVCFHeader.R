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
#' @importFrom VariantAnnotation scanVcfHeader
loadVCFHeader <- function(info, project) {
    header.path <- acquireFile(project, info$path)
    scanVcfHeader(header.path)
}