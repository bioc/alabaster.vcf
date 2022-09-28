#' Stage a VCF object
#'
#' Save the contents of a \linkS4class{VCF} object to file. 
#'
#' @param x Any instance of a \linkS4class{VCF} class or one of its subclasses.
#' @inheritParams alabaster.base::stageObject
#' @param ... Further arguments to pass to \code{\link{stageObject,RangedSummarizedExperiment-method}}.
#'  
#' @details
#' Note that we do \emph{not} save the contents of \code{x} in VCF format.
#' Rather, we re-use the existing machinery for staging SummarizedExperiments from the \pkg{alabaster.se}.
#' This is more amenable for random access by feature/sample and ensures that we are consistent with the expectations of the parent class.
#' Applications requiring actual VCF files can instead use \code{\link{writeVcf}} to generate them from \code{x}.
#'
#' @author Aaron Lun
#'
#' @return The contents of \code{x} are saved to file inside \code{path}.
#' A named list containing metadata is returned. 
#' 
#' @examples
#' fl <- system.file("extdata", "structural.vcf", package="VariantAnnotation")
#' vcf <- readVcf(fl, genome="hg19")
#'
#' tmp <- tempfile()
#' dir.create(tmp)
#' stageObject(vcf, dir=tmp, path="experiment-1")
#'
#' @export
#' @rdname stageVCF
#' @import methods alabaster.base VariantAnnotation
#' @importFrom S4Vectors metadata<- metadata
#' @importMethodsFrom alabaster.string stageObject
#' @importMethodsFrom alabaster.se stageObject
setMethod('stageObject', "VCF", function(x, dir, path, child=FALSE, ...) {
    dir.create(file.path(dir, path), showWarnings=FALSE)

    stopifnot(is(metadata(x)$header, "VCFHeader"))
    vdeets <- tryCatch({
        vmeta <- .stageObject(metadata(x)$header, dir, file.path(path, "header"))
        .writeMetadata(vmeta, dir=dir)
    }, error=function(e) stop("failed to stage 'header(<", class(x)[1], ">)'\n  - ", e$message))
    metadata(x)$header <- NULL # removing it to avoid restaging in SE's metadata.

    fixed.df <- fixed(x)
    fdeets <- tryCatch({
        fstuff <- .stageObject(fixed.df, dir, file.path(path, "fixed"), child=TRUE)
        .writeMetadata(fstuff, dir=dir)
    }, error=function(e) stop("failed to stage 'fixed(<", class(x)[1], ">)'\n  - ", e$message))
    fixed(x) <- fixed(x)[,0] # stop rowRanges() from appending more columns to the mcols of the GRanges.

    info.df <- info(x, row.names=FALSE)
    ideets <- tryCatch({
        istuff <- .stageObject(info.df, dir, file.path(path, "info"), child=TRUE)
        .writeMetadata(istuff, dir=dir)
    }, error=function(e) stop("failed to stage 'info(<", class(x)[1], ">)'\n  - ", e$message))

    info <- callNextMethod()

    info[["$schema"]] <- "vcf_experiment/v1.json"
    info[["vcf_experiment"]] <- list(
        header = list(resource=vdeets),
        fixed = list(resource=fdeets),
        info = list(resource=ideets),
        expanded = is(x, "ExpandedVCF")
    )

    info
})
