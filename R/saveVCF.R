#' Save a VCF object to disk
#'
#' Save a \linkS4class{VCF} object to its on-disk representation, namely a VCF file with the same contents. 
#'
#' @param x Any instance of a \linkS4class{VCF} class or one of its subclasses.
#' @inheritParams alabaster.base::saveObject
#'  
#' @author Aaron Lun
#'
#' @return \code{x} is saved to file inside \code{path}, and \code{NULL} is returned.
#'
#' @seealso
#' \code{\link{readVCF}}, to read a VCF object back to the R session.
#' 
#' @examples
#' fl <- system.file("extdata", "structural.vcf", package="VariantAnnotation")
#' vcf <- readVcf(fl)
#'
#' tmp <- tempfile()
#' saveObject(vcf, tmp)
#'
#' @export
#' @rdname saveVCF 
#' @aliases 
#' stageObject,VCF-method
#' stageObject,VCFHeader-method
#' @import methods alabaster.base VariantAnnotation
setMethod("saveObject", "VCF", function(x, path, ...) {
    dir.create(path, showWarnings=FALSE)

    # Saving as bgzip for now, as VariantAnnotation::readVcf refuses to work
    # with regularly Gzipped files, see Bioconductor/VariantAnnotation#32.
    tmp <- tempfile()
    out <- writeVcf(x, tmp, index=TRUE)
    file.rename(out, file.path(path, "file.vcf.gz"))

    saveObjectFile(path, "vcf_experiment", list(vcf_experiment=list(version="1.0", dimensions=dim(x), expanded=is(x, "ExpandedVCF"))))
    invisible(NULL)
})

########################
#### OLD STUFF HERE ####
########################

#' @export
#' @import alabaster.string
#' @importFrom S4Vectors metadata metadata<-
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
