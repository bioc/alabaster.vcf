#' @export
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
