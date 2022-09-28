#' Stage the VCF headers
#'
#' Save the contents of a \linkS4class{VCFHeader} object to file.
#' This is formatted as a valid VCF file that lacks any entries.
#'
#' @param x A \linkS4class{VCFHeader} object.
#' @inheritParams alabaster.base::stageObject
#'
#' @author Aaron Lun
#'
#' @return The contents of \code{x} are saved to file inside \code{path}.
#' A named list containing metadata is returned. 
#' 
#' @examples
#' fl <- system.file("extdata", "structural.vcf", package="VariantAnnotation")
#' hdr <- scanVcfHeader(fl)
#'
#' tmp <- tempfile()
#' dir.create(tmp)
#' stageObject(hdr, dir=tmp, path="headers")
#'
#' @rdname stageVCFHeader
#' @importFrom Rsamtools bgzip
setMethod("stageObject", "VCFHeader", function(x, dir, path, child=FALSE) { 
    dir.create(file.path(dir, path), showWarnings=FALSE)

    header.lines <- .create_vcf_header(x)
    tmp.file <- tempfile()
    write(file=tmp.file, unlist(header.lines))
    header.path <- file.path(path, "header.bcf")
    bgzip(tmp.file, dest=file.path(dir, header.path))

    list(
        `$schema` = "vcf_file/v1.json",
        path = header.path,
        is_child = child,
        vcf_file = list(
            compression="bgzip",
            header_only=TRUE
        )
    )
})

#' @importFrom VariantAnnotation samples geno
.create_vcf_header <- function(hdr) {
    fileformat <- "fileformat" %in% rownames(meta(hdr)$META)
    if (!fileformat) {
        fileformat <- "fileformat" %in% names(meta(hdr))
    }
    if (fileformat && grepl(fileformat, "v4.2", fixed = TRUE) || !fileformat) {
        if (any(idx <- rownames(geno(hdr)) == "AD")) {
            geno(hdr)[idx, ]$Number <- "G"
        }
    }

    dflist <- header(hdr)
    header.lines <- Map(.format_header, as.list(dflist), as.list(names(dflist)))
    idx <- which(names(header.lines) == "fileformat")
    if (length(idx) && idx != 1) {
        header.lines <- c(header.lines[idx], header.lines[-idx])
    }

    colnms <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
    if (nrow(geno(hdr)) > 0) {
        colnms <- c(colnms, "FORMAT", samples(hdr))
    }

    c(header.lines, paste(colnms, collapse = "\t"))
}

#' @importFrom S4Vectors DataFrame
.format_header <- function (df, nms) {
    if (nms == "META" && ncol(df) == 1L) {
        paste("##", rownames(df), "=", df[, 1], sep = "")
    } else if (nms == "PEDIGREE" || nms == "ALT") {
        if (!is.null(rownames(df))) {
            df <- DataFrame(ID = rownames(df), df)
        }
        if ("Description" %in% colnames(df)) { # Missing from VariantAnnotation's function, causes loss of ALT information.
            if (nrow(df) != 0L) {
                df$Description <- ifelse(is.na(df$Description), "\".\"", paste("\"", df$Description, "\"", sep = ""))
            }
        }
        .paste_multi_field_df(df, nms)
    } else if (ncol(df) == 1L && nrow(df) == 1L) {
        paste("##", nms, "=", df[, 1], sep = "")
    } else {
        if ("Description" %in% colnames(df)) {
            if (nrow(df) == 0L) {
                return(character())
            }
            df$Description <- ifelse(is.na(df$Description), "\".\"", paste("\"", df$Description, "\"", sep = ""))
        }
        df <- DataFrame(ID = rownames(df), df)
        .paste_multi_field_df(df, nms)
    }
}

#' @importFrom S4Vectors unstrsplit
.paste_multi_field_df <- function (df, nms) {
    if (nrow(df) == 0L) {
        character(0L)
    } else {
        prs <- paste(rep(colnames(df), each = nrow(df)), "=", unlist(lapply(df, as.character), use.names = FALSE), sep = "")
        lst <- split(prs, row(df))
        lns <- unstrsplit(lst, ",")
        paste("##", nms, "=<", lns, ">", sep = "")
    }
}
