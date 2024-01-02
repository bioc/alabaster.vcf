# This checks that we can successfully stage and load CollapsedVCF objects.
# library(testthat); library(alabaster.vcf); source('setup.R'); source("test-CollapsedVCF.R")

# Just re-using the file from there.
library(VariantAnnotation)
fl <- system.file("extdata", "structural.vcf", package="VariantAnnotation")
vcf <- readVcf(fl)

test_that("Staging an CollapsedVCF works as expected", {
    tmp <- tempfile()
    dir.create(tmp)

    info <- stageObject(vcf, dir=tmp, path="experiment-1")
    expect_identical(dirname(info[["$schema"]]), "vcf_experiment")
    expect_false(info$vcf_experiment$expanded)
    expect_error(.writeMetadata(info, dir=tmp), NA)

    roundtrip <- loadObject(info, tmp)
    expect_s4_class(assay(roundtrip), "DelayedMatrix")
    expect_identical2(vcf, roundtrip)
})

test_that("Saving a CollapsedVCF works for all files", {
    exdir <- system.file("extdata", package="VariantAnnotation")
    all.files <- list.files(exdir, pattern=".vcf(.gz)?$")
    all.files <- setdiff(all.files, "hapmap_exome_chr22.vcf.gz") # TODO: fix in VariantAnnotation.

    for (d in file.path(exdir, all.files)) {
        vcf <- readVcf(d)

        tmp <- tempfile()
        saveObject(vcf, tmp)
        roundtrip <- readObject(tmp)

        # TODO: fix in VariantAnnotation.
        n <- names(metadata(vcf)$header@header)
        metadata(roundtrip)$header@header <- metadata(roundtrip)$header@header[n]
        metadata(roundtrip)$header@reference <- metadata(vcf)$header@reference
        metadata(roundtrip)$header@header$fileDate <- metadata(vcf)$header@header$fileDate

        expect_identical(vcf, roundtrip)
    }
})

test_that("Staging an CollapsedVCF works with higher-dimensional arrays", {
    assay(vcf, "WHEE", withDimnames=FALSE) <- array(runif(prod(dim(vcf)) * 2), c(dim(vcf), 2))

    tmp <- tempfile()
    dir.create(tmp)
    info <- stageObject(vcf, dir=tmp, path="experiment-1")

    ass.meta <- info$summarized_experiment$assays 
    keep <- which(vapply(ass.meta, function(x) x$name == "WHEE", TRUE))
    target <- acquireMetadata(tmp, ass.meta[[keep]]$resource$path)
    expect_identical(dirname(target[["$schema"]]), "hdf5_dense_array")

    roundtrip <- loadObject(info, tmp)
    expect_s4_class(assay(roundtrip, "WHEE"), "DelayedArray")
    expect_identical2(vcf, roundtrip)
})
