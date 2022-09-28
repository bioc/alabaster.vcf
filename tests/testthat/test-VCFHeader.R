# Staging and loading the VCFHeader works as expected.
# library(testthat); library(alabaster.vcf); source('setup.R'); source("test-VCFHeader.R")

fl <- system.file("extdata", "structural.vcf", package="VariantAnnotation")
hdr <- VariantAnnotation::scanVcfHeader(fl)

test_that("VCFHeader stage/load works as expected", {
    tmp <- tempfile()
    dir.create(tmp)

    info <- stageObject(hdr, dir=tmp, path="header")
    expect_error(.writeMetadata(info, dir=tmp), NA)

    roundtrip <- loadVCFHeader(info, tmp)
    expect_identical(hdr, roundtrip)
})
