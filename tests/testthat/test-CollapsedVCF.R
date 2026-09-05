# library(testthat); library(alabaster.vcf); source('setup.R'); source("test-CollapsedVCF.R")

# Just re-using the file from there.
library(VariantAnnotation)
fl <- system.file("extdata", "structural.vcf", package="VariantAnnotation")
vcf <- readVcf(fl)

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
