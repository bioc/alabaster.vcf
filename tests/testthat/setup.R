expect_identical2 <- function(ref, obs) {
    for (x in assayNames(obs)) {
        assay(obs, x, withDimnames=FALSE) <- as.array(assay(obs, x, withDimnames=FALSE))
    }
    fixed(ref)$REF <- DNAStringSet(as.character(fixed(ref)$REF))
    expect_identical(ref, obs)
}
