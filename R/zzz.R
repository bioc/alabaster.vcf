.onLoad <- function(libname, pkgname) {
    registerReadObjectFunction("vcf_experiment", readVCF)
}

.onUnload <- function(libname, pkgname) {
    registerReadObjectFunction("vcf_experiment", NULL)
}
