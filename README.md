# Save `VCF`s to file

The **alabaster.vcf** package implements methods for saving and loading `VCF` objects (from the [**VariantAnnotation**](https://bioconductor.org/packages/VariantAnnotation) package) under the **alabaster** framework.
Note that it doesn't save the VCF files themselves, but rather, the `SummarizedExperiment`-based representation of those files in the R session.
To get started, install the package and its dependencies from GitHub:

```r
devtools::install_github("ArtifactDB/alabaster.schemas")
devtools::install_github("ArtifactDB/alabaster.base")
devtools::install_github("ArtifactDB/alabaster.ranges")
devtools::install_github("ArtifactDB/alabaster.matrix")
devtools::install_github("ArtifactDB/alabaster.se")
devtools::install_github("ArtifactDB/alabaster.string")
devtools::install_github("ArtifactDB/alabaster.vcf")
```

In the example below, we save a `VCF` object to file:

```r
library(VariantAnnotation)
fl <- system.file("extdata", "structural.vcf", package="VariantAnnotation")
vcf <- readVcf(fl, genome="hg19")
vcf
## class: CollapsedVCF
## dim: 7 1
## rowRanges(vcf):
##   GRanges with 5 metadata columns: paramRangeID, REF, ALT, QUAL, FILTER
## info(vcf):
##   DataFrame with 10 columns: BKPTID, CIEND, CIPOS, END, HOMLEN, HOMSEQ, IMPR...
## info(header(vcf)):
##              Number Type    Description
##    BKPTID    .      String  ID of the assembled alternate allele in the asse...
##    CIEND     2      Integer Confidence interval around END for imprecise var...
##    CIPOS     2      Integer Confidence interval around POS for imprecise var...
##    END       1      Integer End position of the variant described in this re...
##    HOMLEN    .      Integer Length of base pair identical micro-homology at ...
##    HOMSEQ    .      String  Sequence of base pair identical micro-homology a...
##    IMPRECISE 0      Flag    Imprecise structural variation
##    MEINFO    4      String  Mobile element info of the form NAME,START,END,P...
##    SVLEN     .      Integer Difference in length between REF and ALT alleles
##    SVTYPE    1      String  Type of structural variant
## geno(vcf):
##   List of length 4: GT, GQ, CN, CNQ
## geno(header(vcf)):
##        Number Type    Description
##    GT  1      String  Genotype
##    GQ  1      Float   Genotype quality
##    CN  1      Integer Copy number genotype for imprecise events
##    CNQ 1      Float   Copy number genotype quality for imprecise events

library(alabaster.vcf)
tmp <- tempfile()
dir.create(tmp)
meta <- stageObject(vcf, tmp, "vcf")
meta[["$schema"]]
## [1] "vcf_experiment/v1.json"

roundtrip <- loadObject(meta, tmp)
class(roundtrip)
## [1] "CollapsedVCF"
## attr(,"package")
## [1] "VariantAnnotation"
```
