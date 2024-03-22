# Save `VCF`s to file

|Environment|Status|
|---|---|
|[BioC-release](https://bioconductor.org/packages/release/bioc/html/alabaster.vcf.html)|[![Release OK](https://bioconductor.org/shields/build/release/bioc/alabaster.vcf.svg)](http://bioconductor.org/checkResults/release/bioc-LATEST/alabaster.vcf/)|
|[BioC-devel](https://bioconductor.org/packages/devel/bioc/html/alabaster.vcf.html)|[![Devel OK](https://bioconductor.org/shields/build/devel/bioc/alabaster.vcf.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/alabaster.vcf/)|

The **alabaster.vcf** package implements methods for saving and loading `VCF` objects (from the [**VariantAnnotation**](https://bioconductor.org/packages/VariantAnnotation) package) under the **alabaster** framework.
It does so by converting them back into VCF files but with additional metadata to store the number of samples and positions.
To get started, install the package and its dependencies from [Bioconductor](https://bioconductor.org/packages/alabaster.vcf):

```r
# install.packages("BiocManager")
BiocManager::install("alabaster.vcf")
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
saveObject(vcf, tmp)

roundtrip <- readObject(tmp)
class(roundtrip)
## [1] "CollapsedVCF"
## attr(,"package")
## [1] "VariantAnnotation"
```
