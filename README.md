
# RUVIII-PRPS

<!-- badges: start -->
<!-- badges: end -->

RUVIII-PRPS: Removing unwanted variation from large-scale RNA sequencing data with PRPS
Please cite the appropriate article when you use results from the software in a publication. Such citations are the main means by which the authors receive professional credit for their work.

The RUVIII-PRPS software package itself can be cited as:

Molania R, Foroutan M, Gagnon-Bartsch JA, Gandolfo LC, Jain A, Sinha A, Olshansky G, Dobrovic A, Papenfuss AT, Speed TP. Removing unwanted variation from large-scale RNA sequencing data with PRPS. Nat Biotechnol. 2023 Jan;41(1):82-95. doi: 10.1038/s41587-022-01440-w. Epub 2022 Sep 15. PMID: 36109686; PMCID: PMC9849124.

##  RUVIII-PRPS Installation

After installing the dependent libraries, RUVIII-PRPS can be installed by running the following lines

``` r
library(devtools)
devtools::install_github(
    repo = 'RMolania/RUVIIIPRPS',
    force = TRUE,
    build_vignettes = FALSE)
```

## Using RUVIII-PRPS to remove unwanted variation from large-scale RNA sequencing data

RUVIII-PRPS is a novel strategy using pseudo-replicates of pseudo-samples (PRPS) to normalize RNA-seq data in situations when technical replicate are not available or well-designed. 
We provided a vignette Introduction_to_RUVIII-PRPS.Rmd that explains step by step how to load and normalise datasets and also how to visualise the diagnostic plots before and after normalisation.


<img src="MetaData.pdf">

##  Note
The RUVIIIPRPS is under final preparation, a manuscript and comprehensive vignettes will be available soon.