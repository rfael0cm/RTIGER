# Quick Guide to RTIGER

## Introduction
Accurate identification of meiotic crossing-over sites (COs) is essential for correct genotyping of recombining samples. RTIGER is a method for predicting genome-wide COs using allele-counts at pre-defined SNP marker positions. RTIGER trains a Hidden Markov Model (HMM) where genomic states (homozygous parent_1, homozygous parent_2 or heterozygous) correspond to the hidden state and the allele-counts as the observed variable. COs are identified as transitions in the HMM state.

To account for variation in the coverage of sequencing data, RTIGER uses Viterbi Path Algorithm and the `rigidity` parameter. This parameter defines the minimum number of SNP markers required to support a state-transition. This filters out low-confidence state-transitions, improving COs identification performance.
\

<!-- ################################################################################ -->
## Installation
#### Pre-Requisites:
* R: Version > 3.6
* **RECOMENDED**: Julia-1.0.5 (Which versions of Julia can be supported?): Julia needs to be installed and available in the environment [Link to Julia](https://www.geeksforgeeks.org/how-to-setup-julia-path-to-environment-variable/?ref=lbp)
  - Julia 1.4.1 Works as well but slightly slower.
  - REQUIRED R LIBRARIES for Julia will be installed from the R terminal just after installation of the package.

RTIGER is an R package meant to be used for Biological data. For that reason, we rely in some Bioconductor libraries that needs to be installed beofrehand to ensure a smooth installation of the package from CRAN. To install Bioconductor and the packages needed, run on the R terminal:
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.14")

BiocManager::install(c("GenomicRanges", "GenomeInfoDb", "TailRank", "IRanges", "Gviz"))

```
Now we are ready to install RTIGER. As mentioned before, RTIGER is an CRAN package and can be installed with the simple command line:
```
install.packages("RTIGER")
```

<!-- ################################################################################ -->
### Preparing input data:
RTIGER uses the allele-count information at the SNP marker positions. The SNP markers correspond to differences between the two genotypes (i.e. parent_1 vs parent_2). RTIGER requires as input one allele-count file for each sample. The allele-count file should be in tab-separated value format, where each row corresponds to a SNP marker. The format of the file is described below:



