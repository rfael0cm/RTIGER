# Quick Guide to RTIGER

## Introduction
Accurate identification of meiotic crossing-over sites (COs) is essential for correct genotyping of recombining samples. RTIGER is a method for predicting genome-wide COs using allele-counts at pre-defined SNP marker positions. RTIGER trains a Hidden Markov Model (HMM) where genomic states (homozygous parent_1, homozygous parent_2 or heterozygous) correspond to the hidden state and the allele-counts as the observed variable. COs are identified as transitions in the HMM state.

To account for variation in the coverage of sequencing data, RTIGER uses Viterbi Path Algorithm and the `rigidity` parameter. This parameter defines the minimum number of SNP markers required to support a state-transition. This filters out low-confidence state-transitions, improving COs identification performance.

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
Neverhteless, to get the latest and debugged versions, we recommend you install it from github using the devtools package:
```
install.packages("devtools")
library(devtools)
install_github("rfael0cm/RTIGER")
```
Once the RTIGER package is installed, next is to install the JULIA libraries needed to run the analysis. We have simplified this step by just running the following function:
```
library(RTIGER)
setupJulia()
```
Now, we have all the necessary prerequisists to start using RTIGER and obtain the maximum information from our data sets. RTIGER needs specific data format to work. On the next section we will describe how it must be enoced and saved.

<!-- ################################################################################ -->
### Preparing input data:
RTIGER uses the allele-count information at the SNP marker positions. The SNP markers correspond to differences between the two genotypes (i.e. parent_1 vs parent_2). RTIGER requires as input one allele-count file for each sample. The allele-count file should be in tab-separated value format, where each row corresponds to a SNP marker. The format of the file is described below:

|Column | Field | Type | Description |
|---|---|---|----|
|1|SeqID|String| Chromosme ID|
|2|Pos|init(>=0)| Position of the SNP marker|
|3|RefA|char| Reference allele|
|4|RefC| int(>=0)| Number of reads with reference alele|
|5|AltA|char|Alternate allele|
|6|AltF|int(>=0)|Number of reads with alternate allele|

Ther order of the columns is **EXTREMELY IMPORTANT**. RTIGER ensures that the data type of each column is the correct. But the interpretation of **references allele** and **alternate allele** is completely arbitrary and it is the user who defines them. Moreover, the chromosome and position is crucial to run our algorithm since we group together consecutive SNPs from the same chromosome.

The SNPs can be identified using any generic SNP identification pipeline. For example look this [method](https://www.ebi.ac.uk/sites/ebi.ac.uk/files/content.ebi.ac.uk/materials/2014/140217_AgriOmics/dan_bolser_snp_calling.pdf.)

SNPs in repetitive regions should be filtered out. Further, as crossing-over usually takes place in syntenic regions between the two genome, for best results, only SNPs in syntenic regions should be selected as markers. If whole genome assemblies are present for both genomes, then this can be easily achieved using methods like [SyRI](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1911-0).

**NOTE 1**: RTIGER assumes that all samples have similar sequencing coverage and, hence, similar distribution of the allele-count values. It does not check or normalise for sequencing coverage variation.

**NOTE 2**: Crossing-over resolution depends on sequenced marker density. Low sequencing coverage could result in few informative markers, which in turn could decrease resolution CO prediction.

**NOTE 3**: RTIGER is designed to be robust against individual outliers, however, the user should check for "bad" markers, i.e. marker positions that are prone to mismapping. These markers result in high allele-count at that position.<!--  as well as in a high, artificial agreement in this position across samples. -->

## Using RTIGER 

To run RTIGER, first we need to invoke julia and the package. To this to work, be sure that everythink run smoothly on the isntallation section.
```
library(RTIGER)
sourceJulia()
```
The function `sourceJulia()` will need to be run every single time that we load RTIGER. The two previous lines of code should be run always.

#### Creating input objects
The primary input for RTIGER is a data-frame termed `expDesign`. The first column of `expDesign` should have paths to allele-count files for all samples and the second column should have unique samples IDs. The columnames should be "files" and "name" to specify the file path and sample name respectivley. We recommend to use full path to the files to avoid further problems.
An example of how to create the data-frame using the example data provided by our package:

```
# Get paths to example allele count files originating from a
# cross between Col-0 and Ler accession of the A.thaliana
file_paths = list.files(system.file("extdata",  package = "RTIGER"), full.names = TRUE)

# Get sample names
sampleIDs <- basename(file_paths)

# Create the expDesign object
expDesign = data.frame(files=file_paths, name=sampleIDs)

print(expDesign)
```
RTIGER also requires chromosome lengths for the parent_1. These need to be provided as a named vector where the values are chromosome lengths and the names are chromosome ids. For our example, we have included the length of the five autosomal chroomosomes in Arabidopsis Thaliana and it can be invoqued as:

```
# Get chromosome lengths for the example data included in the package
chr_len <- RTIGERJ::ATseqlengths
names(chr_len) <- c('Chr1' , 'Chr2', 'Chr3', 'Chr4', 'Chr5')
print(chr_len)
```
**Note**: It is important the that the names of the vector match exactly with the chromosome names in your data files. If any of the files uses another coding for the chromosomes, RTIGER will rise and error. Be consistent!






