# Quick Guide to RTIGER

## Introduction
Accurate identification of meiotic crossing-over sites (COs) is essential for correct genotyping of recombining samples. RTIGER is a method for predicting genome-wide COs using allele-counts at pre-defined SNP marker positions. RTIGER trains a Hidden Markov Model (HMM) where genomic states (homozygous parent_1, homozygous parent_2 or heterozygous) correspond to the hidden state and the allele-counts as the observed variable. COs are identified as transitions in the HMM state.

To account for variation in the coverage of sequencing data, RTIGER uses Viterbi Path Algorithm and the `rigidity` parameter. This parameter defines the minimum number of SNP markers required to support a state-transition. This filters out low-confidence state-transitions, improving COs identification performance.
\

<!-- ################################################################################ -->
## Installation
#### Pre-Requisites:
* R: Version > X_X_X
* Julia-1.0.5 (Which versions of Julia can be supported?): Julia needs to be installed and available in the environment^[https://www.geeksforgeeks.org/how-to-setup-julia-path-to-environment-variable/?ref=lbp]
* REQUIRED R LIBRARIES WOULD BE INSTALLED AUTOMATICALLY... RIGHT?

\
