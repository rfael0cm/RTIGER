#### Plot functions for each chromosome and sample
#### 27/04/2019
#### by Rafael Campos Martin


#'
#' This funciton creates an RViterbi object. 
#'
#' @param object an RTIGER or RViterbi object with the filtered data
#' @param samp character with the sample name as specified in the experimentalDesign
#' @param chr character with the chromosome to plot
#' @param col color to define segments in P1, P2 and heterozygous regions.
#' @param size size of each  panel (see Gviz documentation)
#' @param main an overall title for the plot. If null, the sample and chromosome are the default.
#' @param showGenAxis boolean parameter whether or not show the genome axis.
#' @param ratio boolean parameter wether or not to plot the ratio counts
#' @param window if ratio TRUE, what is the window for the sliding window to be computed of the ratios.
#' 
#' @return plot with the genomic segmentation for chromosome chr in sample samp.
#' @usage plotGenotype(object, samp, chr, col = c("red", "blue", "mediumorchid4"), size = c(1,1), main = NULL, showGenAxis = FALSE)
#' 
#' @examples 
#' 
#' data("R.ViterbiExample")
#' info = myDat@info
#' plotGenotype(RVit, samp = info$sample_names[1], info$part_names[1])
#' @export plotGenotype
#'
#'

plotGenotype = function(object,
                        samp,
                        chr,
                        col = c("red", "blue", "mediumorchid4"),
                        size = c(1,1),
                        main = NULL,
                        showGenAxis = FALSE,
                        ratio = FALSE,
                        window = 10
){
  DataViterbi_GR = object@Viterbi
  
  FinalRes_chr = DataViterbi_GR[[samp]]
  FinalRes_chr = FinalRes_chr[seqnames(FinalRes_chr) == chr]
  
  if(ratio){
    ratio_V = FinalRes_chr$P1.Allele.Count/(FinalRes_chr$P1.Allele.Count + FinalRes_chr$P2.Allele.Count)
    ratio_V = ((ratio_V - .5)/.5)*100
    ratio_V = myrunningMean(ratio_V, window)
    posRatio = ratio_V
    posRatio[posRatio < 0] = NA
    negRatio = ratio_V
    negRatio[negRatio >= 0 ] = NA
    "P1.Allele.Count"
    # rat io[is.na(ratio)] = 0
    # FinalRes_chr$ratio = ratio_V
    
    myratio = FinalRes_chr
    myratio$P1.Allele.Count = posRatio
    myratio$P2.Allele.Count = negRatio
    myratio = myratio[,c("P1.Allele.Count", "P2.Allele.Count")]
    ratlim = c(-100,100)
    ratgviz = DataTrack(myratio, type = "l", col = col[1:2], name = "Ratio", ylim = ratlim)
    size = c(1,1,1)
    
  }
  
  Viterbi = RTIGER:::Vit2GrangesGen(FinalRes_chr, value = "Viterbi")
  
  
  dat = FinalRes_chr[,c("P1.Allele.Count", "P2.Allele.Count")]
  dat$P2.Allele.Count = - dat$P2.Allele.Count
  
  lim = max(abs(as.matrix(mcols(dat))), na.rm = TRUE) + 2
  
  datgviz = DataTrack(dat, type = "h", name = "Data", col = col[1:2], ylim = c(-lim, lim))
  
  names(col) = c("pat", "mat", "het")
  
  myclass = ifelse(is(object,"RTIGER"), "Viterbi","R-Viterbi")
  if(myclass == "R-Viterbi") myclass = paste(myclass, object@Rparameter)
  vitGviz = AnnotationTrack(Viterbi, name = myclass,
                            fill = col[Viterbi$Viterbi],
                            col = "transparent", stacking = "dense")
  
  mySize = size
  mytracks =  c(datgviz, vitGviz)
  if(ratio){
    mytracks = c(datgviz, ratgviz, vitGviz)
  }
  
  if(showGenAxis){
    gtrack <- GenomeAxisTrack()
    mytracks = c(gtrack, mytracks)
    mySize = c(.5,size)
    
  }
  
  if(is.null(main)) main = paste(samp, "\n", chr)
  plotTracks(mytracks, groups = colnames(mcols(dat)),
             background.title="darkgrey", lwd=2,
             to= seqlengths(FinalRes_chr)[chr] ,
             main = main,
             sizes=mySize,
             showFeatureId=FALSE,
             fontcolor.feature="black",  background.title="darkgrey",
             showId=TRUE)
}

myrunningMean = function(x, winHalfSize = 2) 
{
  out = x
  for (i in 1:length(x)) {
    start = (max(c(1, i - winHalfSize)))
    end = (min(c(length(x), i + winHalfSize)))
    out[i] = mean(x[start:end], na.rm = TRUE)
  }
  out
}
