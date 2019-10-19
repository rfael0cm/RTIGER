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
                        showGenAxis = FALSE
){
  DataViterbi_GR = object@Viterbi

  FinalRes_chr = DataViterbi_GR[[samp]]
  FinalRes_chr = FinalRes_chr[seqnames(FinalRes_chr) == chr]

  Viterbi = Vit2GrangesGen(FinalRes_chr, value = "Viterbi")


  dat = FinalRes_chr[,c("P1.Allele.Count", "P2.Allele.Count")]
  dat$P2.Allele.Count = - dat$P2.Allele.Count

  lim = max(abs(as.matrix(mcols(dat))), na.rm = TRUE) + 2

  datgviz = DataTrack(dat, type = "h", name = "Data", col = col[1:2], ylim = c(-lim, lim))

  names(col) = c("pat", "mat", "het")

  myclass = ifelse(class(object) == "RTIGER", "Viterbi","R-Viterbi")
  if(myclass == "R-Viterbi") myclass = paste(myclass, object@Rparameter)
  vitGviz = AnnotationTrack(Viterbi, name = myclass,
                            fill = col[Viterbi$Viterbi],
                            col = "transparent", stacking = "dense")

  mySize = size
  mytracks = c(datgviz, vitGviz)
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




Vit2GrangesGen = function (FinalRes_chr, value)
{
  chr = seqnames(FinalRes_chr)[1]
  myRle = rle(mcols(FinalRes_chr)[[value]])
  len = myRle$lengths
  vals = myRle$values
  cum = len[1]
  if(length(len) > 1){
    for(i in 2:length(len)) cum[i] = cum[i-1] + len[i]
  }

  starts = start(FinalRes_chr)
  ends = end(FinalRes_chr)
  good.ends = ends[cum]
  good.start = vector(length = length(len))
  good.start[1] = starts[1]
  if(length(len) > 1) good.start[2:length(good.start)] = starts[cum[-length(cum)] + 1]
  myGR = data.frame(seqnames = rep(chr, length(good.start)),
                    start = good.start,
                    end = good.ends)
  myGR[[value]] = vals
  myGR = GRanges(myGR)
  return(myGR)
}
