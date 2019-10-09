#### Plot functions for each chromosome and sample
#### 27/04/2019
#### by Rafael Campos Martin


## TODO: When there are no CO plotting breaks

plotGenotype = function(object,
                        samp,
                        chr,
                        col = c("red", "blue", "mediumorchid4"),
                        size = c(1,1),
                        cex.feature = .7,
                        from = NULL,
                        to = NULL,
                        main = NULL,
                        showId = TRUE,
                        showGenAxis = FALSE
){
  DataViterbi_GR = object@Viterbi

  FinalRes_chr = DataViterbi_GR[[samp]]
  FinalRes_chr = FinalRes_chr[seqnames(FinalRes_chr) == chr]

  Viterbi = Vit2GrangesGen(FinalRes_chr, value = "Viterbi")


  dat = FinalRes_chr[,c("P1.Allele.Count", "P2.Allele.Count")]
  dat$P2.Allele.Count = - dat$P2.Allele.Count

  lim = max(abs(as.matrix(mcols(dat))), na.rm = T) + 2

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


  plotTracks(mytracks, groups = colnames(mcols(dat)),
             cex.feature= cex.feature, background.title="darkgrey", lwd=2,
             # from=1,
             to= seqlengths(FinalRes_chr)[chr] ,
             # to = 1733000,
             main = paste(samp, "\n", chr),
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
