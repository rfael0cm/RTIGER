#'
#' Exports the count and viterbi code to a folder in order to visualize it in IGV.
#'
#' @param object an RTIGER or RViterbi object with the filtered data
#' @param sample character with the sample name as specified in the experimentalDesign
#' @param dir character with the chromosome to plot. If left null, it will create a folder with the name of the sample.
#' @param ratio boolean variable. Whether the ratio P1/total should be created as well.
#' @param newn named vector with the new and shorter names for the samples.
#'
#' @return a folder with files to visualize results in igv.
#' @usage export2IGV(object, sample, dir = NULL, ratio = FALSE, newn = NULL)
#'
#' @examples
#'
#' data("fittedExample")
#' info = myDat@info
#' export2IGV(myDat, samp = info$sample_names[1])
#' @export export2IGV
#'
export2IGV = function( object, sample, dir = NULL, ratio = FALSE, newn = NULL){
  if(!is.null(dir)){
    if(!dir.exists(dir)) dir.create(dir)

    compdir = file.path(dir, sample)
  }
  else{
    compdir = sample
  }

  # if(dir.exists(dir)) stop("A new directory will be created with the results. \n Please change the name of the directory or delete the existing one.")
  if(!dir.exists(compdir)) dir.create(compdir)

  if(!is.null(newn)) sample = newn[sample]

  P1.statesfile = file.path(compdir, paste("P1-state-", sample, ".bed", sep = ""))
  P2.statesfile = file.path(compdir, paste("P2-state-", sample, ".bed", sep = ""))
  Het.statesfile = file.path(compdir, paste("Het-state-", sample, ".bed", sep = ""))
  Comp.statesfile = file.path(compdir, paste("CompleteBlock-state-", sample, ".bed", sep = ""))

  P1.countsfile = file.path(compdir, paste("P1-countperbin-", sample, ".bw", sep = ""))
  P2.countsfile = file.path(compdir, paste("P2-countperbin-", sample, ".bw", sep = ""))

  Viterbi = object@Viterbi[[sample]]
  mycomp = sapply(object@info$part_names, function(x){
    myx = Viterbi[seqnames(Viterbi) == x]
    myx = Vit2GrangesGen(myx, "Viterbi")
    return(myx)
  })
  gr = mycomp[[1]]
  for(i in 2:length(mycomp)){
    gr = c(gr,mycomp[[i]])
  }
  goodanno = c("AA", "AB", "BB")
  names(goodanno) = c("pat", "het", "mat")
  df <- data.frame(seqnames=seqnames(gr),
                   starts=start(gr)-1,
                   ends=end(gr),
                   names=c(rep(".", length(gr))),
                   scores=goodanno[elementMetadata(gr)$Viterbi],
                   strands=strand(gr))
  df = df[,-c(4,6)]


  P1.states = Viterbi[Viterbi$Viterbi == "pat"][, "Viterbi"]
  P2.states = Viterbi[Viterbi$Viterbi == "mat"][, "Viterbi"]
  Het.states = Viterbi[Viterbi$Viterbi == "het"][, "Viterbi"]

  P1.count = Viterbi[,"P1.Allele.Count"]
  P1.count$score = P1.count$P1.Allele.Count
  P1.count$score[is.na(P1.count$score)] = 0
  P2.count = Viterbi[,"total"]
  P2.count$score = Viterbi$total - Viterbi$P1.Allele.Count
  # P2.count = Viterbi[,"P2.Allele.Count"]
  # P2.count$score = P2.count$P2.Allele.Count
  P2.count$score[is.na(P2.count$score)] = 0
  if(ratio){
    ratiofile = file.path(compdir, paste("Count-ratio-", sample, ".bw", sep = ""))
    ratio = Viterbi$P1.Allele.Count/Viterbi$total
    ratio = (ratio - .5)/.5
    ratio[is.na(ratio)] = 0
    Viterbi$score = ratio

    export.bw(Viterbi[,"score"], ratiofile)
  }
  write.table(df, file= Comp.statesfile, quote=F, sep="\t", row.names=F, col.names=F)
  export.bed(P1.states, P1.statesfile)
  export.bed(P2.states, P2.statesfile)
  export.bed(Het.states, Het.statesfile)

  export.bw(P1.count, P1.countsfile)
  export.bw(P2.count, P2.countsfile)
}
