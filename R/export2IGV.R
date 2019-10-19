#'
#' Exports the count and viterbi code to a folder in order to visualize it in IGV.
#'
#' @param object an RTIGER or RViterbi object with the filtered data
#' @param sample character with the sample name as specified in the experimentalDesign
#' @param dir character with the chromosome to plot. If left null, it will create a folder with the name of the sample.
#' @param ratio boolean variable. Whether the ratio P1/total should be created as well.
#' 
#' @usage export2IGV(object, sample, dir = NULL, ratio = FALSE)
#' 
#' @examples 
#' 
#' data("R.ViterbiExample")
#' info = myDat@info
#' export2IGV(RVit, samp = info$sample_names[1])
#' @export export2IGV
#'
export2IGV = function( object, sample, dir = NULL, ratio = FALSE){
  if(!is.null(dir)){
    dir.create(dir)
    compdir = file.path(dir, sample)
  }
  else{
    compdir = sample
  }
  # if(dir.exists(dir)) stop("A new directory will be created with the results. \n Please change the name of the directory or delete the existing one.")
  dir.create(compdir)

  P1.statesfile = file.path(compdir, paste("P1-state-", sample, ".bed", sep = ""))
  P2.statesfile = file.path(compdir, paste("P2-state-", sample, ".bed", sep = ""))
  Het.statesfile = file.path(compdir, paste("Het-state-", sample, ".bed", sep = ""))

  P1.countsfile = file.path(compdir, paste("P1-countperbin-", sample, ".bw", sep = ""))
  P2.countsfile = file.path(compdir, paste("P2-countperbin-", sample, ".bw", sep = ""))

  Viterbi = object@Viterbi[[sample]]

  P1.states = Viterbi[Viterbi$Viterbi == "pat"][, "Viterbi"]
  P2.states = Viterbi[Viterbi$Viterbi == "mat"][, "Viterbi"]
  Het.states = Viterbi[Viterbi$Viterbi == "het"][, "Viterbi"]

  P1.count = Viterbi[,"P1.Allele.Count"]
  P1.count$score = P1.count$P1.Allele.Count
  P1.count$score[is.na(P1.count$score)] = 0
  P2.count = Viterbi[,"P2.Allele.Count"]
  P2.count$score = P2.count$P2.Allele.Count
  P2.count$score[is.na(P2.count$score)] = 0
  if(ratio){
    ratiofile = file.path(compdir, paste("Count-ratio-", sample, ".bw", sep = ""))
    ratio = Viterbi$P1.Allele.Count/(Viterbi$P1.Allele.Count + Viterbi$P2.Allele.Count)
    ratio = (ratio - .5)/.5
    ratio[is.na(ratio)] = 0
    Viterbi$score = ratio

    export.bw(Viterbi[,"score"], ratiofile)
  }

  export.bed(P1.states, P1.statesfile)
  export.bed(P2.states, P2.statesfile)
  export.bed(Het.states, Het.statesfile)

  export.bw(P1.count, P1.countsfile)
  export.bw(P2.count, P2.countsfile)
}


#'
#' Exports the count at the marker level to a bw file to visualize in IGV
#'
#' @param object an RTIGER object with the filtered data
#' @param sample character with the sample name as specified in the experimentalDesign
#' @param dir character with the chromosome to plot. If left null, it will create a folder with the name of the sample.
#' @usage snpCounts2IGV(object, sample, dir = NULL)
#' 
#' @examples 
#' 
#' data("fittedExample.rda")
#' info = myDat@info
#' snpCounts2IGV(myDat, samp = info$sample_names[1])
#' @export snpCounts2IGV
#'
snpCounts2IGV = function(object, sample, dir = NULL){
  if(!is.null(dir)){
    dir.create(dir)
    compdir = file.path(dir, sample)
  }
  else{
    compdir = sample
  }

  P1.countsfile = file.path(compdir, paste("P1-countperSNP-", sample, ".bw", sep = ""))
  P2.countsfile = file.path(compdir, paste("P2-countperSNP-", sample, ".bw", sep = ""))

  Raw = object@RawData[[sample]]

  P1.count = Raw[,"P1.Allele.Count"]
  P1.count$score = P1.count$P1.Allele.Count
  P1.count$score[is.na(P1.count$score)] = 0
  P2.count = Raw[,"P2.Allele.Count"]
  P2.count$score = P2.count$P2.Allele.Count
  P2.count$score[is.na(P2.count$score)] = 0

  export.bw(P1.count, P1.countsfile)
  export.bw(P2.count, P2.countsfile)

}
