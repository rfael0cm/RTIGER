#### 05.06.2019
#### csv output as Laia likes it


#'
#' write.geno.csv prints its required Viterbi path from an object (after converting it to a data frame if it is not one nor a matrix) to a file or connection.
#'
#' @param object A Rtiger or RViterbi object
#' @param file ither a character string naming a file or a connection open for writing. "" indicates output to the console.
#' @param samples a set of samples. If null, all samples will be used
#' @param chr The Chromosome to use.
#' @return List of genomic Ranges object
#' 
#' 
#' @usage write.geno.csv(object, file, samples = NULL, chr = NULL)
#' 
#' @examples 
#' 
#' data("R.ViterbiExample")
#' write.geno.csv(RVit, file = "example.geno.csv", chr = 1)
#' 
#' @export write.geno.csv
#'

write.geno.csv = function(object, file, samples = NULL, chr = NULL){
  if(!is.rtiger(object)) stop("The object is not an R-TIGER or RViterbi object")
  viterbis = object@Viterbi
  if(is.null(samples)) samples = names(viterbis)
  viterbis = viterbis[samples]
  if(is.null(chr)){
    res = sapply(viterbis, function(x) x$Viterbi)
    rn = as.data.frame(granges(viterbis[[1]]))[,1:3]
    # rn[1,2]
    rn$ranges = paste(rn$start, rn$end, sep = ":")
    rn = rn[, c(1,4)]
    res = t(cbind(rn, res))
  }
  else{
    res = sapply(viterbis, function(x) x$Viterbi[as.vector(seqnames(x)) == chr])
    rn = as.data.frame(granges(viterbis[[1]])[as.vector(seqnames(viterbis[[1]])) == chr])[,1:3]
    # rn[1,2]
    rn$ranges = paste(rn$start, rn$end, sep = ":")
    rn = rn[, c(1,4)]
    res = t(cbind(rn, res))
  }
  # cat(dim(res))
  write.csv(res, file = file)

}


is.rtiger = function(object){
  res = class(object) == "RTIGER"
  res2 = class(object) == "RViterbi"
  fres = res | res2

  return(fres)
}

#'
#' write.geno.csv prints its required Viterbi path in the marker level from an object (after converting it to a data frame if it is not one nor a matrix) to a file or connection.
#'
#' @param rtigerobj A Rtiger object.
#' @param rviterbiobj A RViterbi object.
#' @param file ither a character string naming a file or a connection open for writing. "" indicates output to the console.
#' @param samples a set of samples. If null, all samples will be used
#' @param chr The Chromosome to use.
#' @return List of genomic Ranges object
#' 
#' 
#' @usage write.marker.geno.csv(rtigerobj, rviterbiobj, file, samples = NULL, chr = NULL)
#' 
#' @examples 
#' 
#' data("fittedExample")
#' data("R.ViterbiExample")
#' write.marker.geno.csv(myDat, RVit, file = "example.geno.csv", chr = 1)
#' 
#' @export write.geno.csv
#'


write.marker.geno.csv = function(rtigerobj, rviterbiobj, file, samples = NULL, chr = NULL){
  if(class(rtigerobj) != "RTIGER") stop("rtiger is not a RTIGER object")
  if(class(rviterbiobj) != "RViterbi") stop("rviterbiobj is not a RViterbi")
  viterbis = rviterbiobj@Viterbi
  mark.gr = rtigerobj@RawData
  
  if(is.null(samples)) samples = names(viterbis)
  viterbis = viterbis[samples]
  mark.gr = mark.gr[samples]
  
  mark.cum = GenomicRanges::granges(mark.gr[[1]])
  for(i in 2:length(mark.gr)){
    mark.cum = c(mark.cum, GenomicRanges::granges(mark.gr[[i]]))
    mark.cum = unique(mark.cum)
    mark.cum = sort(mark.cum)
  }
  
  if(is.null(chr)){
    res = sapply(viterbis, function(x){
      
    } )
    rn = as.data.frame(granges(viterbis[[1]]))[,1:3]
    # rn[1,2]
    rn$ranges = paste(rn$start, rn$end, sep = ":")
    rn = rn[, c(1,4)]
    res = t(cbind(rn, res))
  }
  else{
    mark.cum = mark.cum[as.vector(seqnames(mark.cum)) == chr]
    res = sapply(viterbis, function(x){
      tx = x[as.vector(seqnames(x)) == chr]
      tx = Vit2GrangesGen(tx, "Viterbi")
      mx = findOverlaps(tx, mark.cum)
      mymark = mark.cum
      mymark$Viterbi = NA
      mymark$Viterbi[mx@to] = tx$Viterbi[mx@from]
      return(mymark$Viterbi)
    } )
    rn = as.data.frame(mark.cum)
    rn$ranges = paste(rn$start, rn$end, sep = ":")
    rn = rn[, c("seqnames","ranges")]
    res = cbind(rn, res)
  }
  
  write.csv(res, file = file)
}

