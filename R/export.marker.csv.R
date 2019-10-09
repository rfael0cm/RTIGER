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
