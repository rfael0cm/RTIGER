#### 05.06.2019
#### csv output as Laia likes it

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
