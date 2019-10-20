#' Checks if the txt files from input have the correct format
#'
#' @keywords internal
#' @noRd
#' 
checkfileColumns = function(f){
  if(!is(f$V2, "integer")){
    stop(cat(paste("Second column in file ", samp, " is not an integer.\nGenomic positions must be integers!" )))
  }

  if(is(f$V4, "integer") | is(f$V6, "integer")){
    stop(cate(paste("Either one or both of the allele counts in file ", samp, " are not integers.\nAllele counts must be integers!")))
  } # If values are integers

  f <- f[rowSums(f[,c(4,6)]) != 0, ]
  f <- f[!duplicated(f$V2),]

  return(f)
}

#' Checks the criterion in datasetimportFromtxt
#'
#' @keywords internal
#' @noRd
#' 

criterion = function(x, quant, min.counts){
  crit = quantile(x[x > 0], quant, na.rm = TRUE) < min.counts
  return(crit)
}

#' Bins the genome and sums the alleles 
#'
#' @keywords internal
#' @noRd
#' 


binningFun = function(myG, bin.length){
  bins <- tileGenome(seqinfo(myG), tilewidth= bin.length,
                     cut.last.tile.in.chrom=TRUE)

  hits = findOverlaps(bins, myG)

  agg = aggregate(myG, hits, P1.Allele.Count = sum(P1.Allele.Count, na.rm = TRUE),
                  total = sum(P1.Allele.Count + P2.Allele.Count, na.rm = TRUE))

  bins$P1.Allele.Count[countQueryHits(hits) > 0L] = agg$P1.Allele.Count

  bins$total[countQueryHits(hits) > 0L] = agg$total

  # bins = bins[!is.na(bins$total)]

  return(bins)
}

