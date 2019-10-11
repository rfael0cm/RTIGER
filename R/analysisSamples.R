#'
#' Library size and sample wise threshold analysis.
#'
#' @param object an RTIGER object 
#' @return A list with the library size of ech sample as total allele count, and a logical vector which specifies wether or not each samples pass the threshold set when the object was created.
#' 
#' @export analysis.samples
#'

analysis.samples = function(object){
  lib.size = sapply(object@FilteredData$as.GR, function(samp){
    lbs = sum(samp$P1.Allele.Count + samp$P2.Allele.Count, na.rm = TRUE)
    return(lbs)
  })

  quantThr = object@FilteringThreshold$quantile
  min.counts = object@FilteringThreshold$min.counts
  passThresh = sapply(object@FilteredData$as.GR, function(samp){
    counts = samp$total
    pass = quantile(counts, quantThr, na.rm = TRUE) >= min.counts
    return(pass)
  })

  return(list(lib.size = lib.size, passThresh = passThresh))

}
