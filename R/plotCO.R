#'
#' Obtain number of Cross-Over events per sample and chromosome.
#'
#' @param object a RViterbi object.
#' @return Matrix m x n. M number of samples and N chromosomes.
#'
#' #' @return a matrix with n chromosomes and m samples (n x m) and the number of CO events.
#' @usage calcCOnumber(object)
#'
#' @examples
#'
#' data("fittedExample")
#' co.num = calcCOnumber(myDat)
#'
#' @export calcCOnumber
#'


calcCOnumber = function(object){
  numCO = sapply(object@Viterbi, function(samp){
    sapply(seqlevels(samp), function(chr){
      newsamp = samp[seqnames(samp) == chr]
      sampGR = try(Vit2GrangesGen(newsamp, "Viterbi"), silent = TRUE)
      COs = length(sampGR)-1
      return(COs)
    })
  })
  return(numCO)

}

plotCOs = function(object, file){
  Cos = calcCOnumber(object = object)
  Cos = melt(Cos)
  colnames(Cos) = c("Chr", "Sample", "value")
  p <- ggplot(Cos, aes( x = factor(Chr), y = value)) +
    geom_boxplot() +
    # theme(legend.position = "none")+
    xlab("chromosome") +
    ylab("Number of double CO")
  pdf(file)
  print(p)
  dev.off()
}
