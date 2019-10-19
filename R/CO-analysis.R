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
#' data("R.ViterbiExample")
#' co.num = calcCOnumber(RVit)
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


#'
#' Haplotypes genomic widths.
#'
#' @param object a RViterbi object.
#' 
#' @export haplotypes.width
#'

haplotypes.width = function(object){
  widthCO = lapply(object@Viterbi, function(samp){
    sapply(seqlevels(samp), function(chr){
      newsamp = samp[seqnames(samp) == chr]
      sampGR = try(Vit2GrangesGen(newsamp, "Viterbi"), silent = TRUE)
      if(class(sampGR) == "try-error") COs = 0
      else{
        COs = width(sampGR)
      }
      # COs = ifelse(class(sampGR) == "try-error", 0,  width(sampGR))
      return(COs)
    })
  })
  return(widthCO)

}

#'
#' Compute number of double CO per sample and chromosome.
#'
#' @param object a RViterbi object.
#' @param samples character vector with the name of the samples.
#' 
#' @return a matrix with n chromosomes and m samples (n x m) and the number of double CO events.
#' @usage calcDoubleCOnumber(object, samples = NULL)
#' 
#' @examples 
#' 
#' data("R.ViterbiExample")
#' dco.num = calcDoubleCOnumber(RVit)
#'
#' @export calcDoubleCOnumber
#
calcDoubleCOnumber = function(object, samples = NULL){
  if(is.null(samples)) samples = names(object@Viterbi)

  numdco = sapply(samples, function(samp){
    samp2 = samp
    samp = object@Viterbi[[samp]]

    dcosamp = sapply(seqlevels(samp), function(chr){
      newsamp = samp[seqnames(samp) == chr]
      vit = try(Vit2GrangesGen(newsamp, "Viterbi"), silent = TRUE)

      if(class(vit) == "try-error") hetpos = 0
      else{
        hetpos = which(vit$Viterbi == "het")
      } # end else

      if( all((length(hetpos) == 1) & (hetpos == 0 | hetpos == 1))) return(dco = 0 )
      else{
        if(hetpos[1] == 1) hetpos = hetpos[-1]
        if(hetpos[length(hetpos)] == length(vit)) hetpos = hetpos[-length(hetpos)]

        dco = unlist(sapply(hetpos, function(i){
          if(vit$Viterbi[i-1] == vit$Viterbi[i + 1]) return(i)
        }) # end sapply dco
        )# end unlist
      }# end else
      return(length(dco))

    }) # dcosamp sapply
    return(dcosamp)
  }) #numdco sapply

  return(numdco)

}


#'
#' genomic widths of the double COs.
#'
#' @param object a RViterbi object.
#' @param samples character vector with the name of the samples.
#' 
#' @export width.DoubleCO
#


width.DoubleCO = function(object, samples = NULL){
  if(is.null(samples)) samples = names(object@Viterbi)
  
  numdco = sapply(samples, function(samp){
    samp2 = samp
    samp = object@Viterbi[[samp]]
    
    dcosamp = sapply(seqlevels(samp), function(chr){
      newsamp = samp[seqnames(samp) == chr]
      vit = try(Vit2GrangesGen(newsamp, "Viterbi"), silent = TRUE)
      
      if(class(vit) == "try-error") hetpos = 0
      else{
        hetpos = which(vit$Viterbi == "het")
      } # end else
      
      if( all((length(hetpos) == 1) & (hetpos == 0 | hetpos == 1))) return(dco = 0 )
      else{
        if(hetpos[1] == 1) hetpos = hetpos[-1]
        if(hetpos[length(hetpos)] == length(vit)) hetpos = hetpos[-length(hetpos)]
        
        dco = unlist(sapply(hetpos, function(i){
          if(vit$Viterbi[i-1] == vit$Viterbi[i + 1]) return(i)
        }) # end sapply dco
        )# end unlist
      }# end else
      res = width(vit[dco])
      return(res)
      
    }) # dcosamp sapply
    return(dcosamp)
  }) #numdco sapply
  numdco = unlist(numdco)
  return(numdco)
  
}

