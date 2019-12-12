#'
#' This funciton creates matrix with the state count ratio 
#'
#' @param rviterbiobj an RVITERBI object.
#' @param samples the samples used to compute the ratio
#' @return returns a nx3 matrix with n being the number of samples
#' @usage state_ratios(RVit)
#' 
#' @examples 
#' 
#' data("R.ViterbiExample")
#' ratios = state_ratios(RVit)
#' 
#' @export state_ratios
#'
#'


state_ratios = function(rviterbiobj, samples = NULL){
  viterbis = rviterbiobj@Viterbi
  if(is.null(samples)) samples = names(viterbis)
  
  states = c("pat", "het", "mat")
  ratios = sapply(states, function(s){
    sapply(samples, function(samp){
      p1 = viterbis[[samp]]$P1.Allele.Count[viterbis[[samp]]$Viterbi == s]
      p2 = viterbis[[samp]]$P2.Allele.Count[viterbis[[samp]]$Viterbi == s]
      return(sum(p1, na.rm = T)/(sum(p2, na.rm = T) + sum(p1, na.rm = T)))
    })#sapply samples
  }) #sapply states
  return(ratios)
}
