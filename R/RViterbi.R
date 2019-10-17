
#' This class is a generic container for R-Viterbi. It contains the modified decoding version of the Viterbi algorithm
#'
#' @slot Rparameter The R parameter. It specifies the minimum number of iterations over a state before a transition is permited.
#' @slot LogLike loglikelihood of the observed viterbi path given the observations.
#' @slot Viterbi List of chromosomes with the viterbi path per sample
#' @rdname RViterbi
#' @exportClass RViterbi

.RViterbi = setClass("RViterbi",
                   representation = representation(
                     Rparameter = "numeric",
                     LogLike = "numeric",
                     Viterbi = "list"
                   ))


#'
#' This funciton creates an RViterbi object. 
#'
#' @param object an RTIGER object with the filtered data
#' @param rigidity rigidity parameter. Number of consecutive iterations in state s before using the transition matrix.
#' #' @return returns and RViterbi object with the new viterbi path with rigidity R.
#' @usage R.Viterbi(object, rigidity = NULL)
#' 
#' @examples 
#' 
#' data("fittedExample")
#' RVit = R.Viterbi(myDat, rigidity = 3)
#' 
#' @export R.Viterbi
#'
#'


R.Viterbi = function(object, rigidity = NULL){
  rigid = rigidity

  info = object@info
  sampleN = info$sample_names
  chrN = info$part_names
  chrL = info$part_lengths
  anno = seqlengths(object@FilteredData$as.GR[[1]])

  if (rigid==1) {break("Standard Viterbi is computing after model fitting.")}

  tiles = width(object@FilteredData$as.GR[[1]])[1]
  goodObs = object@FilteredData$as.mat
  psis = object@Model$psi
  transmat = object@Model$params$a
  pivec = object@Model$params$piv

  newvit = lapply(psis, function(chr){
    chrvit = apply(chr, 1, function(samp) rigid_Viterbi(psimat = t(samp[, c("pat", "het", "mat")]), transmat = transmat, pivec = pivec, rigid = rigid))
    # rownames(chrvit) = dimnames(chr)[[2]]
    return(chrvit)
  }) #lapply pisis

  loglike = sum(unlist(lapply(newvit, function(chr) unlist(lapply(chr, function(samp) samp$loglikelihood)))))

  vits = vector("list", length(chrN))
  names(vits) = chrN
  for(chr in chrN){
    vits[[chr]] = matrix(NA, nrow = dim(psis[[chr]])[2],
                         ncol = length(sampleN),
                         dimnames = list(dimnames(psis[[chr]])[[2]], sampleN))
    for(samp in sampleN){
      vits[[chr]][,samp] = newvit[[chr]][[samp]]$viterbipath
    }

  }
  # loglike = 0 ## TODO!! IMPLEMENT LOGLIKE
  Viterbis = do.call(rbind, vits)

  DataViterbi_GR = lapply(colnames(Viterbis), function(x){
    pos = unlist(strsplit(rownames(Viterbis), split = "_"))
    mat = t(do.call(cbind, goodObs[[x]]))
    res = GRanges(seqnames = pos[seq(1,length(pos), by = 2)],
                  IRanges(start = as.numeric(pos[seq(2,length(pos), by = 2)]), width = tiles),
                  P1.Allele.Count = mat[,1],
                  P2.Allele.Count = mat[,2] - mat[,1],
                  Viterbi = Viterbis[,x],
                  # PostDec = postDec[,x],
                  seqlengths = anno
                  )
    return(res)


  })
  names(DataViterbi_GR) = colnames(Viterbis)



  myObj = .RViterbi(Rparameter = rigidity,
                    LogLike = loglike,
                    Viterbi = DataViterbi_GR)

  return(myObj)


}


summary.states = function(object,
                          main = NULL
                          # probs = seq(0,1,.1),
                          # breaks = 20, xlim = NULL
                          ){
  DataViterbi_GR = object@Viterbi

  state.lengths = unlist(lapply(DataViterbi_GR, function(samp){
    chr = lapply(seqlevels(samp), function(cr){
      FinalRes_chr = samp[seqnames(samp) == cr]
      Viterbi = rle(FinalRes_chr$Viterbi)$lengths
      return(Viterbi)
    })
    names(chr) = seqlevels(samp)
    return(chr)
  }))
  state.lengths = as.vector(state.lengths)

  myA = table(state.lengths)
  myA = myA/sum(myA) * 100
  x = log2(as.numeric(names(myA)))

  plot(x, myA, type = "h", xlab = "log2 state nucleotide length", main = main)
  # if(!is.null(xlim)){
  #     hist(state.lengths,
  #                   probability = TRUE,
  #                   breaks = breaks,
  #                   xlim = xlim
  # ) } else
  # {
  # hist(state.lengths,
  #         probability = TRUE,
  #         breaks = breaks)
  # }


  # quantile(state.lengths, probs)
}

#' Prints description of RViterbi object
#'
#' @keywords internal
#' @noRd
setMethod(f = "show", signature = c("RViterbi"), function(object) {


  cat(" An object of class \"", class(object), "\" \n", sep = "")
  cat(" R-Parameter:  ", object@Rparameter, "\n", sep = "")
  cat(" LogLikelihood: ", object@LogLike, "\n",
      sep = "")

})

