fit = function(rtigerobj, max.iter, eps, trace, all = TRUE, random = FALSE, specific = FALSE, nsamples = 20, post.processing = TRUE){
  params = rtigerobj@params
  obs = rtigerobj@matobs
  info = rtigerobj@info

  obs = lapply(obs, function(samp){
    chr = lapply(samp, function(cr){
      return(t(cr))
    })
    return(chr)
  })
  cat("Inside fit the postprocessing value is:", post.processing, "\n")
  # function fit(Observations,info,initial_parameter,max_iter=100,eps=10^(-5),trace=false)
  myfit = julia_call("fit",obs, info, params, as.integer(max.iter), eps, trace, all, random , as.integer(nsamples), specific, post.processing )
  # function fit(
  #   input_Observations,
  #   info,
  #   initial_parameter,
  #   max_iter = 100,
  #   eps = 10^(-5),
  #   trace = false,
  #   all = true,
  #   random = true,
  #   nsamples=20,
  #   specific = nothing,
  # )
  nstates = myfit$parameterSet$nstates

  myvit = myfit$viterbiPath
  myvit = myvit[names(obs)]
  myvit = lapply(myvit, function(x) x[info$part_names])
  myvit = lapply(myvit, unlist)
  if(nstates == 3){
    vits = c( "pat", "het", "mat")
  } else if(nstates == 2){
    vits = c( 0,0)
    a = myfit$parameterSet$paraBetaAlpha
    b = myfit$parameterSet$paraBetaBeta
    par.dif = as.vector(a-b)
    hetst = which.min(abs(par.dif))
    vits[hetst] = "het"
    vits[-hetst] = ifelse(par.dif[-hetst] > 0, "pat", "mat")
  } else{
    vits = 1:nstates
  }


  for(i in info$sample_names){
    rtigerobj@Viterbi[[i]]$Viterbi = vits[myvit[[i]]]
  }
  rtigerobj@params = myfit$parameterSet
  rownames(rtigerobj@params$logtransition) = colnames(rtigerobj@params$logtransition) = vits
  rownames(rtigerobj@params$paraBetaAlpha) = vits
  rownames(rtigerobj@params$paraBetaBeta) = vits
  rownames(rtigerobj@params$logpi) = vits
  rownames(rtigerobj@params$pi) = vits
  rownames(rtigerobj@params$transition) = colnames(rtigerobj@params$transition) = vits

  rtigerobj@Probabilities = myfit[c("alpha", "beta", "gamma", "psi")]
  rtigerobj@num.iter = myfit$numberofiterations
  return(rtigerobj)

}
