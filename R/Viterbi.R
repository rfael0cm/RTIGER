Viterbi_decode = function(obs, params, psis = NULL){
  if(is.null(psis)){
    vstates = c("pat","het","mat")
    psivec = function(x){
      sapply(vstates,function(emissionstate){
        if( any(is.na(x))) r = 1
        else{
          r = TailRank::dbb(x[1],x[2],
                            params$alpha_s[emissionstate],params$beta_s[emissionstate])
        }
        return(r)
      })
    }
  }
   
  transMat = params$a
  ip = params$piv
  vstates = names(ip)
  nStates = length(vstates)

  nsamples = names(obs)
  nparts = names(obs[[1]])

  p = lapply(params$alpha_s/(params$alpha_s + params$beta_s), identity)
  p = p[vstates]

  Emis = STAN::HMMEmission(type = "Bernoulli",
                     parameters = list(p = p),
                     nStates = nStates)

  hmm = STAN::HMM(initProb = ip,
            transMat = transMat,
            emission = Emis,
            nStates = nStates)

  myVit = vector("list", length(nparts))
  names(myVit) = nparts
  for(part in nparts){
    nMarkers = colnames(obs[[1]][[part]])
    if(is.null(nMarkers)) nMarkers = paste("pos", 1:ncol(obs[[1]][[part]]), sep ="_")
    myVit[[part]] = matrix(NA, ncol = length(nMarkers), nrow = length(nsamples),
                           dimnames = list( nsamples, nMarkers))
    for(samp in nsamples){
      
      if(is.null(psis)){
        psimat = apply(obs[[samp]][[part]],2,psivec)
        rownames(psimat) = vstates
        psimat = t(psimat)
      }
      
    
      op = t(obs[[samp]][[part]])
      op =as.matrix(as.numeric(op[,2]))

      if(!is.null(psis)) psimat = psis[[part]][samp,,vstates]
      STANVit = colnames(transMat)[unlist(STAN::getViterbi(hmm = hmm,  emissionProbs = list(psimat), obs = list(op)))]
      myVit[[part]][samp,] = STANVit
    }
  }

  return(myVit)


}
