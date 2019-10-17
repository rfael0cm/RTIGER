prescreen = function(gammas,    # posterior state probability matrix (states x timepoints)
                     threshold = 0.5, # probability threshold for the gammas, must be >0.5
                     maxwidth = 3  # maximum width of the interval
){
  rt = apply(gammas,2,function(x){
    y = which(x>threshold)
    ifelse(length(y)==0,0,y)} )
  if (all(rt==0)) stop("no candidates found")
  rt = cbind(rt,1:length(rt))
  st = rt[which(rt[,1]>0),,drop=F]
  selected = which(diff(st[,1])!=0)
  changepoints = st[selected,2]
  leftstate = st[selected,1]
  rightstate = st[selected+1,1]
  
  res = sapply(1:length(changepoints),function(x){
    position = changepoints[x]
    fromstate = leftstate[x]
    tostate = rightstate[x]
    leftselect = rt[max(position-maxwidth+1,1):position,1] == fromstate
    leftends = (max(position-maxwidth+1,1):position)[leftselect]
    maxpos = ncol(gammas)
    rightselect = rt[min(maxpos,position+1):min(position+maxwidth,maxpos),1] == tostate
    rightends = (min(maxpos,position+1):min(position+maxwidth,maxpos))[rightselect]
    mat = outer(rightends,leftends,"-")
    indices = which(maxwidth >= mat,arr.ind=T)
    items = nrow(indices)
    if (items==0) return(NULL)
    intervals = cbind(rep(fromstate,items),rep(tostate,items),
                      leftends[indices[,2]],rightends[indices[,1]])
    colnames(intervals) = c("left_state","righ_tstate","left_end","right_end")
    return(t(intervals))
  },simplify=T)
  
  res = matrix(unlist(res),nrow=4)
  res = t(unique(t(res)))
  rownames(res) = c("left_state","righ_tstate","left_end","right_end")
  return(res) 
}


alpha_beta_gamma_psi = function(mat,  params, psis = NULL){
  
  transitions = params$a
  initiation = params$piv
  vstates = colnames(transitions)
  timepoints = ncol(mat)
  sizes = colSums(mat)
  psivec = function(x){
    sapply(vstates,function(emissionstate){
      if( any(is.na(x))) r = 1
      else{
        r = dbb(x[1],x[2],
                params$alpha_s[emissionstate],params$beta_s[emissionstate])
      }
      return(r)
    })
  }
  psis = apply(mat,2,psivec)
  rownames(psis) = vstates
  
  # calculation of the forward probabilities
  alphas = matrix(NA, nrow = 3, ncol = timepoints, dimnames = list(rownames(transitions),1:timepoints))
  normfactors = numeric(timepoints)
  
  normfactors[1] = sum(initiation*psis[,1])
  alphas[,1] = initiation*psis[,1]/normfactors[1] 
  
  sapply(2:timepoints,function(x){
    hilf = t((t(alphas[,x-1]) %*% transitions)) * psis[,x]
    normfactors[x] <<- sum(hilf)
    alphas[,x] <<- hilf / normfactors[x]
  })
  
  # calculation of the backward probabilities
  betas = matrix(NA, nrow = 3, ncol = timepoints, dimnames = list(rownames(transitions),1:timepoints))
  betas[,timepoints] = rep(1,3)
  
  sapply((timepoints-1):1,function(x){
    hilf = transitions %*% (psis[,x+1] * betas[,x+1])
    betas[,x] <<- hilf / sum(hilf)
  })
  
  gammas = alphas * betas
  gammas = t(t(gammas)/colSums(gammas))
  
  return(list(alphas=alphas,betas=betas,gammas=gammas,psis=psis,normfactors=normfactors))
}

intervalTransition <- function(a, # Left boundary of the interval
                               b, # Right boundary of the interval
                               sj, # First state from transition
                               sk,  # Second state of the transition
                               alphas, # Normalized forward probability with normalization Factors
                               betas,  # Normalized backward probability with Normalization Factors
                               normfactors, # scaling factors for the alphas
                               transitions, # a_jk transition probs for moving from j to k
                               psis # Emission Probabilities
){
  transition = numeric(3)
  transition[1] = transitions[sj,sj] 
  transition[2] = transitions[sj,sk]
  transition[3] = transitions[sk,sk]
  
  # fill the transition probability matrix
  d = b - a 
  transitionnumber_mat = diag(2,d) + upper.tri(diag(d))*3 + lower.tri(diag(d))*1
  amat = matrix(transition[transitionnumber_mat],nrow=nrow(transitionnumber_mat))
  
  # fill the emission probability matrix
  psimat = sapply(1:d,function(x){ifelse(transitionnumber_mat[,x]==1,psis[sj,a+x], psis[sk,a+x])})
  psimat = t(t(psimat)/normfactors[(a+1):b]) ## DIVIDE psi by normFactors - yes, that saves time to do it here
  
  # calculate the interval transition probability
  sumofpaths = ifelse((sj!=sk), sum(apply(amat * psimat,1,prod)), prod(psimat[1,]*amat[1,]) )
  if (b<ncol(psis)){
    nominator = alphas[sj,a]*betas[sk,b] * sumofpaths
    denominator = sum(alphas[,b]*betas[,b])
  } else {
    nominator = alphas[sj,a] * sumofpaths
    denominator = sum(alphas[,b])
  }
  res = nominator/denominator
  
  return(res)
}

wigglyTransition <- function(a, # Left boundary of the interval
                             b, # Right boundary of the interval
                             sj, # First state from transition
                             sk,  # Second state of the transition
                             alphas, # Normalized forward probability with normalization Factors
                             betas,  # Normalized backward probability with Normalization Factors
                             normfactors, # scaling factors for the alphas
                             transitions, # a_jk transition probs for moving from j to k
                             psis # Emission Probabilities
){
  d = b - a + 1
  gtilde = matrix(NA,ncol=2,nrow=d)
  colnames(gtilde)=c("sj","sk")
  gtilde[1,]=c(1,0)
  transreduced = transitions[c(sj,sk),c(sj,sk)]
  if (sj != sk){
    sapply(2:d,function(x){
      gtilde[x,] <<- t(gtilde[x-1,,drop=F] %*% transreduced) * psis[c(sj,sk),a+x-1] / normfactors[a+x-1]
    })
    gende = gtilde[d,"sk"]
  }
  if (sj == sk){
    gende = transitions[sj,sj]^(b-a) * prod( psis[sj,(a+1):b] * normfactors[(a+1):b])  
  }
  
  # calculate the interval transition probability
  if (b<ncol(psis)){
    nominator = alphas[sj,a]*betas[sk,b] * gende
    denominator = sum(alphas[,b]*betas[,b])
  } else {
    nominator = alphas[sj,a] * gende
    denominator = sum(alphas[,b])
  }
  res = nominator/denominator
  
  return(res)
}


