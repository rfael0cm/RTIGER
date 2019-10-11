###
# Genotyping HMM part1 data generation.R
###

# contains the functions:
#   generate_params
#   hidden_to_emission
#   generate_Markovchain 
#   generate_data
#   sanity_check

### function that generates a set of parameters
# params is a list containing the entries 
#$a (the transition matrix), $pim (the marker probabilities)
# $piv (the initial hidden genomic state distribution) and $alpha_s, 
# $beta_s (the parameters of the emission distributions)
# randomize tells whether some random variation should be added to the standard parameters
# generate_params = function(params=list(),randomize = FALSE){
#   vstates = c("pat","het","mat")
#   mstates = c("+","badP","badM")
#   estates = c("pat","het","mat","badP","badM")
#   
#   # transition matrix
#   if (is.null(params$a)){
#     # transition matrix a
#     vstates = c("pat","het","mat")
#     a = matrix(0.1,ncol=3,nrow=3) + diag(3)*10
#     # if (randomize) a = a + matrix(runif(9,min = 0, max = 0.5),ncol=3)
#     if (randomize) a = a + matrix(runif(9,min = 0, max = 0.01),ncol=3)
#     a = a / rowSums(a)
#     dimnames(a) = list(vstates,vstates)
#     a["pat", "mat"] = a["mat", "pat"] = 0
#     a = a/rowSums(a)
#     a["het", "het"] = mean(a["pat", "pat"], a["mat", "mat"])
#     a = a/rowSums(a)
#     params$a = a
#     }
#   
#   # beta-binomial emission distributions
#   if (is.null(params$alpha)){
#     psi = vector("list", length(estates))
#     priorstrength = 1 + randomize*runif(1,min = -0.5, max = 0.5)
#     alphavec = c(20,20,1,20,1) + randomize*runif(5,min = -.5, max = .5)
#     alpha_s = alphavec * priorstrength
#     names(alpha_s) = estates
#     betavec = c(1,20,20,1,20) + randomize*runif(5,min = -.5, max = .5)
#     beta_s = betavec * priorstrength
#     names(beta_s) = estates
#     
#     params$alpha_s = alpha_s
#     params$beta_s = beta_s
#     }
# 
#   # starting distribution for hidden states
#   if (is.null(params$piv)){
#     piv = c(1,1,1) + randomize * runif(3,min=-1,max=3)
#     piv = piv / sum(piv)
#     names(piv) = vstates
#     params$piv = piv
#   }
# 
#   # probabilities for the marker states
#   if (is.null(params$pim)){
#     pim = c(8,1,1) + randomize * runif(3,min=-1,max=3)
#     pim = pim / sum(pim)
#     names(pim) = mstates
#     params$pim = pim
#   }
#   
#   return(params)
# } # end generate_params


### The data: observations, info
### observations: A named list, each entry corresponding to one sample.
# Each sample is another list, consisting of chromosomes
# (or consecutive parts, can be analyzed in parallel)
# each chromosome/part must have a unique name, in order to match parts across samples.
# The actual data (for each sample, each part) consists of a 2 x n matrix, 
# with paternal / total counts as rows, and marker positions as columns
# The rows names are "maternal" and "total"
# The columns are named with their unique marker name
# identical markers are matched across samples during the analysis process
### info: A list containing $sample_names, $part_names, $marker_names
# the last one being again a named list (names = part_names), 
# each entry contains the unique marker names for the respective part


### the hidden_toemission function maps hidden states and marker states to an emission state
hidden_to_emission = function(v,m){
  if (m=="+") v else rep(m,length(v))  
}


### function that generates data (for one part of one sample)
# from a Markov chain of length n, given the parameters (params) and the marker names of the part
# coverage is a function that (either deterministically or randomly) 
# creates the total number of reads at each marker position
generate_Markovchain = function(marker_states,params, r=1,coverage=function(){sample(5:50,1)}){
  n = length(marker_states)
  if (n<2) stop("Markov chain too short. Choose larger n")
  if (sum(marker_states == "+") < 2) stop("Too few + markers in the chain")
  if (length(names(marker_states))==0) stop("Marker states have no marker names")
  vstates = names(params$piv)
  hidden = numeric(n)
  names(hidden) = names(marker_states)
  hidden[1] = sample(vstates,1,prob=params$piv)
  myr = 1
  for (j in 2:n){
    if(myr <= r){
      hidden[j] = hidden[j-1]
      myr = myr+1
    }
    else{
      hidden[j] = ifelse(marker_states[j]=="+",
                         sample(vstates,1,prob=params$a[hidden[j-1],]),
                         hidden[j-1])
      if(hidden[j] == hidden[j-1]) r = r+1
      else{
        myr = 1
      }
      
    }
    
  }
  nreads = coverage(n)
  obs = sapply(1:n,function(x){
    emi = hidden_to_emission(hidden[x],marker_states[x])
    rbb(1,nreads[x],params$alpha_s[emi],params$beta_s[emi])
  })
  obs = rbind(obs,nreads)
  rownames(obs) = c("paternal","total")
  colnames(obs) = names(marker_states)
  res = list(obs = obs, hidden = hidden)
  return(res)
}


### function that generates artificial data, given the parameters
# and either given an info object (containing sample_names,part_names,marker_names as specified above)
# or given the number of samples, parts per sample, and the respective part lengths
generate_data = function(params=NULL, info=NULL, r = 1,
                         sample_nr = NULL, part_nr = NULL, part_lengths = NULL, 
                         coverage=function(){sample(5:50,1)}){
  # generation of an the info object
  if (is.null(params)) params = generate_params() 
  if (is.null(info)){
    if (is.null(sample_nr)) stop("sample_nr must be specified.")
    sample_names = paste("sample",1:sample_nr,sep="_")
    if (is.null(part_nr)) stop("part_nr must be specified.")
    part_names = paste("chr",1:part_nr,sep="")
    if (is.null(part_lengths)) stop("part_lengths must be specified.")
    if (length(part_lengths) == 1) part_lengths = rep(part_lengths,part_nr)
    if (length(part_lengths) != part_nr) stop("length of part_lengths does not match part_nr.")
    names(part_lengths) = part_names
    marker_names = vector("list",length(part_names))
    names(marker_names) = part_names
    for (part in part_names){
      marker_names[[part]]=paste(part,1:part_lengths[[part]],sep="_M")
    }
  } else {
    sample_names = info$sample_names
    sample_nr = length(sample_names)
    part_names = info$part_names
    part_nr = length(part_names)
    part_lengths = sapply(info$marker_names,length)
    marker_names = info$marker_names
  }
  
  info = list(sample_names = sample_names, part_names = part_names, marker_names = marker_names,
              sample_nr = sample_nr, part_nr = part_nr, part_lengths = part_lengths)
  
  # drawing of marker states according to pim, for each part
  marker_states = vector("list",length(part_names))
  names(marker_states) = part_names
  for (part in part_names){
    marker_states[[part]] = sample(c("+","badP","badM"),part_lengths[[part]],prob=params$pim,replace = TRUE)
    names(marker_states[[part]]) = marker_names[[part]]
  }
  
  # drawing of observations according to the marker states and the Markov chain parameters,
  # for each sample and each part
  # the "true" hidden states are recorded in the truth object
  obs = vector("list",sample_nr)
  names(obs) = sample_names
  truth = vector("list",part_nr)
  names(truth) = part_names
  for (part in part_names){
    truth[[part]] = matrix(0,nrow=sample_nr,ncol=part_lengths[part])
    dimnames(truth[[part]]) = list(sample_names,marker_names[[part]])
  }
  for (samp in sample_names){
    obs[[samp]] = vector("list",part_nr)
    names(obs[[samp]]) = part_names
    for (part in part_names){
      hilf = generate_Markovchain(marker_states[[part]], r = r, params,coverage=coverage)
      obs[[samp]][[part]] = hilf$obs
      truth[[part]][samp,] = hilf$hidden
    }
  }
  truth = c(truth,marker_states = list(marker_states))
  
  return(list(obs=obs, info=info, truth=truth))
} # end generate_data

# check whether the data has the correct format
sanity_check = function(info,obs){
  sample_names = info$sample_names
  sample_nr = length(sample_names)
  if (sample_nr<2) stop("Need at least 2 samples.")
  if (sample_nr==0) stop("Sample names are missing.")
  rest = setdiff(sample_names,names(obs))
  if (length(rest)>0) stop(paste("The sample(s) (",rest,") are missing.",collapse=""))
  
  part_names = info$part_names
  part_nr = length(part_names)
  if (part_nr==0) stop("Chromosome/parts names are missing.")
  
  marker_names = info$marker_names
  if (length(unique(marker_names)) < length(marker_names)) stop("Marker names are not unique.")
  if ( !identical(sort(names(marker_names)), sort(part_names)) ) stop("marker list does not fit parts names.")
  
  test = sapply(sample_names,function(x){all(names(obs[[x]])==part_names)})
  if (!all(test==TRUE)) stop("Chromosome/parts names are inconsistent across samples.")
  part_lengths = sapply(info$marker_names,function(x){length(x)})

  for (samp in sample_names){
    if (!all(part_names %in% names(obs[[samp]]))) stop("Part names are missing.") 
    for (part in part_names){
      if (!all(marker_names[[part]] %in% colnames(obs[[samp]][[part]]))) stop("Marker names are missing.")
      if (!all(rownames(obs[[samp]][[part]]) == c("paternal","total"))) stop("Row names are not (paternal,total).")
    }
  }
 
  return(TRUE)
} # end sanity_check

