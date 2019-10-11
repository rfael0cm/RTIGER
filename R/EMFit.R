###########################################################
##### Functions for the EM algorithm
##### Final function fitModel(object)
##### Intermediate steps
##############################################################

# contains the functions:
#   forwardbackward
#   wrap_fwbw
#   posterior_markers
#   zeta_function
#   update_transitions
#   tau_function
#   update_emissions
#   update_piv
#   update_pim
#   hidden_function
#   EMalgorithm

#' Compute the forward and backward probabilities
#' @param mat the observation matrix.
#' @param params the parameters to compute the forward and backward
#'
#' @keywords internal
#' @noRd
#'


forwardbackward = function(mat,params){
  vstates = c("pat","het","mat")
  marker_nr = ncol(mat)

  # psi_psistate(O^t) matrix
  psistates = c("pat","het","mat","badP","badM")
  psivec = function(x){
    sapply(psistates,function(emissionstate){
      if( any(is.na(x))) r = 1
      else{
        r = TailRank::dbb(x[1],x[2],
                params$alpha_s[emissionstate],params$beta_s[emissionstate])
      }
      return(r)
    })
  }
  
  psimat = apply(mat,2,psivec)
  
  rownames(psimat) = psistates

  # alphas
  alpha = matrix(0,ncol=marker_nr,nrow=3,
                 dimnames = list(vstates,colnames(mat)))
  nalpha = numeric(marker_nr)
  pim = params$pim
  a = params$a
  alpha[,1] = params$piv * (psimat[vstates,1]*pim["+"] +
                              rep(psimat["badP",1]*pim["badP"] + psimat["badM",1]*pim["badM"],3) )
  nalpha[1] = sum(alpha[,1], na.rm = TRUE)
  alpha[,1] = alpha[,1]/nalpha[1]
  sapply(2:marker_nr,function(t){
    hilf = alpha[,t-1] *
      (psimat["badP",t]*pim["badP"] + psimat["badM",t]*pim["badM"]) +
      psimat[vstates,t]*pim["+"]*(a%*% alpha[,t-1])
    nalpha[t] <<- sum(hilf)
    alpha[,t] <<- hilf/nalpha[t]
  })
  #betas
  beta = matrix(0,ncol=marker_nr,nrow=3,
                dimnames = list(vstates,colnames(mat)))
  nbeta = numeric(marker_nr)
  beta[,marker_nr] = rep(1,3)
  nbeta[marker_nr] = 1
  sapply((marker_nr):2,function(t){
    hilf = (pim["badP"]*psimat["badP",t]+pim["badM"]*psimat["badM",t]) +
      pim["+"]*psimat[vstates,t]*(a%*%beta[,t])
    nbeta[t-1] <<- sum(hilf, na.rm = TRUE)
    beta[,t-1] <<- hilf/nbeta[t-1]
  })
  return(list(alpha=alpha,beta=beta,nalpha=nalpha,nbeta=nbeta,psimat=psimat))
} # end forwardbackward

#' Generates the forward-backward table (for all samples, all parts)
#' @param obs the list of samples with the list of parts observations
#' @param info the information object containing NAs, number samples, etc
#' @param params the parameters to compute the forward and backward
#'
#' @keywords internal
#' @noRd
#'

wrap_fwbw = function(obs,info,params){
  sample_names = info$sample_names
  sample_nr = length(sample_names)
  part_names = info$part_names
  part_nr = length(part_names)

  vstates = c("pat","het","mat")
  psistates = c("pat","het","mat","badP","badM")
  alpha_complete = vector("list",part_nr)
  names(alpha_complete) = part_names
  beta_complete = vector("list",part_nr)
  names(beta_complete) = part_names
  psimat_complete = vector("list",part_nr)
  names(psimat_complete) = part_names
  k_obs = vector("list",part_nr)
  names(k_obs) = part_names
  n_obs = vector("list",part_nr)
  names(n_obs) = part_names
  nalphas = vector("list",part_nr)
  names(nalphas) = part_names

  nbetas = vector("list",part_nr)
  names(nbetas) = part_names

  for (part in part_names){
    marker_names = info$marker_names[[part]]
    marker_nr = length(marker_names)
    alpha_complete[[part]] = array(0,dim=c(sample_nr,ncol=marker_nr,length(vstates)),
                                   dimnames=list(sample_names,marker_names,vstates))
    beta_complete[[part]] = array(0,dim=c(sample_nr,ncol=marker_nr,length(vstates)),
                                  dimnames=list(sample_names,marker_names,vstates))
    psimat_complete[[part]] = array(0,dim=c(sample_nr,ncol=marker_nr,length(psistates)),
                                    dimnames=list(sample_names,marker_names,psistates))
    k_obs[[part]] = array(0,dim=c(sample_nr,ncol=marker_nr),
                          dimnames=list(sample_names,marker_names))
    n_obs[[part]] = array(0,dim=c(sample_nr,ncol=marker_nr),
                          dimnames=list(sample_names,marker_names))
    nalphas[[part]] = matrix(0, nrow = sample_nr, ncol = marker_nr,
                             dimnames = list( sample_names, marker_names))
    nbetas[[part]] = matrix(0, nrow = sample_nr, ncol = marker_nr,
                            dimnames = list( sample_names, marker_names))

    for (samp in sample_names){
      # cat(samp, " ", part, "\n")
      res = forwardbackward(obs[[samp]][[part]],params)
      alpha_complete[[part]][samp,marker_names,vstates] = t(res$alpha[vstates,marker_names])
      beta_complete[[part]][samp,marker_names,vstates] = t(res$beta[vstates,marker_names])
      psimat_complete[[part]][samp,marker_names,psistates] = t(res$psimat[psistates,marker_names])
      k_obs[[part]][samp,] = obs[[samp]][[part]]["paternal",]
      n_obs[[part]][samp,] = obs[[samp]][[part]]["total",]
      nalphas[[part]][samp,] = res$nalpha
      nbetas[[part]][samp,] = res$nbeta
    }
  }
  return(list(alpha = alpha_complete,
              beta = beta_complete,
              psimat = psimat_complete,
              k_obs = k_obs,
              n_obs = n_obs,
              # sample_names = info$sample_names,
              # part_names=info$part_names,
              # marker_names = info$marker_names,
              nalpha = nalphas,
              nbeta = nbetas))
} # end wrap_fwbw

#' Compute posterior probability for the markers
#'
#' @param alpha Alpha array generated from the wrap_fwbw function
#' @param beta Beta array generated from the wrap_fwbw function
#' @param psimat psi matrix generated form wrap_fwbw function
#' @param info the information object containing NAs, number samples, etc
#' @param params the parameters to compute the forward and backward
#'
#' @keywords internal
#' @noRd
#'

posterior_markers = function(alpha,beta,psimat,info,params){
  vstates = c("pat","het","mat")
  markerstates = c("+","badP","badM")

  sample_names = info$sample_names
  sample_nr = length(sample_names)
  part_names = info$part_names
  part_nr = length(part_names)

  vstates = c("pat","het","mat")
  piv = params$piv[vstates]
  pim = params$pim
  a = params$a
  gamma_complete = vector("list",part_nr)
  names(gamma_complete) = part_names
  marker_posterior = vector("list",part_nr)
  names(marker_posterior) = part_names

  for (part in part_names){
    marker_names = info$marker_names[[part]]
    marker_nr = length(marker_names)
    gamma_complete[[part]] = array(0,dim=c(samples=sample_nr,markers=marker_nr,
                                           markerstates=length(markerstates)),dimnames=list(sample_names,marker_names,markerstates))
    for (samp in sample_names){
      for (marker in markerstates){
        # gamma_1
        betas = beta[[part]][samp,1,vstates]
        emis = if (marker=="+") vstates else rep(marker,3)
        psis = psimat[[part]][samp,1,emis]
        gamma_complete[[part]][samp,1,marker] = sum(betas*psis*piv)

        # gamma_t, t>1
        if (marker=="+"){
          betas = beta[[part]][samp,-1,vstates]
          psis = psimat[[part]][samp,-1,vstates]
          alphas = alpha[[part]][samp,-marker_nr,vstates]
          alphas = alphas %*% t(a)
          gamma_complete[[part]][samp,2:marker_nr,marker] = rowSums(betas * psis * alphas)
        } else {
          betas = beta[[part]][samp,-1,vstates]
          psis = psimat[[part]][samp,-1,rep(marker,3)]
          alphas = alpha[[part]][samp,-marker_nr,vstates]
          gamma_complete[[part]][samp,2:marker_nr,marker] = rowSums(betas * psis * alphas)
        } # end if else

      } # end for marker
    } # end for samp

    for(samp in sample_names){
      gamma_complete[[part]][samp, , ] = gamma_complete[[part]][samp, , ]*info$NAs_matrices[[samp]][[part]][1,]
    }

    hilf = t(t(apply(log(gamma_complete[[part]]),2:3,sum, na.rm = TRUE)) + log(pim)) # implement na.rm = T
    hilf = hilf - apply(hilf,1,max)
    hilf = exp(hilf)
    hilf = hilf/rowSums(hilf)
    marker_posterior[[part]] = hilf
  } # end for part

  marker_annotation = vector("list",length(info$part_names))
  names(marker_annotation) = info$part_names
  for (part in info$part_names){
    hilf = apply(marker_posterior[[part]],1,which.max)
    marker_annotation[[part]] = colnames(marker_posterior[[part]])[hilf]
  }

  return(list(mu = marker_posterior,gamma = gamma_complete,marker_annotation = marker_annotation))
} # end posterior_markers


# Maximization Step -------------------------------------------------------

#' Compute the zeta parameter to learn the transitions
#'
#' @param alpha Alpha array generated from the wrap_fwbw function
#' @param beta Beta array generated from the wrap_fwbw function
#' @param psimat psi matrix generated form wrap_fwbw function
#' @param info the information object containing NAs, number samples, etc
#' @param params the parameters to compute the forward and backward
#'
#' @keywords internal
#' @noRd
#'

zeta_function = function(alpha,beta,psimat,info,params){
  vstates = c("pat","het","mat")
  badstates = c("badM","badP")
  nr_states = length(vstates)
  part_names = info$part_names
  nr_parts = length(part_names)
  sample_names = info$sample_names
  nr_samples = length(sample_names)
  zeta = vector("list",nr_parts)
  names(zeta) = part_names
  a = params$a

  for (part in part_names){
    nr_markers = length(info$marker_names[[part]])
    zeta[[part]] = array(0,dim=c(nr_samples,nr_markers-1,nr_states,nr_states))
    dimnames(zeta[[part]]) = list(info$sample_names,info$marker_names[[part]][-nr_markers],vstates,vstates)
    # zeta here equals zeta_t^j(r,s) in the lyx, for one part. dim = (j,t,r,s)
    for (s_state in vstates){
      hilf2 = psimat[[part]][,-1,s_state] # psi_t
      hilf3 = beta[[part]][,-1,s_state] # beta_t
      for( r_state in vstates){
        hilf1 = alpha[[part]][,-nr_markers,r_state] # alpha_t-1
        zeta[[part]][,,r_state,s_state] = hilf1 * a[r_state,s_state] * hilf2 * hilf3
      } # end for r_state
    } # end for s_state

    normfactor = apply(zeta[[part]],1:2,sum)
    for (r_state in vstates){
      for (s_state in vstates){
        zeta[[part]][,,r_state,s_state] = zeta[[part]][,,r_state,s_state] / normfactor
        # zeta[[part]] now equals zeta_t^j(r,s) from the lyx (for one part); dim = (j,t,r,s)
      } # end for s_state
    } # end for r_state
  } # end for part
  return(list(zeta=zeta))
} # end zeta_function

#' Update the transition matrix
#'
#' @param zeta the zeta array generated from the zeta function
#' @param mu the mu array computed in posterior_markers
#' @param info the information object containing NAs, number samples, etc
#'
#' @keywords internal
#' @noRd
#'

update_transitions = function(zeta,mu,info){
  vstates = c("pat","het","mat")
  nr_states = length(vstates)
  badstates = c("badM","badP")
  part_names = info$part_names
  nr_parts = length(part_names)

  a_new = matrix(0,ncol=nr_states,nrow=nr_states)
  dimnames(a_new) = list(vstates,vstates)

  for (part in part_names){
    hilf = apply(zeta[[part]],2:4,sum)
    hilf = as.vector(hilf) * mu[[part]][-1,"+"]
    dim(hilf)=c(dim(zeta[[part]])[2],nr_states,nr_states)
    # dim(zeta[[part]])[2] contains the number of markers in part, -1
    a_new = a_new + apply(hilf,2:3,sum)
  } # end for part

  # a_new["pat", "mat"] = a_new["mat", "pat"] = 0
  a_new = a_new/rowSums(a_new)
  return(a_new)
} # end update_transitions

#' Compute the tau probability to learn the emission distributions
#'
#' @param zeta Zeta prbability array generated from the zeta_function
#' @param beta Beta array generated from the wrap_fwbw function
#' @param psimat psi matrix generated form wrap_fwbw function
#' @param info the information object containing NAs, number samples, etc
#' @param params the parameters to compute the forward and backward
#'
#' @keywords internal
#' @noRd
#'


tau_function = function(zeta,beta,psimat,mu,info,params){
  vstates = c("pat","het","mat")
  badstates = c("badM","badP")
  estates = c(vstates,badstates)
  nr_estates = length(estates)
  part_names = info$part_names
  nr_parts = length(part_names)
  sample_names = info$sample_names
  nr_samples = length(sample_names)
  tau = vector("list",nr_parts)
  names(tau) = part_names

  for (part in part_names){
    marker_names = info$marker_names[[part]]
    nr_markers = length(marker_names)
    tau[[part]] = array(0,dim=c(nr_samples,nr_markers,nr_estates))
    dimnames(tau[[part]]) = list(info$sample_names,marker_names,estates)
    # tau here equals tau_t^j(s) in the lyx, for one part. dim = (j,t,s)
    # update vstates
    for (e_state in vstates){
      tau[[part]][,,e_state] = matrix(rep(mu[[part]][,"+"],nr_samples),nrow=nr_samples,byrow= TRUE)
    }
    tau[[part]][,2:nr_markers,vstates] = tau[[part]][,2:nr_markers,vstates] * apply(zeta[[part]],c(1,2,4),sum)
    hilf = t( t(psimat[[part]][,1,vstates] * beta[[part]][,1,vstates]) * params$piv[vstates] )
    tau[[part]][,1,vstates] = tau[[part]][,1,vstates] * hilf / rowSums(hilf)
    # update badstates
    for (e_state in badstates){
      tau[[part]][,,e_state] = matrix(rep(mu[[part]][,e_state],nr_samples),nrow=nr_samples,byrow = TRUE)
    }
  } # end for part
  # tau[[part]] here is tau_t^j(s) from the lyx (for one part); dim = (j,t,s)
  return(list(tau=tau))
} # end tau_function

#' Update the emission distributions
#'
#' @param tau Tau probability generated from the tau_function
#' @param k_obs allele count for allele of parent 1
#' @param n_obs total allele counts
#' @param info the information object containing NAs, number samples, etc
#' @param iteration the iteration step
#'
#' @keywords internal
#' @noRd
#'

update_emissions = function(tau,k_obs,n_obs,info, iteration){
  part_names = info$part_names
  nr_parts = length(part_names)
  vstates = c("pat","het","mat")
  badstates = c("badP","badM")
  estates = c(vstates,badstates)
  nr_estates = length(estates)

  expectation = numeric(nr_estates)
  names(expectation) = estates
  nominator = numeric(nr_parts)
  names(nominator) = part_names
  denominator = numeric(nr_parts)
  names(denominator) = part_names

  for (e_state in estates){
    for (part in part_names){
      nominator[part] = sum(tau[[part]][,,e_state] * k_obs[[part]], na.rm = TRUE)
      denominator[part] = sum(tau[[part]][,,e_state] * n_obs[[part]], na.rm = TRUE)
    }
    expectation[e_state] = sum(nominator) / sum(denominator)
  }

  # targetfunction calculates the relevant part of Q(Theta,Theta') for a given e_state
  targetfunction = function(testvalue,expected,e_state){
    # create psimat for the parameters alpha_e
    alpha_s = testvalue * expected
    beta_s = testvalue * (1-expected)
    hilf = sapply(part_names,function(part){
      psimat_new = TailRank::dbb(k_obs[[part]],n_obs[[part]],alpha_s,beta_s,log = TRUE)
      sum(tau[[part]][,,e_state]*psimat_new)
    }) # end sapply part
    return(sum(hilf))
  } # end targetfunction

  strength = numeric(nr_estates)
  names(strength) = estates
  alpha_new = numeric(nr_estates)
  names(alpha_new) = estates
  beta_new = numeric(nr_estates)
  names(beta_new) = estates

  for (e_state in estates){
    strength = optimize(f = targetfunction,
                        interval = c(1,1000),
                        expected=expectation[e_state],
                        e_state=e_state,
                        maximum = TRUE)$maximum
    alpha_new[e_state] = expectation[e_state] * strength
    beta_new[e_state] = (1-expectation[e_state]) * strength
  }

  return(list(alpha_new=alpha_new,beta_new=beta_new))
} # end update_emissions


#' Update the initial probability for HMM states
#'
#' @param beta Beta probabilities form the wrap_fwbw function
#' @param psimat emission probabilities from the wrap_fwbw function
#' @param mu posterior probabilites from the posterior HMM function
#' @param info the information object containing NAs, number samples, etc
#' @param params parameters of the model
#' @param initial_parts In case the chromosomes are split into smaller parts, it indicates which of the parts are considered as startin parts ofthe markov chain and have to be considered for piv
#'
#' @keywords internal
#' @noRd
#'

update_piv = function(beta,psimat,mu,info,params,initial_parts=NULL){
  sample_names = info$sample_names
  nr_samples = length(sample_names)
  part_names = info$part_names
  if (is.null(initial_parts)) initial_parts = part_names
  if (!all(initial_parts %in% part_names)) stop("Initial parts are contain undefined part names.")
  nr_parts = length(initial_parts)
  vstates = c("pat","het","mat")
  nr_vstates = length(vstates)
  badstates = c("badP","badM")
  mstates = c("+",badstates)
  nr_mstates = length(mstates)

  hilf = array(0,dim=c(nr_parts,nr_samples,nr_vstates,nr_mstates))
  dimnames(hilf) = list(initial_parts,sample_names,vstates,mstates)
  for (part in initial_parts){
    for (v_state in vstates){
      for (m_state in mstates){
        hilf[part,,v_state,m_state] =
          psimat[[part]][,1,hidden_to_emission(v_state,m_state)] *
          beta[[part]][,1,v_state] * rep(params$piv[v_state],nr_samples)
      }
    }
  }

  normfactor = apply(hilf,c(1,2,4),sum)
  for (part in initial_parts){
    for (v_state in vstates){
      for (m_state in mstates){
        hilf[part,,v_state,m_state] =
          hilf[part,,v_state,m_state] / normfactor[part,,m_state] * params$pim[m_state]
      }
    }
  }

  piv = apply(hilf,3,sum) / (nr_samples*nr_parts)
  return(list(piv=piv))
} #end update_piv


#' Update the state distribution of the good, badP, badM clusters
#'
#' @param mu posterior probabilites from the posterior HMM function
#'
#' @keywords internal
#' @noRd
#'

update_pim = function(mu){
  pim_new = sapply(mu,colSums)
  pim_new = rowSums(pim_new)
  pim_new = pim_new / sum(pim_new)
  return(list(pim_new=pim_new))
} # end update_pim


#'  the hidden_toemission function maps hidden states and marker states to an emission state
#'
#' @param v hidden states
#' @param m marker state
#' @keywords internal
#' @noRd
#'

###
hidden_to_emission = function(v,m){
  if (m=="+") v else rep(m,length(v))
}


#' Maximum posterior decoding for hidden states V_t^j
#'
#' @param alpha Alpha probabilites from the warp_fwbw function
#' @param beta Beta probabilities form the wrap_fwbw function
#' @param info the information object containing NAs, number samples, etc
#'
#' @keywords internal
#' @noRd
#'

hidden_function = function(alpha,beta,info){
  vstates = dimnames(alpha[[1]])[[3]]
  part_names = info$part_names
  hidden_posterior = vector("list",length(part_names))
  names(hidden_posterior) = part_names
  hidden_annotation = vector("list",length(part_names))
  names(hidden_annotation) = part_names
  for (part in info$part_names){
    hilf = alpha[[part]] * beta[[part]]
    hilf = as.vector(hilf) / rep(apply(hilf,1:2,sum),3)
    dim(hilf) = dim(alpha[[part]])
    dimnames(hilf) = dimnames(alpha)
    hidden_posterior[[part]] = hilf
    hilf2 = apply(hilf,1:2,which.max)
    hilf2 = vstates[hilf2]
    dim(hilf2) = dim(hilf)[1:2]
    dimnames(hilf2) = dimnames(alpha[[part]])[1:2]
    hidden_annotation[[part]] = hilf2
  }
  return(list(hidden_annotation = hidden_annotation,
              hidden_posterior = hidden_posterior))
} # hidden_function

# Final EM algorithm ------------------------------------------------------


#' EM algorithm
#'
#' @param obs list of list of chromosome matrices containing the paternal and total allele count
#' @param info the information object containing NAs, number samples, etc
#' @param initial_params list with initial parameters a (transition matrix), piv (itial state probailities) and pim (state distribution of markers)
#' @param ransomize Logical variable wether the parameters should be randomized if initial parameters are null.
#' @param eps Total difference of the parameters after each iteration
#' @param max.iter Maximum number of iterations before stopping the EM algorithm
#' @param trace Logical variable wether the updated parameters should be saved or not.
#'
#' @keywords internal
#' @noRd
#'


EMalgorithm = function(observations, info,
                       initial_params=NULL, randomize = FALSE,
                       eps=0.01, max.iter=50,
                       trace=FALSE,
                       verbose = FALSE){
  params = generate_params(initial_params,randomize)

  if (trace) {
    a_trace = matrix(0,nrow=max.iter+1,ncol=9)
    a_trace[1,] = params$a
    alpha_s_trace = matrix(0,nrow=max.iter+1,ncol=5)
    alpha_s_trace[1,] = params$alpha_s
    beta_s_trace = matrix(0,nrow=max.iter+1,ncol=5)
    beta_s_trace[1,] = params$beta_s
    piv_trace = matrix(0,nrow=max.iter+1,ncol=3)
    piv_trace[1,] = params$piv
    pim_trace = matrix(0,nrow=max.iter+1,ncol=3)
    pim_trace[1,] = params$pim
  }
  iteration = 1

  repeat{
    if(verbose) cat("Iteration: ", iteration, "\n")
    res = wrap_fwbw(observations,info,params)
    res = c(res,posterior_markers(res$alpha,res$beta,res$psimat,info,params))
    res = c(res,zeta_function(res$alpha,res$beta,res$psimat,info,params))
    res = c(res,tau_function(res$zeta,res$beta,res$psimat,res$mu,info,params))

    a_new = update_transitions(res$zeta,res$mu,info)

    # print(a_new)
    psi_new = update_emissions(res$tau,res$k_obs,res$n_obs,info, iteration)
    alpha_s_new = psi_new$alpha_new
    beta_s_new = psi_new$beta_new

    ### currently, we do not learn the bad states:
    badstates = c("badP","badM")
    parentstates = c("pat","mat")
    alpha_s_new[badstates] = alpha_s_new[parentstates]
    beta_s_new[badstates] = beta_s_new[parentstates]

    # initial Probs

    piv_new = update_piv(res$beta,res$psimat,res$mu,info,params)$piv
    pim_new = update_pim(res$mu)$pim
    params_new = list( a = a_new, piv = piv_new, pim = pim_new ,
                       alpha_s = alpha_s_new, beta_s = beta_s_new )
    iteration = iteration+1
    if (trace){
      cat(".")
      a_trace[iteration,] = a_new
      alpha_s_trace[iteration,] = alpha_s_new
      beta_s_trace[iteration,] = beta_s_new
      piv_trace[iteration,] = piv_new
      pim_trace[iteration,] = pim_new
    }
    differenz = max(abs(params$a - a_new)) +
      max(abs(params$piv -piv_new)) +
      max(abs(params$pim -pim_new)) +
      max(abs(params$alpha -alpha_s_new)) +
      max(abs(params$beta-beta_s_new))
    params = params_new

    if(verbose) cat("Absolute difference of the parameters after maximization: ", differenz, "\n")

    marker_annotation = vector("list",length(info$part_names))
    names(marker_annotation) = info$part_names
    for (part in info$part_names){
      hilf = apply(res$mu[[part]],1,which.max)
      marker_annotation[[part]] = colnames(res$mu[[part]])[hilf]
    }
    if(verbose){
      cat("Marker state distribution: \n")

      print(table(unlist(marker_annotation))/length(unlist(marker_annotation)))
      cat("\n")
    }
    if (differenz < eps) break()
    if(iteration > max.iter) {
      cat("\nReached maximum number of iterations (",max.iter,") without convergence.\n",sep="")
      break()
    }
  } # repeat

  # maximum posterior decoding for marker states M_t

  # hiddens = hidden_function(res$alpha,res$beta,info)

  result = list(params=params_new,
                marker_posterior = res$mu,
                marker_annotation = marker_annotation,
                # hidden_posterior = hiddens$hidden_posterior,
                # hidden_annotation = hiddens$hidden_annotation,
                alpha = res$alpha,
                beta = res$beta,
                psi = res$psimat,
                nalpha = res$nalpha,
                nbeta = res$nbeta)

  if (trace) result = c(result,list(iterations = iteration,
                                    a_trace=a_trace[1:iteration,],
                                    alpha_s_trace=alpha_s_trace[1:iteration,],
                                    beta_s_trace=beta_s_trace[1:iteration,],
                                    piv_trace=piv_trace[1:iteration,],
                                    pim_trace=pim_trace[1:iteration,]))
  return(result)
} # end EMalgorithm


#'
#'This function fits our extended HMM to an object with the observations
#'
#' @param object an RTIGER object with the filtered data
#' @param initial_params initial params to start the EM algorithm
#' @param randomize Logical value, wehter to randomize the initial parameters for the EM algorithm
#' @param eps convergence cutoff for EM algorithm. Difference of the the parameters between two consecutive iterations.
#' @param max.iter Maximum number of iterations.
#' @param trace Logical value, wehter the values of the parameters along the EM algorithm should be saved
#' @param verbose logical value, wether the algorithm status should be printed or not.
#'
#'
#'

fitModel = function( object,
                     initial_params = NULL,
                     randomize = FALSE,
                    eps = 0.01,
                    max.iter = 50,
                    trace = FALSE,
                    verbose = FALSE) {
  observations = object@FilteredData$as.mat
  info = object@info
  if(is.null(initial_params)) initial_params = generate_params(randomize = randomize)

  res = EMalgorithm(observations = observations,
                    info = info,
                    initial_params = initial_params,
                    max.iter = max.iter,
                    eps = eps,
                    trace = trace,
                    verbose = verbose)

  anno = seqlengths(object@FilteredData$as.GR[[1]])

  goodObs = lapply(names(observations), function(samp){
    newsamp = lapply(names(observations[[samp]]), function(chr){
      gM = res$marker_annotation[[chr]] == "+"
      gobs = observations[[samp]][[chr]][, gM]
      return(gobs)
    })# chr lapply
    names(newsamp) = names(observations[[samp]])
    return(newsamp)
  }) #samp lapply

  names(goodObs) = names(observations)

  # cat("GoodObs computed \n")
  goodPsis = lapply(names(res$psi), function(chr){
    gM = res$marker_annotation[[chr]] == "+"
    newPsi = res$psi[[chr]][, gM, ]
    return(newPsi)
  })
  names(goodPsis) = names(res$psi)
  # cat("GoodPsis computed\n")

  object@FilteredData$as.mat = goodObs
  res$psi = goodPsis
  cat("Viterbi decoding...\n\n")

  vit = Viterbi_decode(obs = goodObs, params =  res$params, psis = goodPsis)

  Viterbis = t(do.call(cbind, vit))
  tiles = width(object@FilteredData$as.GR[[1]])[1]
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


  object@FilteredData$as.GR = DataViterbi_GR
  object@Model = res
  object@Viterbi = DataViterbi_GR

  return(object)
}
