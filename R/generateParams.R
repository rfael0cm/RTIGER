###########################################
### Auxiliar functions to generate params
###########################################


#' Generate random parameters
#' @param params list of parameters for the HMM.
#' @param randomize logical variable wether the initial values should be taken randomly.
#'
#' @keywords internal
#' @noRd
#'

generate_params = function(params=list(),randomize = FALSE){
  vstates = c("pat","het","mat")
  mstates = c("+","badP","badM")
  estates = c("pat","het","mat","badP","badM")

  # transition matrix
  if (is.null(params$a)){
    # transition matrix a
    vstates = c("pat","het","mat")
    a = matrix(0.1,ncol=3,nrow=3) + diag(3)*10
    # if (randomize) a = a + matrix(runif(9,min = 0, max = 0.5),ncol=3)
    if (randomize) a = a + matrix(runif(9,min = 0, max = 0.01),ncol=3)
    a = a / rowSums(a)
    dimnames(a) = list(vstates,vstates)
    a["pat", "mat"] = a["mat", "pat"] = 1e-16
    a = a/rowSums(a)
    a["het", "het"] = mean(a["pat", "pat"], a["mat", "mat"])
    a = a/rowSums(a)
    params$a = a
  }

  # beta-binomial emission distributions
  if (is.null(params$alpha)){
    psi = vector("list", length(estates))
    priorstrength = 1 + randomize*runif(1,min = -0.5, max = 0.5)
    alphavec = c(20,20,1,20,1) + randomize*runif(5,min = -.5, max = .5)
    alpha_s = alphavec * priorstrength
    names(alpha_s) = estates
    betavec = c(1,20,20,1,20) + randomize*runif(5,min = -.5, max = .5)
    beta_s = betavec * priorstrength
    names(beta_s) = estates

    params$alpha_s = alpha_s
    params$beta_s = beta_s
  }

  # starting distribution for hidden states
  if (is.null(params$piv)){
    piv = c(1,1,1) + randomize * runif(3,min=-1,max=3)
    piv = piv / sum(piv)
    names(piv) = vstates
    params$piv = piv
  }

  # probabilities for the marker states
  if (is.null(params$pim)){
    pim = c(8,1,1) + randomize * runif(3,min=-1,max=3)
    pim = pim / sum(pim)
    names(pim) = mstates
    params$pim = pim
  }

  return(params)
} # end generate_params
