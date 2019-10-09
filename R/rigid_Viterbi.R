
### Efficient rigid Viterbi decoding

# Input: The emission matrix psi_k(o_t)
# The transition matrix A=(a_jk), possibly the initiation vector pi = (pi_k)
# rigidity parameter "rigid": The minimum length of a stretch of the same state
# before a transition may occur

# Output: The r-Viterbi path and its log-likelihood

rigid_Viterbi = function (psimat, transmat, pivec = NULL, rigid = 1){


  # log transformation of the parameters
  states = rownames(psimat) # check if there are rownames
  psimat = log(psimat + 1e-16)
  transmat = log(transmat + 1e-16)
  nstates = nrow(psimat)
  npos = ncol(psimat)
  if (is.null(pivec)) {pivec = rep(0,nstates)} else {pivec = log(pivec + 1e-16)}

  phi = matrix(0,nrow=nstates,ncol=npos)  # the matrix with the max probabilities
  rownames(phi) = states
  back = matrix(0,nrow=nstates,ncol=npos) # the backtracking pointers
  rownames(back) = states

  # hilf will contain the probabilities of the contiguous stretches of length rigid
  hilf = psimat
  hilf[,1] = hilf[,1] + pivec
  hilf[,2:npos] = hilf[,2:npos] + diag(transmat)
  hilf = t(apply(hilf,1,cumsum))
  hilf[,(rigid+1):npos] = hilf[,(rigid+1):npos] - hilf[,1:(npos-rigid)] # nochmal genau pr?fen
  # cat("hilf:\n")
  #print(hilf)

  # Initialize phi
  phi[,1] = hilf[,1]


  # phi for positions 2 : rigid still missing here
  sapply(2:rigid,function(k){

    # recursive calculation of the mat path
      stay = phi["mat",k-1] + transmat["mat","mat"]
      hetmat = hilf["het",k-1] + transmat["het","mat"]
      if (stay>hetmat) {
        phi["mat",k] <<- stay + psimat["mat",k]
        back["mat",k] <<- "mat"
      } else {
        phi["mat",k] <<- hetmat + psimat["mat",k]
        back["mat",k] <<- "het"
      }

    # recursive calculation of the pat path
      stay = phi["pat",k-1] + transmat["pat","pat"]
      hetpat = hilf["het",k-1] + transmat["het","pat"]
      if (stay>hetpat) {
        phi["pat",k] <<- stay + psimat["pat",k]
        back["pat",k] <<- "pat"
      } else {
        phi["pat",k] <<- hetpat + psimat["pat",k]
        back["pat",k] <<- "het"
      }

    # recursive calculation of the het path
      stay = phi["het",k-1] + transmat["het","het"]
      mathet =  hilf["mat",k-1] + transmat["mat","het"]
      pathet =  hilf["pat",k-1] + transmat["pat","het"]
      if (stay>= max(mathet,pathet)) {
        phi["het",k] <<- stay + psimat["het",k]
        back["het",k] <<- "het"
      } else if (mathet>pathet) {
        phi["het",k] <<- mathet + psimat["het",k]
        back["het",k] <<- "mat"
      } else {
        phi["het",k] <<- pathet + psimat["het",k]
        back["het",k] <<- "pat"
      }
    #

  } ) #end sapply


  # phi for positions (rigid+1) until npos
  sapply((rigid+1):npos,function(k){
      # cat("pos = ", k, "\n\n")
    # recursive calculation of the mat path
      stay = phi["mat",k-1] + transmat["mat","mat"]
      hetmat = phi["het",k-rigid] + hilf["het",k-1] + transmat["het","mat"]
      if (stay>hetmat) {
        phi["mat",k] <<- stay + psimat["mat",k]
        back["mat",k] <<- "mat"
      } else {
        phi["mat",k] <<- hetmat + psimat["mat",k]
        back["mat",k] <<- "het"
      }

    # recursive calculation of the pat path
      stay = phi["pat",k-1] + transmat["pat","pat"]
      hetpat = phi["het",k-rigid] + hilf["het",k-1] + transmat["het","pat"]
      if (stay>hetpat) {
        phi["pat",k] <<- stay + psimat["pat",k]
        back["pat",k] <<- "pat"
      } else {
        phi["pat",k] <<- hetpat + psimat["pat",k]
        back["pat",k] <<- "het"
      }

    # recursive calculation of the het path
      stay = phi["het",k-1] + transmat["het","het"]
      mathet = phi["mat",k-rigid] + hilf["mat",k-1] + transmat["mat","het"]
      pathet = phi["pat",k-rigid] + hilf["pat",k-1] + transmat["pat","het"]
      if (stay>= max(mathet,pathet)) {
        phi["het",k] <<- stay + psimat["het",k]
        back["het",k] <<- "het"
      } else if (mathet>pathet) {
        phi["het",k] <<- mathet + psimat["het",k]
        back["het",k] <<- "mat"
      } else {
        phi["het",k] <<- pathet + psimat["het",k]
        back["het",k] <<- "pat"
      }
  } ) #end sapply

  #print(phi)
  #print(back)

  #backtracking
  viterbipath = character(npos)
  currentpos = npos
  currentstate = states[which.max(phi[,currentpos])]
  viterbipath[currentpos] = currentstate
  loglikelihood = phi[currentstate,currentpos]
  repeat{
    zur = back[currentstate,currentpos]
    if (zur == currentstate) {
      currentpos = currentpos - 1
      viterbipath[currentpos] = currentstate
    } else {
      jump = max(1,currentpos-rigid)
      viterbipath[jump:(currentpos-1)] = zur
      currentpos = jump
      currentstate = zur
    }
    if (currentpos == 1) break()
  }

  return(list(viterbipath=viterbipath,loglikelihood=loglikelihood))
  # return(viterbipath)

} # end rigid_Viterbi


