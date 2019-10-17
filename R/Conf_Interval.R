conf.intervals = function(rtigerobj, rviterbiobj, conf = 0.95, max.width = 3, samples = NULL){
  if(class(rtigerobj) != "RTIGER") stop("rtiger is not a RTIGER object")
  if(class(rviterbiobj) != "RViterbi") stop("rviterbiobj is not a RViterbi")
  
  info = rtigerobj@info
  params = rtigerobj@Model$params
  DataViterbi_GR = rviterbiobj@Viterbi
  goodObs = rtigerobj@FilteredData$as.mat
  if(is.null(samples)) samples = info$sample_names 
  
  myconfs = lapply(samples, function(samp){
    cat(samp, "\n\n")
  # for(samp in info$sample_names){
    # final_conf = list()
    # for(chr in info$part_names){
    myres = lapply(info$part_names, function(chr){
      cat(chr, "\n\n")
      FinalRes = DataViterbi_GR[[samp]]
      FinalRes = FinalRes[seqnames(FinalRes) == chr]
      ln = table(FinalRes$Viterbi)
      if(length(ln) > 1){
        myGR = Vit2GrangesGen(FinalRes, "Viterbi")
        mypos = end(myGR)[-length(myGR)]
        mypos = which(end(FinalRes) %in% mypos)
        myres = alpha_beta_gamma_psi(goodObs[[samp]][[chr]], params = params)
        
        gammas = myres$gammas
        vstates = rownames(gammas)
        
        transMat = params$a
        # candidates = prescreen(gammas = gammas, maxwidth = max.width, threshold = 0.5)
        PosCand = sapply(mypos, function(x){
          nmax.width = max.width
          if(max.width > x) nmax.width = x-1
          if(max.width > (ncol(gammas) - x)) nmax.width = ncol(gammas) - x
          cat(x, "\n")
          sj = which(vstates %in% FinalRes$Viterbi[x])
          sk = which(vstates %in% FinalRes$Viterbi[x+1])
          myscore = wigglyTransition(a = x,
                                       b = x +1,
                                       sj = sj,
                                       sk = sk,
                                       alphas = myres$alphas,
                                       betas = myres$betas,
                                       normfactors = myres$normfactors, # scaling factors for the alphas
                                       transitions = transMat, # a_jk transition probs for moving from j to k
                                       psis = myres$psis # Emission Probabilities
                                       )
          i = 1
          newx = x
          newy = x + 1
          score.mat = data.frame(iteration = i, left_end = newx, right_end = newy, score = myscore)
          while(myscore < conf & diff(c(newx, newy)) < nmax.width){
            i = i + 1
            cat(i, "\n")
            if(i%%2) newx = newx -1
            else{ newy = newy + 1}
            myscore = intervalTransition(a = newx,
                                         b = newy,
                                         sj = sj,
                                         sk = sk,
                                         alphas = myres$alphas,
                                         betas = myres$betas,
                                         normfactors = myres$normfactors, # scaling factors for the alphas
                                         transitions = transMat, # a_jk transition probs for moving from j to k
                                         psis = myres$psis # Emission Probabilities
            )
            # cat(newx,",", newy, "\n", myscore, "\n\n")
            score.mat = rbind(score.mat, c(i, newx, newy, myscore ))
          } # while
          if(max(score.mat$score) > conf){
            myval = score.mat[which.max(score.mat$score),]
            return(c(score = myval$score, left_end = myval$left_end, right_end = myval$right_end, left_state = sj, right_state = sk ))
          } 
          else{
            newx = x
            newy = x + 1
            
            while(myscore < conf & diff(c(newx, newy)) < nmax.width){
              i = i + 1
              cat(i, "\n")
              newy = newy + 1
              myscore = intervalTransition(a = newx,
                                           b = newy,
                                           sj = sj,
                                           sk = sk,
                                           alphas = myres$alphas,
                                           betas = myres$betas,
                                           normfactors = myres$normfactors, # scaling factors for the alphas
                                           transitions = transMat, # a_jk transition probs for moving from j to k
                                           psis = myres$psis # Emission Probabilities
              )
              
              score.mat = rbind(score.mat, c(i, newx, newy, myscore ))
            } # while
            if(max(score.mat$score) > conf){
              myval = score.mat[which.max(score.mat$score),]
              return(c(score = myval$score, left_end = myval$left_end, right_end = myval$right_end, left_state = sj, right_state = sk ))
            }
            else{
              newx = x-1
              newy = x + 1
              
              while(myscore < conf & diff(c(newx, newy)) < nmax.width){
                i = i + 1
                cat(i, "\n")
                newx = newx -1
                myscore = intervalTransition(a = newx,
                                             b = newy,
                                             sj = sj,
                                             sk = sk,
                                             alphas = myres$alphas,
                                             betas = myres$betas,
                                             normfactors = myres$normfactors, # scaling factors for the alphas
                                             transitions = transMat, # a_jk transition probs for moving from j to k
                                             psis = myres$psis # Emission Probabilities
                )
                
                score.mat = rbind(score.mat, c(i, newx, newy, myscore ))
              } # while
              
              myval = score.mat[which.max(score.mat$score),]
              return(c(score = myval$score, left_end = myval$left_end, right_end = myval$right_end, left_state = sj, right_state = sk ))
            
            }#else
          } # else
        })
        
        
        PosCand = t(PosCand)
        
        goodTrans = as.data.frame(PosCand)
        goodTrans$left_state = vstates[goodTrans$left_state]
        goodTrans$right_state = vstates[goodTrans$right_state]
        
        
        final_conf = data.frame(seqnames = rep(chr, nrow(goodTrans)), 
                                       start = sapply(as.numeric(goodTrans$left_end), function(x) round(mean(c(start(FinalRes)[x], end(FinalRes)[x])), 0)) , 
                                       end = sapply(goodTrans$right_end, function(x) round(mean(c(start(FinalRes)[x], end(FinalRes)[x])), 0)),
                                       strand = "*",
                                       from = goodTrans$left_state,
                                       to = goodTrans$right_state,
                                       confidence = round(as.numeric(goodTrans$score) * 100, 2))
        return(final_conf)
      } # if there are CO
      else{ 
        final_conf = NULL
        return(final_conf)
      }
      
    }) # For loop chr
    
    final_conf = do.call(rbind, myres)
    final_conf = GRanges(final_conf)
    return(final_conf)
    
    } )# For loop samp
  names(myconfs) = samples
  return(myconfs)
  
}