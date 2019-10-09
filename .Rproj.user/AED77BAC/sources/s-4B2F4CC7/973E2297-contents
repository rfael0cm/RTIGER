
#'  
#' The function is used to fit (bidirectional) Hidden Markov Models, given one or more observation sequence. 
#' 
#' @title Fit a Hidden Markov Model
#' 
#' @param obs The observations. A list of one or more entries containing the observation matrix (\code{numeric}) for the samples (e.g. chromosomes).
#' @param hmm The initial Hidden Markov Model. This is a \code{\linkS4class{HMM}}.
#' @param convergence Convergence cutoff for EM-algorithm (default: 1e-6).
#' @param maxIters Maximum number of iterations.
#' @param dirFlags The flag sequence is needed when a bdHMM is fitted on undirected data (e.g.) ChIP only. It is a \code{list} of \code{character} vectors indication for each position its knwon directionality. U allows all states. F allows undirected states and states in forward direction. R allows undirected states and states in reverse direction. 
#' @param emissionProbs List of precalculated emission probabilities of emission function is of type 'null'.
#' @param effectiveZero Transitions below this cutoff are analytically set to 0 to speed up comptuations.
#' @param verbose \code{logical} for printing algorithm status or not.
#' @param nCores Number of cores to use for computations.
#' @param incrementalEM When TRUE, the incremental EM is used to fit the model, where parameters are updated after each iteration over a single observation sequence.
#' @param updateTransMat Wether transitions should be updated during model learning, default: TRUE.
#' @param sizeFactors Library size factors for Emissions PoissonLogNormal or NegativeBinomial as a length(obs) x ncol(obs[[1]]) matrix.
#' @param clustering Boolean variable to specify wether it should be fit as an HMM or or bdClustering. Please, use function bdClust when bdClust is prefered.
#' 
#' @return A list containing the trace of the log-likelihood during EM learning and the fitted HMM model.
#' @usage fitHMM(obs=list(), hmm, convergence=1e-6, maxIters=1000, dirFlags=list(), emissionProbs=list(), effectiveZero=0, verbose=FALSE, nCores=1, incrementalEM=FALSE, updateTransMat=TRUE, sizeFactors=matrix(1, nrow=length(obs), ncol=ncol(obs[[1]])), clustering = FALSE)
#' 
#' @seealso \code{\linkS4class{HMM}}
#' @examples 
#' 
#' data(example)
#' hmm_ex = initHMM(observations, nStates=3, method="Gaussian")
#' hmm_fitted = fitHMM(observations, hmm_ex)
#'
#' @export fitHMM
fitHMM = function(obs=list(), hmm, convergence=1e-6, maxIters=1000, dirFlags=list(), emissionProbs=list(), effectiveZero=0, verbose=FALSE, nCores=1, incrementalEM=FALSE, updateTransMat=TRUE, sizeFactors=matrix(1, nrow=length(obs), ncol=ncol(obs[[1]])), clustering = FALSE){
    
    myFile = file.path(tempdir(),paste("STAN.temp.params.", Sys.getpid(), ".rda",sep="") )
    myparams = EmissionParams(hmm)
    save(myparams, file=myFile)

    
    originalStateLabel = NULL
    if(class(hmm) == "bdHMM") {
        originalStateLabel = hmm@stateNames
        hmm = reorderbdHMMStates(hmm, originalStateLabel, FALSE)
    }
    if(class(obs) != "list") stop("Observations must be submitted as a list of matrices!")
    if(! all(sapply(obs, class) == "matrix")) stop("Observations must be submitted as a list of matrices!")
    for(i in 1:length(obs)) {
        for(j in 1:ncol(obs[[i]])) {
            if(class(obs[[i]][,j]) != "numeric") {
                warning("Data track ", j, "in observation sequence ", i, " is not numeric.")
                obs[[i]][,j] = as.numeric(obs[[i]][,j])
            }   
        }
    }
    
    if(!all(unlist(lapply(obs, function(x) apply(x,2,class))) == "numeric")) {
        stop("All observations must be of class numeric!")
    }
    
    if(! class(hmm) %in% c("HMM", "bdHMM")) stop("Parameter hmm is not of class HMM or bdHMM!")
    if(length(dirFlags) > 0) {
        if(length(dirFlags) != length(obs)) stop("Flag sequence and observations must have the same length!")
        if(sum(sapply(dirFlags, length)) != sum(sapply(obs, nrow))) stop("Flag sequence and observations must have the same length!")
    }
    if(length(emissionProbs) > 0) {
        if(length(emissionProbs) != length(obs)) stop("Pre-computed emission probabilities and observations must have the same length!")
    #   if(sum(sapply(emissionProbs, length)) != sum(sapply(obs, nrow))) stop("Pre-computed emission probabilities and observations must have the same length!")
        if(! all(unique(sapply(emissionProbs, ncol)) == hmm@nStates)) stop("Number of columns in emissionProbs object must equal number of HMM states.")
    }
    if(length(unique(sapply(obs, ncol))) != 1) stop("All entries of observation sequences must have the same number of columns!")
    if(ncol(obs[[1]]) != sum(sapply(hmm@emission, function(y) y@dim))) stop("Number of columns in observation sequence must match dimensions of observations sequence!")
    
    for(i in 1:length(obs)) {
        if(! all(apply(obs[[i]], 2, is.numeric))) stop("Data tracks are not numeric.")
    }
	for(i in 1:length(hmm@emission)) {
		if(!all(sizeFactors == 1) & (! hmm@emission[[i]]@type %in% c("NegativeBinomial", "PoissonLogNormal"))) stop("Size factor correction is only implemented for NegativeBinomial and PoissonLogNormal distributions!\n")
	}
	
    ## Every HMM object contains JointlyIndependent Emissions
    myParameters = list(emissions=hmm@emission)
#   if(any(sapply(hmm@emission, function(x) x@type) == "null")) {
        jiEmission = HMMEmission(type="JointlyIndependent", parameters=myParameters, nStates=as.integer(hmm@nStates))
        hmm@emission = list(jiEmission)
#   }
    
    emission = hmm@emission[[1]]
    emissionParams = emission@parameters
    emissionParams = prepareEmission(emissionParams, emission@type)
    
    if(emission@type == "null" & length(emissionProbs) == 0) stop("Must supply pre-calculated emission probabilities when using emission function of type 'null'!")
    if(emission@type == "null" & (length(emissionProbs) != length(obs))) stop("Observation sequence must have same length as pre-calculated emissions!")
        
    D = hmm@emission[[1]]@dim
    initProb = hmm@initProb
    transMat = hmm@transMat
    nStates = hmm@nStates
    
    
    bdHMM.settings=list(dirFlags=dirFlags, stateNames=character(), bidirectional.mc=TRUE, optim.method="", rev.operation=NULL, bidirOptimParams=list())
    
    if(class(hmm) == "bdHMM") {
        bdHMM.settings$optim.method = hmm@transitionsOptim
        
        if(length(bdHMM.settings$dirFlags) > 0) {
            if(! all(sapply(bdHMM.settings$dirFlags, class) %in% c("character", "factor"))) {
                stop("Flags must be submitted as list of characters or factors!\n")
            }
            bdHMM.settings$dirFlags = lapply(bdHMM.settings$dirFlags, factor)
        }
        if(length(bdHMM.settings$dirFlags) > 0 & (! all(unlist(bdHMM.settings$dirFlags) %in% c("F", "R", "U")))) {
            stop("dirFlags must be one of c(\"F\", \"R\", \"U\")!\n")
        }
        else {
            flag2int = c(1,-1,0)
            names(flag2int) = c("F","R","U")
            bdHMM.settings$dirFlags = lapply(bdHMM.settings$dirFlags, function(x) as.integer(flag2int[x]))
            #print(dirFlags)
        }
        
        bdHMM.settings$stateNames = hmm@stateNames
        bdHMM.settings$bidirectional.mc = hmm@transitionsOptim %in% c("rsolnp", "analytical")
        bdHMM.settings$rev.operation = 1:sum(emission@dim) ##changed Julia
        for(i in 1:length(hmm@directedObs)) {
            if(hmm@directedObs[i] == 0) {
                bdHMM.settings$rev.operation[i] = i
            }
            else {
                myPair = which(hmm@directedObs == hmm@directedObs[i])
                myPair = myPair[myPair != i]
                bdHMM.settings$rev.operation[i] = myPair
            }
        }
    #   print(bdHMM.settings$rev.operation)
    }

    observationEmissionType = unlist(sapply(emission@parameters$emissions, function(x) rep(x@type, x@dim)))
    if(any(length(observationEmissionType) != sapply(obs, ncol))) stop("Number of data tracks (columns in observation matrix) do not match dimension of (bd)HMM.")
    
    couples = NULL
#   print("1")
    if(class(hmm) == "bdHMM") {
        bdHMM.settings$state2flag=rep(0,nStates)
        bdHMM.settings$couples=NULL
        bdHMM.settings = bdHMM_get_info(bdHMM.settings)
        if(any(bdHMM.settings$rev.operation < 1)) stop("Reverse oparations contain dimensions < 0!\n")
        if(! is.null(bdHMM.settings$rev.operation)) {
            emissionParams = emissionRevOp(emissionParams, emission@type, bdHMM.settings, nStates, hmm@directedObs, NULL, hmm@stateNames)
            bdHMM.settings$rev.operation = bdHMM.settings$rev.operation - 1
        }
        if(bdHMM.settings$optim.method == "analytical") {
            bdHMM.settings$bidirOptimParams=list()
        }
        else {
            bdHMM.settings$bidirOptimParams$c2optimize = c2optimize
        }
    }

#   print("2")
    #emissionPrior = emissionParams$emissionPrior
#   if(is.null(emissionPrior)) {
#       emissionPrior = list(useLongPrior=as.integer(FALSE))
#   }
#   else {
#       emissionPrior$calldiwish = calldiwish
#       emissionPrior$useLongPrior = FALSE
#   }
    
    mySplit = list()
    if(emission@type == "JointlyIndependent") {
        # print("I am in JointlyInd")
        if(length(observationEmissionType) == 0) stop("Please specify observationEmissionType to ensure that JointlyIndependent emission match observation data types!")
        if(length(observationEmissionType) != dim(obs[[1]])[2]) stop("Length of observationDataType does not match number of observation cols (data tracks)!")
        emissionTypes = unlist(sapply(emission@parameters$emissions, function(x) rep(x@type, x@dim)))
        matchEmission2Data(obs, emissionTypes, verbose)
        
        if(! all(emissionTypes == observationEmissionType)) stop("Types of emission functions does not match emission types of observations!")
        ## TODO: check what happens when dimensionality does not match
        for(currType in c("NegativeBinomial", "PoissonLogNormal")) {
            
            myD = which(observationEmissionType == currType)
            if(length(myD) > 0) {
                for(i in 1:length(emissionParams$emissions)) {
                    emissionParams$emissions[[i]]@parameters$sizeFactor = lapply(1:nStates, function(x) sizeFactors[,i])
                    #print(emissionParams$emissions[[i]]@parameters)
                }
            }
            
            if(length(myD) > 0) {
                if(verbose) cat("Preparing unique count map for ", paste(unique(observationEmissionType[myD]), collapse=" & "), "\n", sep="")
                tempObs = lapply(obs, function(x) x[,myD])
                if(length(myD) == 1) {
                    tempObs = lapply(tempObs, function(x) matrix(x, ncol=1))
                }
            #   print(bdHMM.settings$rev.operation)
                if(class(hmm) == "bdHMM") {
                    # print(bdHMM.settings$rev.operation)
                    if(any( (bdHMM.settings$rev.operation[myD]+1) != (1:length(bdHMM.settings$rev.operation))[myD])) {
                        tempObs = c(tempObs, tempObs)
                        #print(bdHMM.settings$rev.operation)
                        tempObs[(length(obs)+1):(2*length(obs))] = lapply(tempObs[(length(obs)+1):(2*length(obs))], function(x) x[,bdHMM.settings$rev.operation[myD]+1])
                        #print("1")
                    }
                }
                
                    
                #print(myD)
                myCurrSplit = lapply(tempObs, function(myMat) lapply(1:ncol(myMat), function(d) as.list(tapply(1:nrow(myMat), INDEX=myMat[,d], identity))))
                if(length(myCurrSplit) > 1) {
                    for(i in 2:length(myCurrSplit)) {
                        addThis = max(unlist(myCurrSplit[[i-1]]))
                        for(d in 1:length(myD)) {
                            for(j in 1:length(myCurrSplit[[i]][[myD[d]]])) {
                                myCurrSplit[[i]][[d]][[j]] = myCurrSplit[[i]][[myD[d]]][[j]]+addThis
                            }
                        }
                    }
                }
                for(i in 1:length(myCurrSplit)) {
                    tempList = rep(list(), length(myCurrSplit[[i]]))
                    tempList[myD] = myCurrSplit[[i]]
                    myCurrSplit[[i]] = tempList
                }
                
            #   reorderObs = unlist(lapply(1:length(obs), function(x) c(x, x+length(obs))))
            #   myCurrSplit = myCurrSplit[reorderObs]
                # browser()
                # print(length(myCurrSplit))
                mySplit[[currType]] = list(countSplit=myCurrSplit, optimFct=ifelse(currType=="PoissonLogNormal",optimizePoiLog,optimizeNB))
            }
            
        }
        
        emissionParams$mySplit = mySplit
    }
    # trace("myQNBinom", quote(if(any(is.nan(myGammas))) {recover()}), at = 3, print = F)
    hmm_out = .Call("RHMMFit", SEXPobs=obs, SEXPpi=initProb, SEXPA=transMat, SEXPemission=emissionParams, SEXPtype=as.character(emission@type), SEXPdim=D, SEXPregularize=as.numeric(0), SEXPk=as.integer(nStates), SEXPmaxIters=as.integer(maxIters), SEXPparallel=as.integer(nCores), SEXPflags=lapply(bdHMM.settings$dirFlags, as.integer), SEXPstate2flag=as.integer(bdHMM.settings$state2flag), SEXPcouples=as.integer(bdHMM.settings$couples), SEXPrevop=as.integer(bdHMM.settings$rev.operation), SEXPverbose=as.integer(verbose), SEXPupdateTransMat=as.integer(updateTransMat), SEXPfixedEmission=emissionProbs, SEXPbidiroptim=bdHMM.settings$bidirOptimParams, SEXPemissionPrior=sizeFactors, SEXPeffectivezero=as.numeric(effectiveZero), SEXPconvergence=as.numeric(convergence), SEXPincrementalEM=as.integer(incrementalEM), SEXPclustering = as.integer(clustering), PACKAGE="STAN") 
    
    if(class(hmm) == "bdHMM") {
        hmm@dirScore = hmm_out$dirScore 
        #print(hmm@dirScore)
    }
    
    hmm@transMat = matrix(hmm_out$transMat, nrow=nStates, ncol=nStates, byrow=TRUE)
    hmm@initProb = hmm_out$initProb
    #print(hmm_out$emission)
    fixed.emission = list()
    if(length(emissionProbs) == 0) {
        hmm@emission = list(reformatEmissionAfterFitting(hmm_out$emission, emission@type, D, nStates))#hmm_out$emission#
    }
    # print(hmm@emission)
    if(! is.null(bdHMM.settings$rev.operation)) {
        bdHMM.settings$rev.operation = bdHMM.settings$rev.operation + 1
        hmm@emission[[1]]@parameters = emissionRevOp(hmm@emission[[1]]@parameters, emission@type, bdHMM.settings, nStates, hmm@directedObs, colnames(obs[[1]]), hmm@stateNames)
    }
    
    if(emission@type != "null") {
        hmm@emission = hmm@emission[[1]]@parameters$emissions
    }
    else {
        hmm@emission = myEmission = HMMEmission(type="null", parameters=list(dim=as.integer(D)),nStates=nStates)
    }
    hmm@status = paste("EM-fit completed: ", Sys.time(), sep="")
    #hmm@converged = (fit$loglik[length(fit$loglik)] - fit$loglik[length(fit$loglik)-1]) < 1e-6
    if(class(hmm) == "bdHMM") {
        hmm = reorderbdHMMStates(hmm, originalStateLabel, TRUE)
        isvalid = validbdHMM(hmm)
    #   print(hmm@dirScore)
    }
    else {
        isvalid = validHMM(hmm)
    }
    
    hmm@LogLik = hmm_out$loglik
    system(paste("rm ", myFile, sep=""))
    hmm
}


#' Internal function matches type of emission function to data type.
#' 
#' @keywords internal
#' @noRd
matchEmission2Data = function(obs, emissionTypes, verbose=FALSE) {
    allData = do.call("rbind", obs)
    dataTypes = getDataType(allData)
    emission2Type = list("null"="null", "Gaussian"="continuous", "Bernoulli"="binary", "NegativeBinomial"="discrete", "PoissonLogNormal"="discrete", "Poisson"="discrete", "Multinomial"=c("discrete", "binary"))
    if (verbose) cat("Matching data to emission distributions:\n")
    if(length(colnames(allData)) == 0) {
        colnames(allData) = 1:ncol(allData)
    }
    myL = max(nchar(colnames(allData)))
    myLemission = max(nchar(emissionTypes))
    
    for(i in 1:ncol(allData)) {
        if(! emissionTypes[i] %in% names(emission2Type)) stop("Unknown emission type specified. Not able to match emission type to data type.")
        myI = i
        if(! is.null(colnames(allData))) myI = paste(c(colnames(allData)[i], rep(" ", myL-nchar(colnames(allData)[i]))), collapse="")
        if(verbose) cat("dim=", myI, "\t", paste(c(emissionTypes[i], rep(" ", myLemission-nchar(emissionTypes[i])))), "\t", dataTypes[i], "\n", sep="")
    #   if(! dataTypes[i] %in% emission2Type[[emissionTypes[i]]]) warning("Using dubious emission model: ", emissionType[[emissionTypes[i]]], " models ", emission2Type[[emissionType]], " data in column ", i, ".", sep="")
        if(any(allData[,i] < 0, na.rm=TRUE) & emissionTypes[i] %in% c("NegativeBinomial", "PoissonLogNormal")) stop("Negative values in discrete data track for emission ", emissionTypes[i], ".", sep="")
    }
}


#' Internal function get data type of observations.
#' 
#' @keywords internal
#' @noRd
getDataType = function(mydata) {
    type = c()
    for(i in 1:ncol(mydata)) {
        if(all(mydata[!is.na(mydata[,i]),i] %in% c(0,1))) {
            type[i] = "binary"
        }
        else if(all(mydata[!is.na(mydata[,i]),i] == as.integer(mydata[!is.na(mydata[,i]),i]))) {
            type[i] = "discrete"
        }
        else if(is.numeric(mydata[!is.na(mydata[,i]),i])) {
            type[i] = "continuous"
        }
    }
    type
}


#' Prepares and reformats emission for model fitting. 
#' 
#' @keywords internal
#' @noRd
prepareEmission = function(emissionParams, type) {
    if(type == "Gaussian" & (all(c("mu", "cov") %in% names(emissionParams)))) {
        emissionParams$cov = lapply(emissionParams$cov, function(x) as.matrix(x, mode="numeric"))
        emissionParams$mu = lapply(emissionParams$mu, as.numeric)
        if(! "updateCov" %in% names(emissionParams)) {
            emissionParams$updateCov = as.integer(1)
        }
        else {
            emissionParams$updateCov = as.integer(emissionParams$updateCov)
        }
        if(! ("sharedCov" %in% names(emissionParams))) {
            emissionParams[["sharedCov"]] = FALSE
        }
        emissionParams$sharedCov = as.integer(emissionParams$sharedCov)
        
    }
    if(type == "JointlyIndependent") {
        for(i in 1:length(emissionParams$emissions)) {
            emissionParams$emissions[[i]]@parameters = prepareEmission(emissionParams$emissions[[i]]@parameters, emissionParams$emissions[[i]]@type)
        }
    }
    emissionParams
}

#' Processing and reformatting of emission after model fitting. 
#' 
#' @keywords internal
#' @noRd
reformatEmissionAfterFitting = function(emissionParams, type, D, nStates) {
    
    myEmission = list()
    if(type == "JointlyIndependent" & (length(emissionParams) != 0)) {
    types = emissionParams$types
    for(i in 1:length(emissionParams$emissions)) {
        if(types[i] == "Bernoulli") {
            emissionParams$emissions[[i]] = unlist(emissionParams$emissions[[i]], recursive=FALSE)
            emissionParams$emissions[[i]] = list(p=emissionParams$emissions[[i]])
        }
        if(types[i] == "Poisson") {
            emissionParams$emissions[[i]] = unlist(emissionParams$emissions[[i]], recursive=FALSE)
            emissionParams$emissions[[i]] = list(lambda=emissionParams$emissions[[i]])
        }
        if(types[i] == "Multinomial") {
            emissionParams$emissions[[i]] = unlist(emissionParams$emissions[[i]], recursive=FALSE)
            emissionParams$emissions[[i]] = list(p=emissionParams$emissions[[i]][seq(1, length(emissionParams$emissions[[i]])-1, by=2)],
                    reverseComplementary=emissionParams$emissions[[i]][seq(2, length(emissionParams$emissions[[i]]), by=2)][[1]])
            for(k in 1:length(emissionParams$emissions[[i]])) {
                emissionParams$emissions[[i]][[k]] = lapply(emissionParams$emissions[[i]][[k]], unlist)
            }
            # emissionParams$emissions[[i]] = lapply(emissionParams$emissions[[i]], unlist)              
        }
        if(types[i] == "NegativeBinomial") {
            emissionParams$emissions[[i]] = unlist(emissionParams$emissions[[i]], recursive=FALSE)
            emissionParams$emissions[[i]] = list(mu=emissionParams$emissions[[i]][seq(1, length(emissionParams$emissions[[i]])-3, by=4)], 
                                                 size=emissionParams$emissions[[i]][seq(2, length(emissionParams$emissions[[i]])-2, by=4)], 
                                                 pi=emissionParams$emissions[[i]][seq(4, length(emissionParams$emissions[[i]]), by=4)])
            for(k in 1:length(emissionParams$emissions[[i]])) {
                emissionParams$emissions[[i]][[k]] = lapply(emissionParams$emissions[[i]][[k]], unlist)
            }
        }
        if(types[i] == "PoissonLogNormal") {
            emissionParams$emissions[[i]] = unlist(emissionParams$emissions[[i]], recursive=FALSE)
            emissionParams$emissions[[i]] = list(mu=emissionParams$emissions[[i]][seq(1, length(emissionParams$emissions[[i]])-2, by=3)], 
                    sigma=emissionParams$emissions[[i]][seq(2, length(emissionParams$emissions[[i]])-1, by=3)])#,
            for(k in 1:length(emissionParams$emissions[[i]])) {
                emissionParams$emissions[[i]][[k]] = lapply(emissionParams$emissions[[i]][[k]], unlist)
            }
        }
        if(types[i] == "Gaussian") {
            emissionParams$emissions[[i]] = list(mu=lapply(emissionParams$emissions[[i]], function(x) x[[1]][[1]]), 
                    cov=lapply(emissionParams$emissions[[i]], function(x) x[[2]][[1]]))
            emissionParams$emissions[[i]]$cov = lapply(emissionParams$emissions[[i]]$cov, function(x) matrix(x, nrow=sqrt(length(x)), ncol=sqrt(length(x)), byrow=TRUE))
            if(length(unique(emissionParams$emissions[[i]]$cov)) == 1) {
                emissionParams$emissions[[i]][["sharedCov"]] = TRUE
            }
            else {
                emissionParams$emissions[[i]][["sharedCov"]] = FALSE
            }
        }
        
    }
    myEmissionList = lapply(1:length(types), function(x) HMMEmission(type=types[x], parameters=emissionParams$emissions[[x]], nStates=nStates))
    myEmission = HMMEmission(type="JointlyIndependent", parameters=list(emissions=myEmissionList), nStates=as.integer(nStates))
        
    }
    if(type == "JointlyIndependent" & (length(emissionParams) == 0)) {
        myEmission = HMMEmission(type="null", parameters=list(dim=as.integer(D)),nStates=nStates)
    }
    myEmission
    
}


#' Reverse operation on observations before and after model fitting.
#' 
#' @keywords internal
#' @noRd
emissionRevOp = function(emissionParams, type, bdHMM.settings, nStates, directedObs, dimNames, stateNames) {
    if(type == "JointlyIndependent") {
        myDimSets = list()
        revop = bdHMM.settings$rev.operation
        Ds = sapply(emissionParams$emissions, function(x) x@dim)
        myDims = cumsum(sapply(emissionParams$emissions, function(x) x@dim))
        currOffSet = 0
        for(i in 1:length(Ds)) {
            myDimSets[[i]] = (currOffSet+1):(currOffSet+Ds[i])
            currOffSet = currOffSet + length(myDimSets[[i]])
        }
        #print(Ds)
        #print(myDims)
        #print(myDimSets)
        #print(revop)
        
        
        emissionParamsTemp = emissionParams
        for(i in 1:length(myDimSets)) {
            for(p in 1:length(emissionParams$emissions[[i]]@parameters)) {
                if(! names(emissionParams$emissions[[i]]@parameters)[p] %in% c("sharedCov", "updateCov", "reverseComplementary")) {
                    for(k in 1:nStates) {
                    #   print(emissionParamsTemp$emissions[[i]]@parameters)
                        if(class(emissionParamsTemp$emissions[[i]]@parameters[[p]][[k]]) == "matrix") {
                            colnames(emissionParamsTemp$emissions[[i]]@parameters[[p]][[k]]) = rownames(emissionParamsTemp$emissions[[i]]@parameters[[p]][[k]]) = dimNames[myDimSets[[i]]]
                        }
                        else {
                            names(emissionParamsTemp$emissions[[i]]@parameters[[p]][[k]]) = dimNames[myDimSets[[i]]]
                        }
                        
                    }
                }
                
                
            }
        }
        for(i in 1:length(myDimSets)) {
            myDimSets[[i]] = revop[myDimSets[[i]]]
        }
    #   print(myDimSets)
        currMax = 0
        for(i in 1:length(myDimSets)) {
            #print(myDimSets)
            #print(currMax)
            if(length(myDimSets[[i]]) == 1) {
                if(myDimSets[[i]] != myDims[i]) {
                    for(p in 1:length(emissionParams$emissions[[i]]@parameters)) {
                        if(! names(emissionParams$emissions[[i]]@parameters)[p] %in% c("sharedCov", "updateCov", "reverseComplementary")) {
                            emissionParamsTemp$emissions[[i]]@parameters[[p]][grep("R", stateNames)] = emissionParams$emissions[[myDimSets[[i]]-currMax]]@parameters[[p]][grep("R", stateNames)]
                            
                        }
                    }
                }
            }
            else if(length(myDimSets[[i]]) > 1) {
                myObsPairs = directedObs[myDimSets[[i]]]
                currRevOp = bdHMM.settings$rev.operation[myDimSets[[i]]]
                currRevOp = currRevOp - currRevOp[1] + 1
                currRevOp = currRevOp[order(myDimSets[[i]])]
                if(any(table(myObsPairs[myObsPairs != 0]) != 2)) stop("Inconsistency of directedObs. No pairing possible within multi-dimensional emission ", emissionParams$emissions[[i]]@type, ".", sep="")
            #   print(currRevOp)
                emissionParamsTemp$emissions[[i]] = revOpMultiDim(emissionParamsTemp$emissions[[i]], currRevOp, nStates, bdHMM.settings$state2flag)
            }
            currMax = currMax+length(myDimSets[[i]])-1
        }
        emissionParams = emissionParamsTemp     
    }
        
    emissionParams
}


#' Internal function to preform reverse operation operator on multid-dimensional emissions.
#' 
#' @keywords internal
#' @noRd
revOpMultiDim = function(myEmission, revop, nStates, state2flag) {
    
    nStates = myEmission@nStates
    if(myEmission@type == "Gaussian" & (all(c("mu", "cov") %in% names(myEmission@parameters)))) {
        for(k in 1:nStates) {
            if(state2flag[k] == 1) {
                currNames = names(myEmission@parameters$mu[[k]])
                myEmission@parameters$mu[[k]] = myEmission@parameters$mu[[k]][revop]
                myEmission@parameters$cov[[k]] = myEmission@parameters$cov[[k]][revop,revop]
                names(myEmission@parameters$mu[[k]]) = currNames
                rownames(myEmission@parameters$cov[[k]]) = colnames(myEmission@parameters$cov[[k]]) = currNames
            }
        }
    }
    if(myEmission@type == "Multinomial" & (all(c("p") %in% names(myEmission@parameters)))) {
        #print("yo")
        for(k in 1:nStates) {
            if(state2flag[k] == 1) {
                currNames = names(myEmission@parameters$p[[k]])
        #       print(myEmission@parameters$p[[k]])
        #       print(myEmission@parameters$reverseComplementary)
                myEmission@parameters$p[[k]] = myEmission@parameters$p[[k]][revop][myEmission@parameters$reverseComplementary]
        #       print(myEmission@parameters$p[[k]])
                names(myEmission@parameters$p[[k]]) = currNames
            }
        }
    }
    myEmission
}



#' Internal function to define state flags for C++ code
#' 
#' @keywords internal
#' @noRd
getStateFlags = function(stateNames) {
    state2flag = rep(100, length(stateNames))
    state2flag[grep("F", stateNames)] = 1
    state2flag[grep("R", stateNames)] = -1
    state2flag[grep("U", stateNames)] = 100
    state2flag = -state2flag
    state2flag
}


#' Prepares settings for bdHMM transition optimization. 
#' 
#' @keywords internal
#' @noRd
bdHMM_get_info = function(bdHMM.settings) {
    
    stateNames = bdHMM.settings$stateNames
    Forward = sapply(strsplit(stateNames[grep("F", stateNames)], "F"), function(x) x[2])
    Reverse = sapply(strsplit(stateNames[grep("R", stateNames)], "R"), function(x) x[2])
    U = sapply(strsplit(stateNames[grep("U", stateNames)], "U"), function(x) x[2])
    nStates = length(stateNames)
    dirStates = (nStates - length(U))/2
    undirStates = length(U)
    myU = ""
    if(length(U) > 0) {
        myU = paste("U", 1:undirStates, sep="")
    }
    if(length(Forward) != length(Reverse)) {stop("Unequal numbers of states with forward and reverse orientation!")}
    
    bdHMM.settings$couples = (1:nStates)-1
    if(dirStates > 0) {
        bdHMM.settings$couples[1:(2*dirStates)] = c((1:dirStates)+dirStates, 1:dirStates)-1
    }
    
    
    bdHMM.settings$state2flag = getStateFlags(stateNames)
    

    if(bdHMM.settings$bidirectional.mc) {
            
        eqB = c(rep(0,nStates),0)
        statD = rep(1/nStates, nStates)
        xsi_sum = matrix(as.numeric(NA), nrow=nStates, ncol=nStates)
        initGamma = as.numeric(rep(0, nStates))
        constraints = getConstraints(stateNames)
        x0 = rep(1/nStates, length(constraints)+((nStates-length(U))/2)+length(U))
        LB=rep(0, length(x0))
        UB=rep(1, length(x0))
        control=list(trace=0, tol=1e-6) ## => Likelihood may decrease if estimate is not accurate enough ( solution: choose relative convergence criterion for optimizer )
        transMat = matrix(1/nStates, nrow=nStates, ncol=nStates)
        bdHMM.settings$bidirOptimParams = list(x0=x0, eqB=eqB, xsi_sum=xsi_sum, initGamma=initGamma, statD=statD, constraints=constraints, nStates=nStates, LB=LB, UB=UB, control=control, transMat=transMat, stateNames=stateNames, method=bdHMM.settings$optim.method, objective=as.numeric(Inf), nrm=as.integer(0), couples=bdHMM.settings$couples, update=as.integer(0))
        
    }
    bdHMM.settings
}



#' Internal function which reorders states before and after optimization
#' 
#' @keywords internal
#' @noRd
reorderbdHMMStates = function(bdhmm, originalLabels, reorderBack=FALSE) {
    if(class(bdhmm) != "bdHMM") stop("Object bdhmm must be of class bdHMM.")
    stateNames = bdhmm@stateNames
    stateNames2 = originalLabels
    reordStates = order(stateNames2)
    reordStatesBack = as.vector(tapply(1:length(reordStates), INDEX=reordStates, identity))
    reorderTo = reordStates
    if(reorderBack) {
        reorderTo = reordStatesBack
    }
#   print(reorderTo)
    if(!all(stateNames2[reordStates][reordStatesBack] == stateNames2)) stop("Inconsistency in labeling state directionality.")
    
    bdhmm[reorderTo,]
}
