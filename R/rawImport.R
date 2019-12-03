RawDataSetImport = function(
  experimentDesign = NULL,
  observations = NULL,
  GenRanges = NULL,
  min.samples = 1,
  seqlengths = NULL
){
  
  if( is.null(experimentDesign) & is.null(observations)){
    stop("No file information found!")
  }
  
  # Check expDesign ---------------------------------------------------------
  
  if(!is.null(experimentDesign)){
    # cat("I am in experimentDesign\n")
    myCol = experimentDesign[, c("files", "name")]
    
    rawGR = sapply(1:nrow(myCol), function(i){
      samp = as.character(myCol$files[i])
      print(paste("Loading file: ", myCol$files[i]))
      
      f <- read.delim(file =samp, header = FALSE)
      
      f = checkfileColumns(f, samp)
      
      myG = GRanges(seqnames =  f$V1,
                    ranges = IRanges(start = f$V2, end = f$V2),
                    P1.Allele = f$V3,
                    P1.Allele.Count = f$V4,
                    P2.Allele = f$V5,
                    P2.Allele.Count = f$V6,
                    seqlengths = seqlengths
      )
      
      return(myG)
    }) # end sapply rawGR
    
    names(rawGR) = myCol$name
    
    
  }
  
  if(all(names(seqlengths) != seqlevels(rawGR[[1]]))) stop("Names of the seqlengths differ from the seqnames in the files. \n")
  
  if(!is.null(observations)){
    
    rawGR = sapply(observations, function(obs){
      obs = t(do.call(cbind, obs))
      
      myG = GRanges(GenRanges,
                    # P1.Allele = NA,
                    P1.Allele.Count = obs[,1],
                    # P2.Allele = f$V5,
                    P2.Allele.Count = obs[,2] - obs[,1]
                    # seqlengths = seqlengths
      )
      
      return(myG)
    }) # end sapply rawGR
    
  }
  
 # !is.null(bin.length)
  listGR = rawGR
  
  chrPos = unlist(sapply(listGR, function(myG){
    vals = paste(seqnames(myG[!is.na(myG$total) ]), start(myG[!is.na(myG$total) ]), sep = "_")
    return(vals)
  }) #sapply
  )
  goods = which(table(chrPos) >= min.samples) # Filter positions by the number of samples that have sequences that position
  
  goodNam = names(goods)
  
  newGoodGR = lapply(listGR, function(myG){
    m = as.data.frame(myG)
    m$chrPos = paste(m$seqnames, m$start, sep = "_")
    m = m[m$chrPos %in% goodNam, ]
    m = GRanges(m)
    seqlengths(m) = seqlengths(myG)
    m = sort(m)
    return(m)
  })
  
  obs = lapply(newGoodGR, function(samp){
    chrs = lapply(seqlevels(samp), function(chr){
      myG = samp[seqnames(samp) == chr]
      p = as.matrix(mcols(myG[, c("P1.Allele.Count", "total")]))
      colnames(p) = c("paternal", "total")
      p[p[,"total"] == 0,] = NA
      
      p = t(p)
      colnames(p) = myG$chrPos
      
      return(p)
    })
    names(chrs) = seqlevels(samp)
    
    return(chrs)
  })
  
  # total = sapply(newGoodGR, function(x) as.numeric(x$total))
  # quantile(total[total > 0], quant, na.rm = TRUE)
  
  info = create_info(obs)
  
  filteringThreshold = list(min.counts = min.counts, quantile = quant, min.samples = min.samples)
  
  myObj <- .RTIGER(RawData = rawGR, FilteredData = list(as.GR = newGoodGR, as.mat = obs), info = info, FilteringThreshold = filteringThreshold, Model =list(), Viterbi = list())
  return(myObj)
}

