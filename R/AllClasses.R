################################################################
####### RTIGER class
####### Constructor and show methods
################################################################

## TODO!!! If min.samples > total number of samples stop(error!!)

#' This class is a generic container for RTIGER analysis
#'
#' @slot RawData Raw data
#' @slot FilteredData List containing the filtered data based on the minimum number of samples supporting the each position and summing over the bin size.
#' @slot FilteringThreshold List of the parameters used to decide the bin size.
#' @slot info List with information about the samples.
#' @slot Model fitted model to the data.
#' @slot Viterbi List of chromosomes with the viterbi path per sample.
#' @rdname RTIGERDataSet
#' @exportClass RTIGER

.RTIGER = setClass("RTIGER",
         representation = representation(
           RawData = "list",
           FilteredData = "list",
           FilteringThreshold = "list",
           info = "list",
           Model = "list",
           Viterbi = "list"
         ))

#'
#'This function creates an RTIGER Object 
#'
#' @param experimentDesign Data frame with minimum two columns: files and name
#' @param observations A list of samples, which each sample is a list of matrix. Each matrix is a 2 by N markers. The first row is the parent 1 allele count and the second row is the total allele count.
#' @param GenRanges A GenomicRanges object with the positions of each marker.
#' @param bin.length Length for which start binning the genome and by which will be increased on each iteration
#' @param min.samples Minumum number of samples supporting an observation to be kept. Default is 1 sample.
#' @param min.counts Total allele count to select a criterion to keep increasing the bin size.
#' @param quant Quantile of observations that support the minimum counts. This quantile is used to decide when to stop increasing the bin size
#' @param seqlengths The genomic lengths of each chromosome (As in GenomicRanges).
#' 
#' 
#' #' @return A RTiger object with Genomic Ranges for the filtered and Raw data
#' @usage DataSetImportFromtxt(experimentDesign = NULL, observations = NULL, GenRanges = NULL, bin.length = NULL, min.samples = 1, min.counts = 10, quant = .2, seqlengths = NULL)
#' 
#'
#' @examples 
#' 
#' data("ATseqlengths")
#' path = system.file("extdata",  package = "RTIGER")
#' files = list.files(path, full.names = T)[1:3]
#' expDesign = data.frame(files = files, name = list.files(path)[1:3])
#' myDat = DataSetImportFromtxt(experimentDesign = expDesign, bin.length = 100, seqlengths = ATseqlengths)
#' 
#' 
#' @export DataSetImportFromtxt
#'

DataSetImportFromtxt = function(
  experimentDesign = NULL,
  observations = NULL,
  GenRanges = NULL,
  bin.length = NULL,
  min.samples = 1,
  min.counts = 10,
  quant = .2,
  seqlengths = NULL
){

  if( is.null(experimentDesign) & is.null(observations)){
    stop("No file information found!")
  }

# Check path and patterns -------------------------------------------------


  # if(!is.null(path)){
  #   if(!dir.exists(path)){
  #     stop("path to the files does not exist")
  #   }
  # 
  #   if(!length(list.files(path))){
  #     stop("No files in path directory")
  #   }
  # 
  #   if(!is.null(pattern)){
  #     if(!list.files(path = path, pattern = pattern)){
  #       stop("Do not exist files with this pattern")
  #     }
  #   }
  # 
  #   if(is.null(pattern)) pattern = ".txt"
  #   files = list.files(path, pattern, full.names = TRUE)
  #   myCol = data.frame(files = files, name = sapply(files, function(filepath) sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(filepath))))
  # }

# Check expDesign ---------------------------------------------------------

  if(!is.null(experimentDesign)){
    # cat("I am in experimentDesign\n")
    myCol = experimentDesign[, c("files", "name")]



# Load files, check rows and columns and create GR object -----------------
  # cat(myCol)



    rawGR = sapply(1:nrow(myCol), function(i){
      samp = as.character(myCol$files[i])
      print(paste("Loading file: ", myCol$files[i]))

      f <- read.delim(file =samp, header = FALSE)

      f = checkfileColumns(f)

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

  if(!is.null(bin.length)){
    listGR = lapply(rawGR, function(myG){
      binningFun(myG, bin.length = bin.length)
    })

    total = sapply(listGR, function(x) as.numeric(x$total))
    print(class(total))
    # if(class(total) == "list") total = unlist(total)

    old.bin.length = bin.length
    while(criterion(total, quant, min.counts)){
      bin.length = bin.length + old.bin.length


      cat(paste("New bin length = ", bin.length, "\n"))
      listGR = lapply(rawGR, function(myG){
        binningFun(myG, bin.length = bin.length)
      })

      total = sapply(listGR, function(x) as.numeric(x$total))
      # if(class(total) == "list") total = unlist(total)
    } # while loop

  } # !is.null(bin.length)

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

#' Create info list for observations
#'
#' @keywords internal
#' @noRd
#'

create_info <- function(obs){
  info = list()

  info$sample_names = names(obs)
  info$part_names = names(obs[[1]])
  info$marker_names = lapply(obs[[1]], colnames)
  info$sample_nr = length(obs)
  info$part_nr = length(obs[[1]])
  info$part_lengths = sapply(obs[[1]], ncol)
  info$NAs_matrices = lapply(obs, function(x) lapply(x, function(y){
    res = apply(y, 2, function(r) ifelse(is.na(r), NA, 1))
  }))
  return(info)
}

#' Prints description of RTIGER object
#'
#' @keywords internal
#' @noRd
#' 

setMethod(f = "show", signature = c("RTIGER"), function(object) {

  scat <- function(fmt, vals=character(), exdent=2, ...)
  {
    vals <- ifelse(nzchar(vals), vals, "''")
    lbls <- paste(S4Vectors:::selectSome(vals), collapse=" ")
    txt <- sprintf(fmt, length(vals), lbls)
    cat(strwrap(txt, exdent=exdent, ...), sep="\n")
  }

  cat("An object of class \"", class(object), "\" \n", sep = "")
  cat(" Number of samples:  ", object@info$sample_nr, "\n", sep = "")
  # cat(" Sample names: ", paste(as.character(object@info$sample_names), collapse = " "), "\n", sep = "")
  scat(" Sample names(%d): %s \n",object@info$sample_names)
  cat(" Number of chromosomes per sample: ", object@info$part_nr, "\n",
      sep = "")
  cat(" Chromosome names: ", paste(as.character(object@info$part_names), collapse = " "), "\n", sep = "")
  cat(" Number of observations per chromosome: ", paste(as.character(object@info$part_lengths), collapse = " "), "\n", sep = "")

})

