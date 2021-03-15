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
                     # Grangesobs = "list",
                     matobs = "list",
                     params = "list",
                     # FilteringThreshold = "list",
                     info = "list",
                     # Model = "list",
                     Viterbi = "list",
                     Probabilities = "list",
                     num.iter = "numeric"
                   ))

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


