#' This class is a generic container for RTIGER analysis
#'
#' @slot matobs Nested lists. the first level is a list of samples. For each sample there are 5 matrices that contains the allele counts for each position.
#' @slot params a list with the parameters after training.
#' \itemize{
#' #' \item nstates: number of states.
#' \item logtransition: log transformed transition matrix.
#' \item paraBetaAlpha: Alpha value for the beta distribution.
#' \item paraBetaBeta: Beta value for the beta distribution.
#' \item rigidity: Rigity selected by the user.
#' \item logpi log transfrom initial probabilities.
#' \item transition: transition probability matrix.
#' \item pi: Initial probabilities.
#' }
#'
#'
#' @slot info List with phenotipic data of the samples.
#' \itemize{
#' #' \item sample_names.
#' \item part_names names of the chromosomes. Obtained from the raw data.
#' \item marker_names: Specific names for the markers. Chr_pos.
#' \item sample_nr Number of samples on the data set.
#' \item part_nr number of chromosomes
#' \item part_lenghts length of each chromsome
#' \item NAs_matrices Matrix that specifies where an NA is in each sample.
#' \item expDesing data frame introduced by the user.
#' }
#'
#' @slot Viterbi List of chromosomes with the viterbi path per sample.
#' @slot Probabilities Computed probabilites for the EM algorithm.
#' @slot num.iter Number of iterations needed to stop the EM algorithm.
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


