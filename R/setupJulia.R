#' Installs the needed packages in JULIA to run the EM algorithm for rHMM.
#'
#' @export setupJulia
setupJulia = function(){
  julia_install_package_if_needed("Optim")
  julia_install_package_if_needed("Distributions")
  julia_install_package_if_needed("LinearAlgebra")
  julia_install_package_if_needed("CSV")
  julia_install_package_if_needed("DelimitedFiles")
  julia_install_package_if_needed("DataFrames")
  # julia_install_package_if_needed("Plots")
}

#' Function needed before using RTIGER() function. It loads the scripts in Julia that fit the rHMM.
#'
#'
#' @export sourceJulia
#'

sourceJulia = function(){
  julia_source(paste(system.file("julia", package = "RTIGER"),"/AuxilaryFunctions.jl", sep = ""))
  julia_source(paste(system.file("julia", package = "RTIGER"), "/rHMM_methods.jl", sep =""))
  # julia_source(paste(system.file("julia", package = "RTIGER"), "/rHMM_methods_old.jl", sep =""))

}
