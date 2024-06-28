#' List of the 2D example for MaxWiK methods}
#'
#' A list containing input and output data for 2D example for Approximate Bayesian Computation,
#' including sampling scheme, meta-sampling, and prediction. To understand all details of the dataset, please, 
#' be kind to see vignette of the package. \cr 
#'
#' @format A list of:
#' \describe{
#'   \item{X}{ Input data frame of the model }
#'   \item{Y}{ Ouput data frame of the model }
#'   \item{observation}{ Data frame with observation info }
#'   \item{ABC}{ List of hyperparameters, the matrix of Voronoi sites, posteriori distribution, and results of MaxWiK algorithm }
#'   \item{metasampling}{ List of results of meta-sampling algorithm, and the network of points during meta-sampling }
#'   \item{sampling}{ List of object which are necessary for sampling algorithm like 
#'   function for simulation, parameters of the model, MSE (mean squared error), and X12 - generated points }
#'   \item{predictor}{ List of object which are necessary for predictor algorithm like
#'   posteriori.MaxWiK, result of the algorithm, and network of points during meta-sampling }
#' }
#'
'Data.2D'
