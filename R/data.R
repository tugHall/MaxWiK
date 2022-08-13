#' Data frame with results of function \code{experiment_models()}
#'
#' A dataset containing results of function \code{experiment_models()} with next input: \cr
#' \code{ Results_toy_experiments  =  experiment_models( file_name = '../Results_ALL.txt', } \cr
#' \code{                           models = c( 'Gaussian', 'Linear' ), } \cr
#' \code{                           dimensions = (1:10)*2, } \cr
#' \code{                           stochastic_terms  =  c( 0, 1, 5, 10, 20, 30 ), } \cr
#' \code{                           rng  =  c( 0,10 ), } \cr
#' \code{                           restrict_points_number = 300 ) ) } \cr
#'
#' @format A data frame with 12960 rows and 8 variables:
#' \describe{
#'   \item{method_name}{ Name of a method for calculation }
#'   \item{kernel_name}{ Name of using kerel or '' if not applicapable }
#'   \item{model_name}{ Name of a model, can be 'Gaussian' or 'Linear' }
#'   \item{dimension}{ Dimension of the toy task }
#'   \item{stochastic_term}{ Stochastic term in a function of a model }
#'   \item{MSE}{ Mean square error }
#'   \item{running_time}{ Running time in seconds }
#'   \item{iteration}{ Iteration number }
#' }
'Results_toy_experiments'
