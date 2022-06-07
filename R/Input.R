

# MODELS ------------------------------------------------------------------


#' The model of simulations that is based on Gaussian functions for each dimension 
#'
#' @description The function \code{Gaussian_model()} allows to generate 
#' parameters and statistics of simulations that are based on Gaussian function for each dimension: \cr
#' \code{f( x ) = { exp( - (x1-x01)^2 / 2 ), ..., exp( - (xn-x0n)^2 / 2 )  } = { y1, y2, ..., yn } } is a vector of output data
#' 
#'
#' @param d Dimension of the parameter and model space
#' @param x0 Numeric vector with length of dimensionality of data frame, 
#' that contents the truth value of parameter. Each number in the vector should be within the range r
#' @param probability Logical, if TRUE then apply uneven distribution for parameters generation
#' @param n Integer number of points in data frames
#' @param r Range \code{r = c(min, max)}, by default  \code{r = range(0,10)}
#'
#' @return The function \code{Gaussian_model()} returns list of two objects: \cr
#' - stat.sim - data frame of simulations statistics,
#' - par.sim - data frame of parameters.
#' 
#' @export
#'
#' @examples
#' NULL
Gaussian_model  <-  function( d = 1, x0 = 3, probability = TRUE,
                                    n = 1000, r = range(0,10) ) {
    # d is dimension of the parameter space x as well as output space y
    # x0 is a vector of truth observation and max position of exp function
    # n is a number of simulations
    # r is a range for all x
    
    # Define and generate the parameters:
    par.sim  =  data.frame( matrix( NA, nrow = n, ncol = d ) )
    names( par.sim )  =  paste0( 'x', 1:d )
    for( i in 1:d ){
        rnd  =  runif( 10 * n , min = r[1], max = r[2] )
        if ( probability ) { 
            prob  =  1 - exp( x = - ( rnd - x0[ i ] ) ** 2 / 2 )
        } else {
            prob  =  rep( 1, length( rnd ))
        }
        par.sim[ , i ]  =  sample( x = rnd, size = n, 
                                   replace = FALSE, prob = prob )
    }
    
    stat.sim  =  data.frame( matrix( NA, nrow = n, ncol = d ) )
    names( stat.sim )  =  paste0( 'Y', 1:d )
    for( i in 1:d ){
        stat.sim[ , i ]  = exp( x = - ( par.sim[ , i ] - x0[ i ] ) ** 2 / 2 )
    }
    
    stat.obs  =  data.frame( NULL )
    for ( i in 1:d ){
        stat.obs[ 1, paste0( 'Y', i ) ]  =  1
    }
    
    return( list( stat.sim = stat.sim, par.sim = par.sim, stat.obs = stat.obs )  )
}


#' @describeIn Gaussian_model The model of simulations that is based on linear functions for each dimension 
#'
#' @description The function \code{linear_model()} allows to generate 
#' parameters and statistics of simulations that are based on linear function for each dimension: \cr
#' \code{f( x ) = { 1  +  ( x1 - x01 ) / x01 + noise_1, ..., 1  +  ( xn - x0n ) / x0n + noise_n  } = { y1, y2, ..., yn } } is a vector of output data
#' 
#' @param noise Noise factor, implemented as coefficient in \code{f(x) = ... + noise * runif(1)}
#'
#' @return The function \code{linear_model()} returns list of two objects: \cr
#' - stat.sim - data frame of simulations statistics,
#' - par.sim - data frame of parameters.
#' 
#' @export 
#'
#' @examples
#' NULL
linear_model  <-  function( d = 1, x0 = 3, probability = TRUE, noise = 0.2,
                            n = 1000, r = range(0,10) ) {
    # d is dimension of the parameter space x as well as output space y
    # x0 is a vector of truth observation and max position of exp function
    # n is a number of simulations
    # r is a range for all x
    
    # Define and generate the parameters:
    par.sim  =  data.frame( matrix( NA, nrow = n, ncol = d ) )
    names( par.sim )  =  paste0( 'x', 1:d )
    for( i in 1:d ){
        rnd  =  runif( 10 * n , min = r[1], max = r[2] )
        if ( probability ) { 
            prob  =  1 - exp( x = - ( rnd - x0[ i ] ) ** 2 / 2 )
        } else {
            prob  =  rep( 1, length( rnd ))
        }
        par.sim[ , i ]  =  sample( x = rnd, size = n, 
                                   replace = FALSE, prob = prob )
    }
    
    if ( min( x0 ) <= 0 ) stop( 'Each element of x0 vector should be positive. ')
    stat.sim  =  data.frame( matrix( NA, nrow = n, ncol = d ) )
    names( stat.sim )  =  paste0( 'Y', 1:d )
    for( i in 1:d ){
        stat.sim[ , i ]  = 1 + ( par.sim[ , i ] - x0[ i ] ) / x0[ i ]
    }
    
    
    stat.obs  =  data.frame( NULL )
    for ( i in 1:d ){
        stat.obs[ 1, paste0( 'Y', i ) ]  =  1
    }
    
    return( list( stat.sim = stat.sim, par.sim = par.sim, stat.obs = stat.obs )  )
}


