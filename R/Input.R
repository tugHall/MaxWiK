
#' Function to get dataset of simulations that is based on Gaussian functions for each dimension 
#'
#' @description The function \code{get_dataset_of_Gaussian_model()} allows to generate 
#' parameters and statistics of simulations that are based on Gaussian function for each dimension: \cr
#' \code{f( x ) = { exp( - (x1-x01)^2 / 2 ), ..., exp( - (xn-x0n)^2 / 2 )  } = { y1, y2, ..., yn } } is a vector of output data
#' 
#'
#' @param d Dimension of the parameter and model space
#' @param x0 Numeric vector with length of dimensionality of data frame, 
#' that contents the truth value of parameter. Each number in the vector should be within the range r
#' @param probability Logical, if TRUE then apply uneven distribution for parameters generation
#' @param n Integer number of points in data frame
#' @param r Range \code{r = c(min, max)}, by default  \code{r = range(0,10)} for all the dimensions
#' @param noise Number corresponding to noise level in the dataset. 0 means the absent of noise, 
#' 1 means noise level is the same as amplitude, more than 1 means that a noise is greater than 'signal'. 
#' @param A Vector of exponent multiplication factors
#' @param sigma Vector of sigmas in Gaussian function. Length of vector should be equal d.
#'
#' @return The function \code{get_dataset_of_Gaussian_model()} returns list of three objects: \cr
#' - stat.sim - data frame of simulations statistics,
#' - par.sim - data frame of parameters,
#' - stat.obs - data frame of an observation point.
#' 
#' @export
#' 
#' @examples
#' NULL
get_dataset_of_Gaussian_model  <-  function( d = 1, x0 = 3, r = range(0,10), 
                                             noise = 0, A = 1, sigma = 1, 
                                             probability = TRUE, n = 1000 ) {

    # Define and generate the parameters:
    par.sim  =  data.frame( matrix( NA, nrow = n, ncol = d ) )
    names( par.sim )  =  paste0( 'X', 1:d )
    for( i in 1:d ){
        rnd  =  runif( 10 * n , min = r[1], max = r[2] )
        if ( probability ) { 
            prob  =  lapply( rnd, FUN = function( x ){ 
                ( A[ i ] - Gauss_function( d = 1, x0 = x0[ i ], 
                                           r = r, noise = 0, A = A[ i ], sigma = sigma[ i ], 
                                           par.sim1 = as.data.frame(x) ) ) / A[ i ]
                            } )
            prob  =  unlist( prob )
        } else {
            prob  =  rep( 1, length( rnd ))
        }
        par.sim[ , i ]  =  sample( x = rnd, size = n, 
                                   replace = FALSE, prob = prob )
    }
    
    stat.sim  =  data.frame( matrix( NA, nrow = n, ncol = d ) )
    names( stat.sim )  =  paste0( 'Y', 1:d )
    for( i in 1:d ){
        stat.sim[ , i ]  = unlist( lapply( X = par.sim[ , i ], FUN = function( x ){ 
            Gauss_function( d = 1, x0 = x0[ i ], 
                                       r = r, noise = noise, A = A[ i ], sigma = sigma[ i ], 
                                       par.sim1 = x )
                            } )   )
            # A * exp( x = - ( par.sim[ , i ] - x0[ i ] ) ** 2 / 2 ) + 
            # runif( n = n, min = -0.5, max = 0.5 ) * noise
    }
    
    stat.obs  =  data.frame( NULL )
    for ( i in 1:d ){
        stat.obs[ 1, paste0( 'Y', i ) ]  =  A[ i ]
    }
    
    return( list( stat.sim = stat.sim, par.sim = par.sim, stat.obs = stat.obs )  )
}


#' @describeIn get_dataset_of_Gaussian_model Function to get output of Gaussian function for d dimension case
#'
#' @description The function \code{Gauss_function()} is used as 
#' a model based on Gaussian function for multidimensional case.
#' 
#' @return \code{Gauss_function()} returns output of Gaussian function for d dimensions 
#' 
#' @export
#'
#' @examples
#' NULL
Gauss_function  <-  function( d = 1, x0 = 3, r = range(0,10), noise = 0, 
                              A = 1, sigma = 1, par.sim1 ){
    
    par.sim  =   as.numeric( par.sim1 )
    sim1  =  data.frame( matrix( NA, nrow = 1, ncol = d ) )
    names( sim1 )  =  paste0( 'Y', 1:d )
    for( i in 1:d ){
        sim1[ 1, i ]  = A[ i ] * ( exp( x = - ( par.sim[ i ] - x0[ i ] ) ** 2 / 2 / sigma[ i ] / sigma[ i ] ) + 
            runif( n = 1, min = -0.5, max = 0.5 ) * noise )
    }
    
    return( sim1 )
}



#' Function to get dataset of simulations that is based on Linear functions for each dimension 
#'
#' @description The function \code{get_dataset_of_Linear_model()} allows to generate 
#' parameters and statistics of simulations that are based on Linear function for each dimension: \cr
#' \code{f( x ) = { A1 * x1 / x01, ..., An * xn / x0n  } = { y1, y2, ..., yn } } is a vector of output data
#' 
#'
#' @param d Dimension of the parameter and model space
#' @param x0 Numeric vector with non-zero values. The vector's length should be equal d.  
#' \code{x0} contents the truth value of parameter. Each number in the vector should be within the range r
#' @param probability Logical, if TRUE then apply uneven distribution for parameters generation
#' @param n Integer number of points in data frame
#' @param r Range \code{r = c(min, max)}, by default  \code{r = range(0,10)} for all the dimensions
#' @param A Vector of multiplication factors.
#'
#' @return The function \code{get_dataset_of_Linear_model()} returns list of three objects: \cr
#' - stat.sim - data frame of simulations statistics,
#' - par.sim - data frame of parameters,
#' - stat.obs - data frame of an observation point.
#' 
#' @export
#' 
#' @examples
#' NULL
get_dataset_of_Linear_model  <-  function( d = 1, x0 = 3, r = range(0,10), 
                                             noise = 0, A = 1, 
                                             probability = TRUE, n = 1000 ) {
    
    # Define and generate the parameters:
    par.sim  =  data.frame( matrix( NA, nrow = n, ncol = d ) )
    names( par.sim )  =  paste0( 'X', 1:d )
    for( i in 1:d ){
        rnd  =  runif( 10 * n , min = r[1], max = r[2] )
        if ( probability ) { 
            prob  =  lapply( rnd, FUN = function( x ){ 
                ( A[ i ] - Gauss_function( d = 1, x0 = x0[ i ], 
                                           r = r, noise = 0, A = A[ i ], sigma = (r[2]-r[1])/5, 
                                           par.sim1 = as.data.frame(x) ) ) / A[ i ]
            } )
            prob  =  unlist( prob )
        } else {
            prob  =  rep( 1, length( rnd ))
        }
        par.sim[ , i ]  =  sample( x = rnd, size = n, 
                                   replace = FALSE, prob = prob )
    }
    
    stat.sim  =  data.frame( matrix( NA, nrow = n, ncol = d ) )
    names( stat.sim )  =  paste0( 'Y', 1:d )
    for( i in 1:d ){
        stat.sim[ , i ]  = unlist( lapply( X = par.sim[ , i ], FUN = function( x ){ 
            Linear_function( d = 1, x0 = x0[ i ], 
                            r = r, noise = noise, A = A[ i ], 
                            par.sim1 = x )
        } )   )
        # A * exp( x = - ( par.sim[ , i ] - x0[ i ] ) ** 2 / 2 ) + 
        # runif( n = n, min = -0.5, max = 0.5 ) * noise
    }
    
    stat.obs  =  data.frame( NULL )
    for ( i in 1:d ){
        stat.obs[ 1, paste0( 'Y', i ) ]  =  A[ i ]
    }
    
    return( list( stat.sim = stat.sim, par.sim = par.sim, stat.obs = stat.obs )  )
}


#' @describeIn get_dataset_of_Linear_model Function to get output of Linear function for d dimension case
#'
#' @description The function \code{Linear_function()} is used as 
#' a model based on Linear function for multidimensional case.
#' 
#' @return \code{Linear_function()} returns output of Linear function for d dimensions 
#' 
#' @export
#'
#' @examples
#' NULL
Linear_function  <-  function( d = 1, x0 = 3, r = range(0,10), noise = 0, 
                              A = 1, par.sim1 ){
    if ( any(x0 == 0 ) )  stop( 'X0 parameters should be non-zeros for Linear model')
    par.sim  =   as.numeric( par.sim1 )
    sim1  =  data.frame( matrix( NA, nrow = 1, ncol = d ) )
    names( sim1 )  =  paste0( 'Y', 1:d )
    for( i in 1:d ){
        sim1[ 1, i ]  = A[ i ] * ( par.sim[ i ] / x0[ i ]  +
                                  runif( n = 1, min = -0.5, max = 0.5 ) * noise )
    }
    return( sim1 )
}



#' @describeIn iKernelABC Function to restrict data in the size to accelerate the calculations 
#' 
#' @description \code{restrict_data()} is based on rejection ABC method to restrict original dataset
#'
#' @param size Integer number of points to leave from original dataset
#'
#' @return \code{restrict_data()} returns the list of: \cr
#' par.sim - restricted parameters which are close to observation point \cr
#' stat.sim - restricted stat.sim which are close to observation point
#' 
#' @export
#'
#' @examples
#' NULL
restrict_data  <-  function( par.sim, stat.sim, stat.obs, size = 300 ){
    l  =  nrow( par.sim )
    if ( l != nrow( stat.sim ) ) stop( "The parameters and statistics of simulations have different number of rows" )
    rej_abc  =  abc( target = stat.obs, param = par.sim, sumstat = stat.sim, 
                     tol = size / l, method = 'rejection' )
    # rej_abc$region
    new_par  =   par.sim[ rej_abc$region, ]
    new_sim  =  stat.sim[ rej_abc$region, ]
    return( list( par.sim   =  new_par, 
                  stat.sim  =  new_sim ) )
}


