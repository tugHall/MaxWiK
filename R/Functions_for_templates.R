
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
                                           x = as.data.frame(x) ) ) / A[ i ]
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
                            x = x )
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
                              A = 1, sigma = 1, x ){
    
    par.sim  =   as.numeric( x )
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
                                           x = as.data.frame(x) ) ) / A[ i ]
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
                             x = x )
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
                               A = 1, x ){
    if ( any(x0 == 0 ) )  stop( 'X0 parameters should be non-zeros for Linear model')
    par.sim  =   as.numeric( x )
    sim1  =  data.frame( matrix( NA, nrow = 1, ncol = d ) )
    names( sim1 )  =  paste0( 'Y', 1:d )
    for( i in 1:d ){
        sim1[ 1, i ]  = A[ i ] * ( par.sim[ i ] / x0[ i ]  +
                                       runif( n = 1, min = -0.5, max = 0.5 ) * noise )
    }
    return( sim1 )
}




# MSE is mean square error:
# for parameters (if true parameter is known):


#' The function to get the mean square error values for statistics of simulations
#'
#' @description The function \code{MSE_sim()} allows to get 
#' the mean square error values for statistics of simulations
#'
#' @param stat.sim Summary statistics of the simulations (model output)
#' @param stat.obs Summary statistics of the observation point
#'
#' @return The function \code{MSE_sim()} returns numeric vector of
#' the mean square error values for statistics of simulations
#' 
#' @export
#'
#' @examples
#' NULL
MSE_sim   <-   function( stat.obs, stat.sim ){
    if ( nrow( stat.obs ) > 1 ) stop( 'Please, use stat.obs for each row data independently. 
                                            stat.obs should be single row data.frame')
    # diff between stat.obs and stat.sim
    df   =  sapply( X = 1:ncol( stat.obs ), FUN = function( x )  stat.sim[,x] - stat.obs[1,x]  )
    if ( !is.data.frame( df ) & !is.matrix( df ) ) df = as.data.frame( matrix( df, nrow = 1 ) )
    mse  =  sapply( X = 1:nrow(df) , FUN = function( x ) sum( df[x,] ** 2 ) )
    
    return( mse ) # mean( mse ) )
}


#' @describeIn MSE_sim The function calculates mean square error (MSE) value 
#' for parameters as differences between them and already the known truth parameter 
#'
#' @description The function \code{MSE_parameters()} allows to get MSE for parameters if the truth parameter is known
#'
#' @param par.truth The truth parameter
#' @param par.top Parameters from the top of similarities of \code{iKernelABC()} algorithm
#' @param par.best The best parameter from \code{iKernelABC()} algorithm
#'
#' @return The function \code{MSE_parameters()} returns list of two numbers: \cr
#' - mean of MSE values for all the points from par.top; \cr
#' - MSE value for the point of par.best 
#' 
#' @export
#' 
#' @examples
#' NULL 
MSE_parameters   <-   function( par.truth, par.top = NULL, par.best ){
    names( par.truth )  =  names( par.best )
    if ( nrow(par.best) > 1 ) stop( 'Please, use par.top for multi row data. par.best should be single row data.frame')
    # diff between par.top and par.truth
    if ( !is.null( par.top ) ){
        df  =  sapply( X = 1:ncol( par.truth ), FUN = function( x ) par.top[,x] - par.truth[,x] )
        mse_top  =  sapply( X = 1:nrow(df) , FUN = function( x ) sum( df[x,] ** 2 ) )
    } else { 
        mse_top = NULL 
    }
    
    # diff between par.best and par.truth
    df  =  sapply( X = 1:ncol( par.truth ), FUN = function( x ) ( par.best[,x] - par.truth[,x] ) )
    mse_best  =  sum( df ** 2 ) 
    
    if ( is.null( mse_top ) ) {
        return( list( mse_top = NULL, mse_best = mse_best ) ) 
    } else {
        return( list( mse_top = mean( mse_top ), mse_best = mse_best ) )
    }
}







#' @describeIn Get_call Function to make hypersurface 
#' 
#' @description Function to make 'hypersurface' or data set reduced from the bigger one with
#' selected points around observation point but 
#' the points are away from each other
#' 
#' @param blend_part Part of blended sets: the hypersurface set and 
#' set from rejection ABC 
#'
#' @export
#'
#' @return \code{make_hypersurface()} returns the list: \cr
#' - stat.sim  -  dataset of statistics of simulations;
#' - par.sim  -  dataset of corresponding parameters of simulations.
#' 
make_hypersurface  <-  function(  stat.obs,
                                  stat.sim,
                                  par.sim, 
                                  size, 
                                  blend_part = 0.5 ){
    
    stp  =  round( nrow( par.sim ) / size )
    if ( stp < 2 ) {
        print( 'The size is similar to the dataset size, so return whole dataset.' )
        return( list( stat.sim =  stat.sim,
                      par.sim  =  par.sim  ) )
    }
    ### Make ordering in data sets from close to far points  to observation
    rej      =  abc( target = stat.obs, param = par.sim, 
                     sumstat = stat.sim, tol = 0.1, method = 'rejection' ) 
    w2   =   order(  rej$dist, decreasing = FALSE )
    
    stat_sim  =  stat.sim[ w2, ] 
    par_sim   =  par.sim[ w2,  ]
    row.names( stat_sim )  =  1 : nrow( stat_sim )
    row.names( par_sim  )  =  1 : nrow( par_sim  )
    
    ### distances between all the points:
    dst       =  as.matrix( dist( x = stat_sim ) )
    
    # Choose the farthest points from the previous points 
    #######  to get maximal volume of hypersurface:
    for( i in 1 : size ){
        
        if (i == 1 ) {
            w = stp
            next
        }
        
        dstncs  =  unlist( 
            sapply( X    =  ( 1 : (i*stp) )[ -w ], 
                    FUN  =  function( x ){
                        sum( dst[ w , x ]  )
                    } )
        )
        w_m   =   which.max( dstncs )[ 1 ]
        w_mm  =   ( 1 : (i*stp) )[ -w ][ w_m ]
        w     =   c( w, w_mm )
        
    }
    
    # Define indexes of dataset for output:
    w3   =  (1:nrow( par_sim  ) )[ -w ]
    prt  =  round( blend_part * size )
    w3   =  w3[ 1 : prt ]
    w    =  w[ 1  : (size  -  prt) ]
    
    w    =  c( w, w3 )
    
    # Save to output:
    stat.sim_red  =  stat_sim[ w, ]
    par.sim_red   =  par_sim[  w, ]
    row.names( par.sim_red  )  =  1:nrow( par.sim_red  )
    row.names( stat.sim_red )  =  1:nrow( stat.sim_red )
    
    
    return( list( stat.sim =  stat.sim_red,
                  par.sim  =  par.sim_red   ) )
}









