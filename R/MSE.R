
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


#' The function to get Maximum A Posteriori from numeric data frame
#' 
#' @description The function \code{Get_MAP} returns the Maximum A Posteriori (MAP) of data frame.
#' 
#' @param DF Data frame of the integer or numeric numbers or characters 
#'
#' @return The function \code{Get_MAP} returns the MAP for each dimension in data frame DF of vector
#' 
#' @export
#' 
#' @examples
#' NULL 
Get_MAP <- function( DF ) {
    
    # Check Data frame:
    if ( !is.data.frame( DF ) ) {
        print('Input is ')
        print( DF )
        print('str of input:')
        print( str(DF) )
        stop( 'The input of Get_MAP function is NOT a data frame' )
    }
    if ( nrow( DF ) == 0 ){
        print( DF )
        stop(' Data frame is empty ')
    }
    if ( any( is.na.data.frame( DF ) ) ){
        print( DF )
        stop( 'Data frame has NA value')
    }
    
    # Define the dimension and format of result with data frame:
    d = ncol( DF )
    res = data.frame( matrix( data = NA, ncol = d ) )
    names( res )  =  names( DF )
    
    if ( nrow( DF ) > 1 ){
        for( i in 1:d ){
            # Get vector from data frame:
            v  =  DF[ , i ]
            # Check that is numeric vector:
            if ( !is.numeric( v ) ){
                print( v )
                print( DF )
                stop( paste( 'Column ', i, 'in data frame is NOT numeric '))
            }
            if ( length( unique( v ) ) < 2 ){
                res[ 1, i ] = unique( v )
            } else {
                res[ 1, i ] =  point_estimate( v + runif( n = length( v ) ) * 1E-20 )$MAP
            }
        }
    } else {
        res  =  DF
    }
    
    return(  res )
}


