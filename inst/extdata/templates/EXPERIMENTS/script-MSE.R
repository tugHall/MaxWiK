
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


