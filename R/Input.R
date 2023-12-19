
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


