# K2-ABC -------------------------------------------
# The 'kernlab' package is used for standard kernels and Gram matrix calculations:
# for example:
# kernlab::kernelMatrix(kernel,x,y), where x and y are matrices, and kernel is a kernel function:
# - rbfdot Radial Basis kernel function:
# - polydot Polynomial kernel function
# - vanilladot Linear kernel function
# - tanhdot Hyperbolic tangent kernel function
# - laplacedot Laplacian kernel function
# - besseldot Bessel kernel function
# - anovadot ANOVA RBF kernel function
# - splinedot the Spline kernel

#' Function to get parameter estimation and weights using K2-ABC method 
#' 
#' @description \code{K2_ABC()} function allows to get parameter 
#' estimation and weights using K2-ABC method described in 
#'  the paper 
#' Mijung Park, Wittawat Jitkrittum, Dino Sejdinovic, 
#' Proceedings of the 19th International Conference on Artificial Intelligence and Statistics, 
#' [PMLR 51:398-407, 2016](https://proceedings.mlr.press/v51/park16.html). 
#'
#' @param G Gram matrix between sim.stat and obs.stat data sets/matrices that can be obtained by \cr
#' \code{ krnl = rbfdot(sigma = 1) } \cr 
#' \code{G = kernelMatrix( kernel = krnl, x = as.matrix(sim.stat), y = as.matrix(obs.stat) ) }
#' @param epsilon adjust parameter in the weight estimation: \cr
#' \code{w_ij = exp(- ( 1 - k_(y_i, y_obs )/epsilon )}, where \code{y_obs } is observation value
#' and \code{y_i} is a i-th point from sim.stat, by default epsilon = 0.5; \cr
#' also it is possible to adjust epsilon parameter by \code{adjust_K2_ABC()} function 
#' to get the best estimation of a parameter
#' @param par.sim dataset/matrix of parameters for simulation
#' 
#' @return \code{K2_ABC()} returns the list of: \cr 
#' 1) weights for \code{par.sim} related to observation point based on Gram matrix \cr
#' 2) parameter estimation par.est 
#' 
#' @export 
#' 
#' @examples 
#' NULL
K2_ABC  <-  function( G, epsilon = 0.5, par.sim ){
    
    weights  =  as.numeric( exp( (G - 1) / epsilon ) )
    sum_wghts  = sum( weights )
    weights  =  weights / sum_wghts 
    
    par.est  =  par.sim[1, ]
    for( i in 1:ncol(par.est) ){
        par.est[1, i]  =  sum( weights * par.sim[ , i ] )
    }
    
    return( list( weights = weights, par.est = par.est ) )
}


#' @describeIn K2_ABC Function to adjust epsilon parameter for K2-ABC method
#' 
#' @description \code{adjust_K2_ABC()} allows to adjust epsilon parameter for K2-ABC method
#' using numeric vector of epsilon, find parameter estimation for each epsilon and choose the best one.
#'
#' @param epsilon Numeric vector of possible values of epsilon,
#' by default \code{epsilon = (0.05 * 1:20) }
#' @param stat.sim Matrix of statistics of simulations
#' @param stat.obs Matrix of statistics of an observation
#' @param kernel Kernel function of class kernel from kernlab package
#'
#' @return \code{adjust_K2_ABC()} returns the best parameter estimation using K2-ABC method varying epsilon
#' 
#' @export 
#'
#' @examples
#' NULL
adjust_K2_ABC  <-  function(epsilon = c(0.01, 0.02, 0.03, 0.04, (0.05 * 1:20) ), par.sim, stat.sim, stat.obs, kernel ){
    
    ### Get the nearest point to stat.obs
    dst  =  as.matrix(dist( x = rbind( stat.sim, stat.obs ) ) )
    id   =  as.numeric( which.min( dst[ nrow(dst), 1:(nrow(dst)-1)]) )[ 1 ]
    
    ### Get new sets and truth parameter to check hyper parameter epsilon
    x    =  as.matrix( stat.sim[ -id, ] )
    y    =  as.matrix( stat.sim[  id, ]  )
    par.truth  =  as.matrix( par.sim[ id, ] )
    
    ### Get adjusted Gram matrix for x and y
    G = adjust_Gram( kernel, sigma = ( 2**(1:20) ) * 1E-3, x, y )[[ 'G' ]]
    
    ### Get the best epsilon for K2_ABC method based on par.truth
    dlt  =  sapply( epsilon, 
                    FUN = function( eps ) sum( ( par.truth - K2_ABC( G, epsilon = eps, par.sim = par.sim[-id, ] )$par.est   ) ** 2  )
    )
    epsilon  =  epsilon[ which.min( dlt ) ]
    
    G = adjust_Gram( kernel, sigma = ( 2**(1:20) ) * 1E-3, 
                     x = as.matrix( stat.sim ), 
                     y = as.matrix( stat.obs ) )[[ 'G' ]]
    
    K2  =  K2_ABC( G, epsilon = epsilon, par.sim )
    
    return( K2$par.est )
}



#' @describeIn K2_ABC Function to adjust epsilon parameter for K2-ABC method
#' 
#' @description \code{adjust_K2_ABC_iKernel()} allows to adjust epsilon parameter for K2-ABC method
#' using numeric vector of epsilon, find parameter estimation for each epsilon and choose the best one
#' based on matrix of isolation kernel.
#'
#' @param G Kernel matrix \code{G} containers similarities between simulation points and observation point based on isolation kernel
#'
#' @return \code{adjust_K2_ABC_iKernel()} returns the best parameter estimation using K2-ABC method varying epsilon and based on isolation kernel
#' 
#' @export 
#'
#' @examples
#' NULL
adjust_K2_ABC_iKernel  <-  function(epsilon = c(0.01, 0.02, 0.03, 0.04, (0.05 * 1:20) ), par.sim, stat.sim, stat.obs, G ){
    
    ### Get the nearest point to stat.obs
    dst  =  as.matrix(dist( x = rbind( stat.sim, stat.obs ) ) )
    id   =  as.numeric( which.min( dst[ nrow(dst), 1:(nrow(dst)-1)]) )[ 1 ]
    
    ### Get new sets and truth parameter to check hyper parameter epsilon
    x    =  as.matrix( stat.sim[ -id, ] )
    y    =  as.matrix( stat.sim[  id, ]  )
    par.truth  =  as.matrix( par.sim[ id, ] )
    
    ### Get Gram matrix for x and y
    G_xy = G[ -id, ]
    
    ### Get the best epsilon for K2_ABC method based on par.truth
    dlt  =  sapply( epsilon, 
                    FUN = function( eps ) sum( ( par.truth - K2_ABC( G_xy, epsilon = eps, par.sim = par.sim[-id, ] )$par.est   ) ** 2  )
    )
    epsilon  =  epsilon[ which.min( dlt ) ]
    
    # G = adjust_Gram( kernel, sigma = ( 2**(1:20) ) * 1E-3, x = as.matrix( stat.sim ), y = as.matrix( stat.obs ) )
    
    K2  =  K2_ABC( G, epsilon = epsilon, par.sim )
    
    return( K2$par.est )
}



#' Function to get Gram matrix after adjusting of sigma parameter of a kernel 
#'
#' @param kernel Function of kernel with parameter sigma, class from kernel lab package
#' @param sigma numeric vector of possible values of the sigma parameter for a kernel function,
#' by default \code{sigma = ( 2**(1:20) ) * 1E-3}
#' @param x Matrix of stat.sim
#' @param y Matrix of stat.obs
#'
#' @return \code{adjust_Gram} function returns list of \cr 
#' - Gram matrix after adjusting of sigma; \cr
#' - adjusted sigma value.
#' 
#' @export
#'
#' @examples
#' NULL 
adjust_Gram  <-  function( kernel, sigma = ( 2**(1:20) ) * 1E-3, x, y ){
    
    get_dlt = function( kernel, sigma, x, y ){
        krnl  =  kernel( sigma )
        G = kernelMatrix( kernel = krnl, x, y)
        dlt = max( G ) - min( G )
        return( dlt )
    }
    
    dlt  =  sapply( sigma, FUN = function(sgm ) get_dlt(kernel = kernel, sigma = sgm, x, y ) )
    
    sigma  =  sigma[ which.max( dlt ) ][ 1 ]
    krnl  =  kernel( sigma )
    G = kernelMatrix( kernel = krnl, x, y)
    
    return( list(G = G, sigma = sigma ) )
}






# ABC_standard  -----------------------------------------------------------

#' @describeIn K2_ABC Function to adjust tolerance parameter for rejection ABC method
#' 
#' @description \code{adjust_ABC_tolerance()} allows to adjust tolerance parameter for rejection ABC method
#' using numeric vector of tolerance, find parameter estimation for each tolerance and choose the best one.
#'
#' @param tolerance Vector of tolerance values for rejection ABC method to get the best one, 
#' by default \code{tolerance = c(0.001, 0.002, 0.005, (0.01 * 1:20)}
#'
#' @return \code{adjust_ABC_tolerance()} returns the best parameter estimation using rejection ABC method varying tolerance and tolerance value
#' 
#' @export 
#'
#' @examples
#' NULL
adjust_ABC_tolerance  <-  function( tolerance = c(0.001, 0.002, 0.005, (0.01 * 1:20) ), par.sim, stat.sim, stat.obs ){
    
    ### Get the nearest point to stat.obs
    dst  =  as.matrix(dist( x = rbind( stat.sim, stat.obs ) ) )
    id   =  as.numeric( which.min( dst[ nrow(dst), 1:(nrow(dst)-1)]) )[ 1 ]
    
    ### Get new sets and truth parameter to check hyper parameter epsilon
    x    =  as.matrix( stat.sim[ -id, ] )
    y    =  as.matrix( stat.sim[  id, ]  )
    par.truth  =  as.matrix( par.sim[ id, ] )
    
    ### Get the best epsilon for K2_ABC method based on par.truth
    dlt  =  sapply( tolerance, 
                    FUN = function( tlr ) {
                        rej      =  abc( target = y, param = par.sim[-id, ], 
                                         sumstat = x, tol = tlr, method = 'rejection' ) 
                        par.est  =  Get_MAP( as.data.frame( rej$unadj.values ) )
                        
                        # if (nrow( unique.data.frame( l ) ) > 1 ){
                        #    par.est  =  point_estimate( l )$MAP
                        # } else {
                        #     par.est  =  unique.data.frame( l )
                        # }
                        return ( sum( ( par.truth - par.est ) ** 2  ) )
                    }
    )
    tol      =  tolerance[ which.min( dlt ) ][ 1 ]
    
    rej      =  abc( target = stat.obs, param = par.sim, 
                     sumstat = stat.sim, tol = tol, method = 'rejection' ) 
    par.est  =  Get_MAP( as.data.frame( rej$unadj.values ) )
    
    # l = as.data.frame( rej$unadj.values )
    # if (nrow( unique.data.frame( l ) ) > 1 ){
    #     par.est  =  point_estimate( l )$MAP
    # } else {
    #     par.est  =  unique.data.frame( l )
    # }
    
    return( list( par.est   =  par.est, 
                  tolerance =  tol ) )
}



# Get hyperparameters -----------------------------------------------------

Get_hyperparameters  <-  function( par.sim, stat.sim, stat.obs ){
    
    
}

