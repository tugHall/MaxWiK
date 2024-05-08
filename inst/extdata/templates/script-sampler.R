#' @describeIn get.iKernelABC  function to get parameter estimation based on isolation kernel
#'
#' @param iKernelABC Result of function \code{get.iKernelABC}
#'
#' @return \code{Get_iKernel_estimation()} returns list of: \cr
#' - iKernel_ABC - parameter estimation based on isolation kernel / weighted sum; \cr
#' - K2_ABC_iKernel - parameter estimation based on K2-ABC method with matrix of isolation kernel.
#' 
#' @export
#'
#' @examples
#' NULL
Get_iKernel_estimation  <-  function( iKernelABC, par.sim, stat.sim, stat.obs ){
    
    ### Isolation kernel estimation:
    weights  =  iKernelABC$similarity
    sum_wghts  = sum( weights )
    weights  =  weights / sum_wghts 
    
    par.est  =  par.sim[1, ]
    for( i in 1:ncol(par.est) ){
        par.est[1, i]  =  sum( weights * par.sim[ , i ] )
    }
    
    ### K2-ABC based on iKernel:
    G = matrix( iKernelABC$similarity, ncol = 1 )
    
    K2_ABC_iKernel  =  adjust_K2_ABC_iKernel(epsilon = c(0.01, 0.02, 0.03, 0.04, (0.05 * 1:20) ), 
                                             par.sim, stat.sim, stat.obs, G )[[ 'par.est' ]]
    
    return( list (K2_ABC_iKernel  =  K2_ABC_iKernel, 
                  iKernel_ABC     =  par.est ) )
}


#' @describeIn get.iKernelABC Function to adjust hyper parameters \code{psi} and \code{t} for isolation kernel ABC
#'
#' @param psi_t Initial data.frame of  \code{psi} and \code{t}, by default \cr
#' \code{psi_t = data.frame( psi = as.numeric( sapply( X = c(2:8)*2, FUN = function( x ) rep(x, 8) ) ), t = rep( c(4,6,8,10,12,14,16,20), 7) )}
#' @param cores Number of available cores for parallel computing
#' @param n_best Number of the best adjusted values of psi_t pairs regarding to MSE 
#' @param par.sim Data frame of parameters 
#' 
#' @return \code{adjust_psi_t() } returns adjusted hyper parameters \code{psi} and \code{t} as a data.frame with set of pair \code{psi_t} 
#' 
#' @export
#'
#' @examples
#' NULL
adjust_psi_t  <-  function(par.sim, stat.sim, stat.obs, talkative = FALSE, check_pos_def = FALSE, n_best = 10,
                           psi_t = data.frame( psi = as.numeric( sapply( X = c(2:8)*2, FUN = function( x ) rep(x, 8) ) ), 
                                               t = rep( c(4,6,8,10,12,14,16,20), 7) ),
                           cores = 4 ){
    
    get_dlt  =  function( psi, t, par.sim, stat.sim, stat.obs, par.truth,
                          talkative = talkative, check_pos_def = check_pos_def  ){
        
        iKernABC  =  get.iKernelABC( psi = psi, t = t, param = par.sim, stat.sim = stat.sim, 
                                        stat.obs = stat.obs, talkative = talkative, 
                                        check_pos_def = check_pos_def )
        ### Isolation kernel estimation:
        wghts  =  as.numeric( iKernABC$similarity )
        sum_wghts  = sum( wghts )
        wghts  =  wghts / sum_wghts 
        
        par.est  =  par.sim[1, ]
        for( i in 1:ncol(par.est) ){
            par.est[1, i]  =  sum( wghts * par.sim[ , i ] )
        }
        
        dlt  =  sum( ( par.truth - par.est ) ** 2  )
        
        return( dlt )
    }
    
    ### Get the nearest point to stat.obs and it's id
    dst  =  as.matrix( dist( x = rbind( stat.sim, stat.obs ) ) )
    id   =  as.numeric( which.min( dst[ nrow(dst), 1:(nrow(dst)-1)]) )[ 1 ]
    
    ### Get new sets and truth parameter to check hyper parameters psi_t:
    x    =  stat.sim[ -id, ] 
    y    =  stat.sim[  id, ] 
    par.truth  =  par.sim[ id, ]
    
    psi_t$dlt = 0
    
    dlt  =  mclapply( 1:nrow( psi_t ), FUN = function(i){
        get_dlt( psi = psi_t$psi[ i ], 
                 t = psi_t$t[ i ], 
                 par.sim  = par.sim[-id, ], 
                 stat.sim = x, 
                 stat.obs = y, 
                 par.truth = par.truth,
                 talkative = talkative, check_pos_def = check_pos_def  )}, 
        mc.cores = cores )
    
    for( i in 1:nrow( psi_t ) ){
        if(FALSE){
            dlt  =  get_dlt( psi = psi_t$psi[ i ], 
                             t = psi_t$t[ i ], 
                             par.sim  = par.sim[-id, ], 
                             stat.sim = x, 
                             stat.obs = y, 
                             par.truth = par.truth,
                             talkative = talkative, check_pos_def = check_pos_def  )
        }
        psi_t$dlt[ i ]  =  dlt[[ i ]] 
    }
    
    rdr  =  order( psi_t$dlt, decreasing = TRUE )[1:n_best]
    
    return( psi_t[rdr, c( 1,2 ) ] )
}



# GET_Spider_MAP ---------------------------------------------------------


#' Function to get MAP of SpiderWeb algorithm based on different psi / t hyperparameters
#'
#' @param cores Number of cores for parallel calculation
#' @param par.sim - data frame of parameters of the model
#' @param stat.sim Summary statistics of the simulations (model output)
#' @param stat.obs Summary statistics of the observation point
#' @param restrict_points_number Maximal number of points in the data sets to get MAP
#' @param n_best Number of best psi_t pairs adjusted for MaxWiK algorithm
#'
#' @return Maximum A Posteriori of meta-sampling distribution of parameters
#' 
#' @export
#'
#' @examples
#' NULL
#' # it takes a time for a simulation and then it will demonstrates results, \cr
#' # so, please, wait for a while
#' sim = simulation_example_many_psi_t( verbose = FALSE , to_plot = FALSE )
get_Spider_MAP  <-  function( stat.sim, par.sim, stat.obs, 
                              restrict_points_number = 300, cores = 4,
                              n_best = 8 ){
    
    tol = restrict_points_number / nrow( stat.sim ) 
    
    rej = abc::abc( target = stat.obs, param = par.sim, sumstat = stat.sim,
                    method = 'rejection', tol = tol )
    
    stat.sim  =  stat.sim[ rej$region, ]
    par.sim   =   par.sim[ rej$region, ]    
    
    psi_t  =  adjust_psi_t( par.sim = par.sim, stat.sim = stat.sim, stat.obs = stat.obs,
                            n_best = n_best )
    
    webnet  =  list( )
    SIM  =  function( j ){
        psi =  psi_t$psi[ j ]
        t   =  psi_t$t[ j ]
        web  =  meta_sampling( psi = psi, t = t, param = par.sim, 
                                stat.sim = stat.sim, stat.obs = stat.obs, 
                                talkative = FALSE, check_pos_def = FALSE ,
                                n_bullets = 10, n_best = 20, halfwidth = 0.5, 
                                epsilon = 0.001 )
        return( web )
    } 
    
    webnet  =  mclapply( 1:nrow( psi_t ) , FUN = SIM, mc.cores = cores )
    
    simnet  =  list(  stat.obs  =  stat.obs, 
                      stat.sim  =  stat.sim,
                      par.sim   =  par.sim,
                      psi_t     =  psi_t,
                      webnet  = webnet )
    
    simnet$networks  =  get_network_from_simnet( simnet = simnet )
    
    MAP  =  Get_MAP( as.data.frame( simnet$networks ) )  # point_estimate( simnet$networks )$MAP
    
    return( MAP )
}






#' The function returns the weighted mean of the parameter based on Isolation Kernel ABC method
#'
#' @description The function \code{Mean_iKernel_parameters()} returns the weighted mean of the parameter 
#' that was calculated with \code{KernelABC()} function that represents Isolation Kernel ABC method 
#' 
#' @param param Data frame of parameters
#' @param sm Numeric vector of weights gotten from \code{get.iKernelABC()} function
#' 
#' @return The function \code{Mean_iKernel_parameters()} returns the weighted mean of the parameter 
#' 
#' @examples 
#' NULL 
Mean_iKernel_parameters  <- function( param, sm ){
    Mean_iKrnl = rep( 0, length( param ) )
    sum_sm  =  sum( sm ) 
    for (i in 1:length( param )) {
        Mean_iKrnl[i] = sum( sm * param[,i]) / sum_sm 
    }
    names(Mean_iKrnl) <- names( param )
    return( Mean_iKrnl )
}


