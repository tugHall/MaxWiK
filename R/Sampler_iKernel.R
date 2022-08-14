# This is a sampler for ABC method that based on maxima weighted isolation kernel method

#' Model definition for sampler to calculate one simulation
#'
#' @param name Name of a model, can be either \code{ 'Gaussian' or 'Linear' }
#' @param parameter Value of a parameter
#' @param x0 True value of parameter
#' @param stat.obs Data frame of statistics of observation point
#' @param noise Value of stochastic term 
#'
#' @return \code{model()} returns data frame with a result of a simulation based on the model
#' 
#' @export
#'
#' @examples
#' NULL
model  <-  function( name = c( 'Gaussian', 'Linear' )[1], 
                     parameter, x0, stat.obs, noise = 0 ){
    
    d = length( x0 )
    sim = as.data.frame( matrix( data = 0, nrow = 1, ncol = d ) ) 
    names( sim )  =  names( stat.obs )
    
    if ( name == 'Gaussian' ){
        for( i in 1:d ){
            sim[ 1, i ]  = exp( x = - ( parameter[1 , i ] - x0[ i ] ) ** 2 / 2 ) + noise * runif(1)
        }
    } else {
        if ( name == 'Linear' ){
            for( i in 1:d ){
                sim[ 1, i ]  =  1  +  ( parameter[1 , i ] - x0[ i ] ) / x0[ i ] + noise * runif(1)
            }
        } else stop( 'This function is not defined' )
    }    
    
    return( sim )
}


#' @describeIn iKernelABC Function to restrict data in the size to accelerate the calculations 
#' 
#' @description \code{restricrt_data()} is based on rejection ABC method to restrict original dataset
#'
#' @param size Integer number of points to leave from original dataset
#'
#' @return \code{restricrt_data()} returns the list of: \cr
#' par.sim - restricted parameters which are close to observation point \cr
#' stat.sim - restricted stat.sim which are close to observation point
#' 
#' @export
#'
#' @examples
#' NULL
restricrt_data  <-  function( par.sim, stat.sim, stat.obs, size = 300 ){
    l  =  nrow( par.sim )
    if ( l != nrow( stat.sim ) ) stop( "The parameters and statistics of simulations have different number of rows" )
    rej_abc  =  abc( target = stat.obs, param = par.sim, sumstat = stat.sim, 
                     tol = size / l, method = 'rejection' )
    # rej_abc$region
    
    return( list( par.sim   =   par.sim[ rej_abc$region, ], 
                  stat.sim  =  stat.sim[ rej_abc$region, ] ) )
}



#' @describeIn iKernelABC Function to generate parameters and simulate a model based on MaxWiK algorithm 
#'
#' @param model is a function to get output of simulation during sampling 
#' @param arg0 is a list with arguments for a model function, so that arg0 is NOT changed during sampling
#' @param size Number of point to restrict original dataset
#' @param nmax is maximal number of iterations
#' @param include_top Logical to include top points from \code{spider_web()} function to simulate or do not
#' @param slowly Logical for two algorithms: slow and fast seekers in sampling
#' @param rate Rate value in the range \code{[0,1]} to define 
#' the rate of changing in the original top of sampled points 
#'
#' @return \code{sampler_MaxWiK()} returns the list: \cr
#' results - results of simulations; \cr 
#' best - the best value of parameter; \cr
#' MSE_min - minimum of MSE; \cr
#' number_of_iterations - number of iterations; \cr
#' time - time of sampling in seconds.
#' 
#' @export
#'
#' @examples
#' NULL
sampler_MaxWiK  <-  function( stat.obs, stat.sim, par.sim, model, 
                                             arg0 = list(),  size = 500, 
                                             psi_t, epsilon, nmax = 100, 
                                             include_top = FALSE,
                                             slowly = FALSE, rate = 0.2 ){ 
    # epsilon is a criterion to stop simulation
    # nmax is maximal number of iterations
    # psi_t is a data.frame of psi and t with top values of MSE of length 20
    # size is a size of data frame for a single simulation
    # i  =  sample( 1:nrow( psi_t ), 1, replace = TRUE, prob = psi_t$prob )
    
    # model is a function to get output of simulation during sampling 
    # arg0 is a list with arguments for a model function, so that arg0 is NOT changed during sampling
     
    start_time = Sys.time()
    # Redefine input:
    stat_sim  =  stat.sim     #  initial data for summary statistics of simulations
    par_sim   =  par.sim      #  initial data for parameters of simulations
    stat_obs  =  stat.obs     #  
    
    # Get dimensionalities of parameter and simulation datasets:
    dim_par  =  ncol( par_sim )
    dim_sim  =  ncol( stat_sim )
    
    results  =  data.frame( NULL )
    results_ALL  =  results
    combine  =  data.frame( NULL )
    err_previous  =  1E3
    pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                         max = nmax, # Maximum value of the progress bar
                         style = 3,    # Progress bar style (also available style = 1 and style = 2)
                         width = 50,   # Progress bar width. Defaults to getOption("width")
                         char = ":")   # Character used to create the bar
    
    for( i in 1:nmax ){
        
        # Progress BAR
        setTxtProgressBar(pb, i)
        
        results  =  data.frame( NULL )
        
        ### Get new parameters and corresponding results of simulations based on each set of psi and t
        # source( './lib_iKernel.R' )
        for ( j in 1:10 ){
            
            web = spiderweb( param = par_sim, stat.sim = stat_sim, stat.obs = stat_obs, 
                             psi = psi_t$psi[ j ], t = psi_t$t[ j ], talkative = FALSE )
            
            res_1  =  web$par.best
            if ( include_top ){ 
                res_1  =  rbind( res_1, web$par.top )
            }
            
            new_best.sim  =  do.call( what = model, args = c( arg0, list( parameter =  web$par.best ) ) )
            if ( include_top ){
                new_top.sim  =  NULL
                for( p in 1:nrow( web$par.top ) ){
                    new_top_1.sim  =  do.call( what = model, 
                                               args = c( arg0, list( parameter =  web$par.top[ p, ] ) ) )
                    new_top.sim  =  rbind( new_top.sim, new_top_1.sim ) 
                }
                
                res_1[ , (ncol(res_1) + 1) : (ncol(res_1) + ncol(new_top.sim) ) ]  =  
                    rbind( new_best.sim, new_top.sim )
            } else {
                res_1[ , (ncol(res_1) + 1) : (ncol(res_1) + ncol(new_best.sim) ) ]  =  
                    new_best.sim
            }
            
            mse_1   =  MSE_sim( stat.obs = stat.obs, stat.sim = new_best.sim )
            mse_1   =  as.data.frame( mse_1 )
            names( mse_1)  = 'mse'
            if ( include_top ){ 
                mse_1[2:20, 1]  = MSE_sim( stat.obs = stat.obs, stat.sim = new_top.sim )
            }
            
            res_1$mse    =  as.numeric( mse_1$mse )
            
            if ( include_top ){ 
                res_1$comm     =  c( 'Best', rep('Top', 19 ) )
            } else {
                res_1$comm     =     'Best'
            }
            if ( include_top ){
                res_1$sim      =  c( web$sim.best, web$sim.top )
            } else{
                res_1$sim      =     web$sim.best
            }
            results  =  rbind( results, res_1 )
        }                
        
        results_ALL  =  rbind( results_ALL, results )
        
        ### Combine results with stat_sim and par_sim:
        row_res  =  nrow( results )
        row_sim  =  nrow( stat_sim )
        
        combine  =  results
        combine[ (row_res+1):(row_res+row_sim), 1:dim_par  ]  =  par_sim
        combine[ (row_res+1):(row_res+row_sim), (dim_par+1):(dim_par+dim_sim)  ]  =  stat_sim
        combine[ (row_res+1):(row_res+row_sim), 'comm'  ]  =  'previous'
        combine[ (row_res+1):(row_res+row_sim), 'sim'   ]  =  web$iKernelABC$similarity
        combine[ (row_res+1):(row_res+row_sim), 'mse'   ]  =  MSE_sim( stat.obs = stat_obs, stat.sim = stat_sim )    
        # Make the sorting of data regarding MSE:
        combine   =   unique.data.frame( combine )
        combine   =   combine[ order( combine$mse , decreasing = FALSE), ]
        row.names( combine ) = 1:nrow( combine )
        # Define new dataset for next step with a new sampling:
        
        # slct      =  c( 1:round( 0.8*size ), ( nrow( combine ) - round( 0.2* size) + 1 ):nrow( combine ) )
        if ( !slowly ){
            slct      =  c( 1:round( size ) )
            par_sim   =  combine[ slct, 1:dim_par  ]
            stat_sim  =  combine[ slct, (dim_par+1):(dim_par+dim_sim)  ]
        } else {
            # if SLOWLY == TRUE then add new best points with slow rate 
            w  =  order( as.numeric( web$iKernelABC$similarity ) )[ 1:round( size*rate) ] 
            slct   =   1:length( w )
            
            par_sim[ w, ]   =  combine[ slct, 1:dim_par  ]
            stat_sim[ w, ]  =  combine[ slct, (dim_par+1):(dim_par+dim_sim)  ]
        }
        
        err       =  combine$mse[ 1 ]  #  minimum of MSE in accordance with sorting
        if ( abs(err - err_previous) < epsilon ) break 
        err_previous  =  err
    }
    
    best  =  results[ which.min(results$mse ) , ]
    end_time = Sys.time()
    
    #Stop progress BAR:
    close(pb)
    
    return( list( results = results_ALL, best = best, MSE_min = err, number_of_iterations = i, time = end_time - start_time ) )
}


if ( FALSE ){
    smpl_1  =  sampler_MaximaWeightediKernel( stat.obs, stat.sim, par.sim,  model = model, 
                                         arg0 = list(  name = c( 'Gauss', 'Linear' )[1],
                                                       x0 = x0, 
                                                       stat.obs = stat.obs, 
                                                       noise = 0 ), 
                                         size = 200, psi_t, epsilon = 1E-12, 
                                         nmax = 30, include_top = TRUE,
                                         slowly = TRUE, rate = rat[ r ] )
    if ( r == 1 ){ smpl = smpl_1 } else { smpl = Map( c, smpl, smpl_1 ) }
}





