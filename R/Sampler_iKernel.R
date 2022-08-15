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


#' Function to call a method to get parameter estimation and MSE for each used method
#'
#' @param method_name Name of a method
#' @param kernel_name Name of kernel function
#' @param stat.obs Data frame of statistics of observation point
#' @param stat.sim Data frame of statistics of simulations
#' @param par.sim Data frame of parameters
#' @param G Matrix of similarities for K2-ABC based on isolation kernel
#' @param par.truth Truth parameter value to check result of estimation
#'
#' @return \code{par.truth,()} returns the list: \cr
#' method_name = method_name, \cr
#' - kernel_name, \cr
#' - model_name, \cr
#' - stochastic_term, \cr
#' - MSE, \cr
#' - running_time, \cr
#' - iteration.
#' @export
#' 
#' @examples
#' NULL
Get_parameter  <-  function( method_name, kernel_name = '',  
                        stat.obs, stat.sim, par.sim, G = NULL, par.truth ){
    
    n_min  =  100 
    time_start  =  Sys.time( )
    par.est  =  NA 
    
    if ( method_name == 'K2-ABC' & kernel_name != 'iKernel' ){
        if ( kernel_name == 'Gaussian'){
            kernel  =  rbfdot
        } else {
            kernel  =  laplacedot
        }
        par.est  =  adjust_K2_ABC( par.sim = par.sim, stat.sim = stat.sim, 
                                   stat.obs = stat.obs, kernel = kernel )
    }
    
    if ( method_name == 'K2-ABC' & kernel_name == 'iKernel' ){
        par.est  =  adjust_K2_ABC_iKernel( par.sim = par.sim, stat.sim = stat.sim, 
                                           stat.obs = stat.obs, G = G )
    }
    
    if ( method_name == 'Rejection' ){
        res      =  adjust_ABC_tolerance( par.sim = par.sim, stat.sim = stat.sim, 
                                          stat.obs = stat.obs )
        par.est  =  res$par.est
    }
    
    if ( method_name == 'Loclinear' ){
        
        res      =  adjust_ABC_tolerance( par.sim = par.sim, stat.sim = stat.sim, 
                                          stat.obs = stat.obs )
        tol   =   res$tolerance
        if ( tol * nrow(stat.sim) < n_min ) tol = n_min / nrow( stat.sim )
        
        loclin   =  abc(   target = stat.obs, param = par.sim, sumstat = stat.sim, 
                           tol = tol, method  =  'loclinear', hcorr   =  FALSE, 
                           transf=c( "none" ) )
        
        lc = round( loclin$adj.values + runif( n= length(loclin$adj.values) )/1E6, digits = 9 )
        
        if (nrow( unique.data.frame(lc) ) > 1 ){
            par.est  =  point_estimate( lc )$MAP
        } else {
            par.est  =  unique.data.frame(lc)
        }
    }
    
    if ( method_name == 'Neuralnet' ){
        
        res      =  adjust_ABC_tolerance( par.sim = par.sim, stat.sim = stat.sim, 
                                          stat.obs = stat.obs )
        tol   =   res$tolerance
        if ( tol * nrow(stat.sim) < n_min ) tol  =  n_min / nrow( stat.sim )
        
        nn   =  abc(   target = stat.obs, param = par.sim, sumstat = stat.sim, 
                       tol = tol, method  =  'neuralnet', hcorr   =  TRUE, 
                       transf=c("none","log"), lambda = 0.0001, trace = FALSE )
        
        nn_adj = round( nn$adj.values + runif( n= length( nn$adj.values ) ) / 1E6, digits = 9 )
        par.est  =  point_estimate( nn_adj )$MAP
    }
    
    if ( method_name == 'Ridge' ){
        
        res      =  adjust_ABC_tolerance( par.sim = par.sim, stat.sim = stat.sim, 
                                          stat.obs = stat.obs )
        tol   =   res$tolerance
        if ( tol * nrow(stat.sim) < n_min ) tol = n_min / nrow( stat.sim )
        
        rdg   =  abc(   target = stat.obs, param = par.sim, sumstat = stat.sim, 
                        tol = tol, method  =  'ridge', hcorr   =  FALSE, 
                        transf=c("none","log"), kernel = 'epanechnikov' )
        
        rdg_adj = round( rdg$adj.values + runif( n= length( rdg$adj.values ) )/1E6, digits = 9 )
        
        if ( nrow( unique.data.frame( rdg_adj ) ) > 1 ){
            par.est  =  point_estimate( rdg_adj )$MAP
        } else {
            par.est  =  unique.data.frame(rdg_adj)
        }
    }
    
    if ( method_name == 'MaxWiK_MAP' ){
        par.est  =  get_Spider_MAP( par.sim = par.sim, stat.sim = stat.sim, stat.obs = stat.obs )
    }
    
    if ( method_name == 'MaxWiK' ){
        par.est  =  get_Spider_MAP( par.sim = par.sim, stat.sim = stat.sim, 
                                    stat.obs = stat.obs, n_best = 1 )
    }
    
    ### Get MSE 
    MSE = NULL
    if ( !is.na( par.est )[1] ) MSE  =  sum( ( par.truth - par.est ) ** 2  )
    
    running_time  =  as.numeric( difftime(Sys.time(), time_start, units = "secs")[[1]] )
    
    return( list( stat = data.frame(    method_name = method_name,
                                        kernel_name = kernel_name,
                                        MSE = MSE, 
                                        running_time = running_time), 
                  par.est = par.est ) )
}

#' @describeIn Get_call Function to generate parameters and simulate a model based on MaxWiK algorithm 
#'
#' @param model is a function to get output of simulation during sampling 
#' @param arg0 is a list with arguments for a model function, so that arg0 is NOT changed during sampling
#' @param size Number of point to restrict original dataset
#' @param nmax is maximal number of iterations
#'
#' @return \code{sampler_method()} returns the list: \cr
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
sampler_method  <-  function( stat.obs, stat.sim, par.sim, model, 
                              method_name, kernel_name, 
                              arg0 = list(),  size = 500, 
                              nmax = 30, 
                              model_name, dimension, stochastic_term, par.truth ){
    
    stat.sim_itt =  stat.sim
    par.sim_itt  =  par.sim
    
    res = NULL
    for( itt in 1:nmax ){

        new_par = Get_parameter( method_name = method_name, 
                            kernel_name = kernel_name, 
                            stat.obs   =  stat.obs, 
                            stat.sim   =  stat.sim_itt, 
                            par.sim    =  par.sim_itt, 
                            G = G, 
                            par.truth = par.truth )
        new_par$stat$iteration = itt
        res = rbind( res, new_par$stat )
        
        par.est  =  data.frame( matrix( new_par$par.est, ncol = dimension ) )
        new_sim  =  do.call( what = model, args = c( arg0, list( parameter =  par.est ) ) )
        names( par.est )  =  names( par.sim_itt )
        names( new_sim )  =  names( stat.sim_itt )
        par.sim_itt  =  rbind( par.sim_itt, par.est )
        stat.sim_itt =  rbind( stat.sim_itt, new_sim )
    }
    
    return( list( results = res, 
                  par.sim_itt = par.sim_itt,
                  stat.sim_itt = stat.sim_itt ) )
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
    
    check_packages()
     
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
            
            web = spiderweb_slow( param = par_sim, stat.sim = stat_sim, stat.obs = stat_obs, 
                             psi = psi_t$psi[ j ], t = psi_t$t[ j ], talkative = FALSE )
            
            res_1  =  web$par.best
            
            new_best.sim  =  do.call( what = model, args = c( arg0, list( parameter =  web$par.best ) ) )

            res_1[ , (ncol(res_1) + 1) : (ncol(res_1) + ncol(new_best.sim) ) ]  =  new_best.sim
            
            mse_1   =  MSE_sim( stat.obs = stat.obs, stat.sim = new_best.sim )
            mse_1   =  as.data.frame( mse_1 )
            names( mse_1)  = 'mse'
            
            res_1$mse    =  as.numeric( mse_1$mse )
            res_1$comm     =     'Best'


            res_1$sim      =     web$sim.best

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
    
    return( list( results = results_ALL, 
                  best = best, 
                  MSE_min = err, 
                  number_of_iterations = i, 
                  time = as.numeric( difftime( end_time, start_time, units = "secs")[[1]] )
                ) )
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





