# This is a sampler for ABC method that based on maxima weighted isolation kernel method

#' Function to generate parameters and simulate a model based on MaxWiK algorithm 
#'
#' @param stat.sim Summary statistics of the simulations (model output)
#' @param stat.obs Summary statistics of the observation point
#' @param par.sim Data frame of parameters of the model
#' @param psi_t Vector of psi and t hyperparameters.
#' @param model Function to get output of simulation during sampling 
#' @param arg0 List with arguments for a model function, so that arg0 is NOT changed during sampling
#' @param size Number of points in the simulation based on MaxWiK algorithm 
#' @param nmax Maximal number of iterations
#' @param epsilon Criterion to stop simulation when \code{MSE_current - MSE_previous < epsilon}
#' @param check_err Logical parameter to check epsilon or do not
#' @param include_top Logical to include top points (network) from \code{spider_web()} function to simulate or do not
#' @param slowly Logical for two algorithms: slow and fast seekers in sampling
#' @param rate Rate value in the range \code{[0,1]} to define 
#' the rate of changing in the original top of sampled points for slow scheme (if slowly = TRUE)
#' @param n_simulation_stop Maximal number of simulations to stop sampling. 
#' If \code{n_simulation_stop = NA} then there is no restriction (by default)
#' @param include_web_rings Logical to include or do not include the cobweb rings to the simulations
#' @param number_of_nodes_in_ring Number of points/nodes between two points in the web ring. By default \code{number_of_nodes_in_ring = 2}
#'
#' @returns \code{sampler_MaxWiK()} returns the list: \cr
#' - results: results of all the simulations; \cr 
#' - best: the best value of parameter; \cr
#' - MSE_min: minimum of MSE; \cr
#' - number_of_iterations: number of iterations; \cr
#' - time: time of sampling in seconds, \cr
#' - n_simulations: the total number of simulations.
#' 
#' @export
#'
#' @examples
#' MaxWiK::MaxWiK_templates(dir = tempdir()) # See the template 'MaxWiK.Sampling.R' 
#' # and vignettes for usage.
sampler_MaxWiK  <-  function( stat.obs, stat.sim, par.sim, model, 
                                             arg0 = list(),  size = 500, 
                                             psi_t, epsilon, nmax = 100, 
                                             include_top = FALSE,
                                             slowly = FALSE, rate = 0.2, 
                                             n_simulation_stop = NA, 
                                             check_err  =  TRUE, 
                                             include_web_rings  =  TRUE,
                                             number_of_nodes_in_ring = 2 ){ 
    # epsilon is a criterion to stop simulation
    # nmax is maximal number of iterations
    # psi_t is a data.frame of psi and t with top values of MSE of length 20
    # size is a size of data frame for a single simulation
    # i  =  sample( 1:nrow( psi_t ), 1, replace = TRUE, prob = psi_t$prob )
    
    # model is a function to get output of simulation during sampling 
    # arg0 is a list with arguments for a model function, so that arg0 is NOT changed during sampling
    
    rep_check  =  FALSE  # Logical check to repeat if err is same
    
    check_packages()
     
    start_time = Sys.time()
    # Redefine input:
    stat_sim  =  stat.sim     #  initial data for summary statistics of simulations
    par_sim   =  par.sim      #  initial data for parameters of simulations
    stat_obs  =  stat.obs     #  
    
    ### Restrict data size in accordance with size parameter:
    rstrct  =  restrict_data( par.sim = par_sim, stat.sim = stat_sim, stat.obs = stat_obs, size = size )
    stat_sim  =  rstrct$stat.sim     #  Restricted initial data for summary statistics of simulations
    par_sim   =  rstrct$par.sim      #  Restricted initial data for parameters of simulations
    
    
    # Get dimensions of parameter and simulation datasets:
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
    
    ############### Function to get data.frame of results:
    make_res  =  function( web, use_network  =  FALSE, model, arg0, n_simulations ){
        
        results  =  NULL
        
        if ( !use_network ) {
            res_1  =  web$par.best
            
            new_best.sim  =  do.call( what = model, args = c( arg0, list( x =  web$par.best ) ) )
            n_simulations  <-  n_simulations  +  1
            
            res_1[ , (ncol(res_1) + 1) : (ncol(res_1) + ncol(new_best.sim) ) ]  =  new_best.sim
            
            mse_1   =  MSE_sim( stat.obs = stat.obs, stat.sim = new_best.sim )
            mse_1   =  as.data.frame( mse_1 )
            names( mse_1)  = 'mse'
            
            res_1$mse    =  as.numeric( mse_1$mse )
            res_1$comm     =     'Best'
            
            res_1$sim      =     web$sim.best
            res_1$sim_ID   =     n_simulations
            results  =  res_1
        } else {
            
            for( i in 1:nrow( web$network ) ){
                
                res_1  =  web$network[ i, ]
                
                new.sim  =  do.call( what = model, args = c( arg0, list( x =  res_1 ) ) )
                n_simulations  <-  n_simulations  +  1
                
                res_1[ , (ncol(res_1) + 1) : (ncol(res_1) + ncol(new.sim) ) ]  =  new.sim
                
                mse_1   =  MSE_sim( stat.obs = stat.obs, stat.sim = new.sim )
                mse_1   =  as.data.frame( mse_1 )
                names( mse_1)  = 'mse'
                
                res_1$mse    =  as.numeric( mse_1$mse )
                res_1$comm     =     'Network'
                
                res_1$sim      =     web$sim_network[ i ]
                res_1$sim_ID   =     n_simulations
                
                results  =  rbind( results, res_1 )
            }
        }
        
        return( list( results       =  results, 
                      n_simulations =  n_simulations ) )
    }
    
    make_web_rings  =  function( web, model, arg0, number_of_nodes  =  1, n_simulations ){
        
        results  =  NULL
        
        res_net  =  web$network
        dst    =  as.matrix( dist( x = res_net ) )
        
        if ( nrow( web$network ) < 2 ) return( NULL )
        for( i in 1:( nrow( web$network ) - 2 ) ){
            
            pnt_1  =  web$network[ i, ]
            dst_1  =  dst[ -1,    ]
            dst    =  dst[ -1, -1 ]
            res_net  =  res_net[ -1, ]
            pnt_2  =  res_net[ as.integer( which.min( dst_1[ , 1 ] )[ 1 ] ), ]

            pnts  =  generate_points_between_two_points( pair = rbind( pnt_1, pnt_2 ), n = number_of_nodes + 2 )
            pnts  =  pnts[ -1, ]
            pnts  =  pnts[ -nrow( pnts ), ]
            
            for( j in 1:nrow( pnts ) ){
                
                res_1    =   pnts[ j, ]
            
                new.sim  =  do.call( what = model, args = c( arg0, list( x =  res_1 ) ) )
                n_simulations  <-  n_simulations  +  1
                
                res_1[ , (ncol(res_1) + 1) : (ncol(res_1) + ncol(new.sim) ) ]  =  new.sim
                
                mse_1   =  MSE_sim( stat.obs = stat.obs, stat.sim = new.sim )
                mse_1   =  as.data.frame( mse_1 )
                names( mse_1)  = 'mse'
                
                res_1$mse    =  as.numeric( mse_1$mse )
                res_1$comm     =     'Cobweb_Ring'
                
                res_1$sim      =     web$sim_network[ i ]
                res_1$sim_ID   =     n_simulations
                
                results  =  rbind( results, res_1 )
            }
        }
        
        return( list( results       =  results, 
                      n_simulations =  n_simulations ) )
    }
    
    n_simulations  =  0
    for( i in 1:nmax ){
        
        # Progress BAR
        setTxtProgressBar(pb, i)
        
        results  =  data.frame( NULL )
        
        
        ### Get new parameters and corresponding results of simulations based on each set of psi and t
        # source( './lib_iKernel.R' )
        for ( j in 1:nrow( psi_t ) ){
            
            web = meta_sampling( param = par_sim, stat.sim = stat_sim, stat.obs = stat_obs, 
                             psi = psi_t$psi[ j ], t = psi_t$t[ j ], talkative = FALSE )
            
            res_1  =  make_res( web = web, 
                                use_network  = include_top, 
                                model =  model, 
                                arg0  =  arg0, 
                                n_simulations = n_simulations ) 
            
            n_simulations = res_1$n_simulations
            res_1  =  res_1$results
            
            if ( include_top & include_web_rings ) {
                res_2  =  make_web_rings( web    =  web, 
                                          model  =  model, 
                                          arg0   =  arg0, 
                                          number_of_nodes = number_of_nodes_in_ring, 
                                          n_simulations = n_simulations )
                n_simulations = res_2$n_simulations
                res_2  =  res_2$results
                
            } else { 
                res_2  =  NULL
            }
            
            if ( FALSE ){
                        res_1  =  web$par.best
                        
                        new_best.sim  =  do.call( what = model, args = c( arg0, list( x =  web$par.best ) ) )
                        n_simulations  =  n_simulations  +  1
            
                        res_1[ , (ncol(res_1) + 1) : (ncol(res_1) + ncol(new_best.sim) ) ]  =  new_best.sim
                        
                        mse_1   =  MSE_sim( stat.obs = stat.obs, stat.sim = new_best.sim )
                        mse_1   =  as.data.frame( mse_1 )
                        names( mse_1)  = 'mse'
                        
                        res_1$mse    =  as.numeric( mse_1$mse )
                        res_1$comm     =     'Best'
                        res_1$sim_ID   =     n_simulations
            
            
                        res_1$sim      =     web$sim.best
            }

            results  =  rbind( results, res_1, res_2 )
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
        if ( check_err & ( abs(err - err_previous) < epsilon ) ) {
            if ( rep_check ){
                break 
            } else { rep_check  =  TRUE }
        } else rep_check  =  FALSE
            
        err_previous  =  err
        if ( !is.na( n_simulation_stop ) & 
             ( n_simulation_stop <= n_simulations ) ) break
        
    }
    
    best  =  results_ALL[ which.min(results_ALL$mse ) , ]
    end_time = Sys.time()
    
    #Stop progress BAR:
    close(pb)
    
    return( list( results = results_ALL, 
                  best = best, 
                  MSE_min = err, 
                  number_of_iterations = i, 
                  time = as.numeric( difftime( end_time, start_time, units = "secs")[[1]] ),
                  n_simulations  =  n_simulations
                ) )
}


# This is a parallel implementation of sampler sampler_MaxWiK

#' @describeIn sampler_MaxWiK Function to generate parameters and simulate a model based on MaxWiK algorithm 
#'
#' @param cores Number of cores for parallel calculations of a model (4 by default)
#'
#' @returns the same as in \code{sampler_MaxWiK()}.
#' 
#' @export
#'
#' @examples
#' MaxWiK::MaxWiK_templates(dir = tempdir()) # See the template 'MaxWiK.Sampling.R' 
#' # and vignettes for usage. For parallel implementation 
#' # change the function 'sampler_MaxWiK()' to 'sampler_MaxWiK_parallel()'.
sampler_MaxWiK_parallel  <-  function(    stat.obs, stat.sim, par.sim, model, 
                                          arg0 = list(),  size = 500, 
                                          psi_t, epsilon, nmax = 100, 
                                          include_top = FALSE,
                                          slowly = FALSE, rate = 0.2, 
                                          n_simulation_stop = NA, 
                                          check_err  =  TRUE, 
                                          include_web_rings  =  TRUE,
                                          number_of_nodes_in_ring = 2, 
                                          cores = 4 ){ 
    # epsilon is a criterion to stop simulation
    # nmax is maximal number of iterations
    # psi_t is a data.frame of psi and t with top values of MSE of length 20
    # size is a size of data frame for a single simulation
    # i  =  sample( 1:nrow( psi_t ), 1, replace = TRUE, prob = psi_t$prob )
    
    # model is a function to get output of simulation during sampling 
    # arg0 is a list with arguments for a model function, so that arg0 is NOT changed during sampling
    
    rep_check  =  FALSE  # Logical check to repeat if err is same
    
    check_packages()
    
    n_simulations  <-  0
    
    start_time = Sys.time()
    # Redefine input:
    stat_sim  =  stat.sim     #  initial data for summary statistics of simulations
    par_sim   =  par.sim      #  initial data for parameters of simulations
    stat_obs  =  stat.obs     #  
    
    ### Restrict data size in accordance with size parameter:
    rstrct  =  restrict_data( par.sim = par_sim, stat.sim = stat_sim, stat.obs = stat_obs, size = size )
    stat_sim  =  rstrct$stat.sim     #  Restricted initial data for summary statistics of simulations
    par_sim   =  rstrct$par.sim      #  Restricted initial data for parameters of simulations
    
    
    # Get dimensions of parameter and simulation datasets:
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
    
    ############### Function to get data.frame of results:
    make_res  =  function( web, use_network  =  FALSE, model, arg0, cores, n_simulations ){
        
        results  =  NULL
        
        if ( !use_network ) {
            res_1  =  web$par.best
            
            new_best.sim  =  do.call( what = model, args = c( arg0, list( x =  web$par.best ) ) )
            n_simulations  <-  n_simulations  +  1
            
            res_1[ , (ncol(res_1) + 1) : (ncol(res_1) + ncol(new_best.sim) ) ]  =  new_best.sim
            
            mse_1   =  MSE_sim( stat.obs = stat.obs, stat.sim = new_best.sim )
            mse_1   =  as.data.frame( mse_1 )
            names( mse_1)  = 'mse'
            
            res_1$mse    =  as.numeric( mse_1$mse )
            res_1$comm     =     'Best'
            
            res_1$sim      =     web$sim.best
            results  =  res_1
        } else {
            
            ## Parallel implementation:
            ################
            # for( i in 1:nrow( web$network ) )
              fun_par =  function( i, web, model, arg0, stat.obs ){
                
                    res_1  =  web$network[ i, ]
                    
                    new.sim  =  do.call( what = model, args = c( arg0, list( x =  res_1 ) ) )
                    # update_n( 1 )
                    
                    res_1[ , (ncol(res_1) + 1) : (ncol(res_1) + ncol(new.sim) ) ]  =  new.sim
                    
                    mse_1   =  MSE_sim( stat.obs = stat.obs, stat.sim = new.sim )
                    mse_1   =  as.data.frame( mse_1 )
                    names( mse_1)  = 'mse'
                    
                    res_1$mse    =  as.numeric( mse_1$mse )
                    res_1$comm     =     'Network'
                    
                    res_1$sim      =     web$sim_network[ i ]
                    
                    return( res_1 )
                    #  results  =  rbind( results, res_1 )
            }
              fun_par_run =  function( x ) fun_par( i = x, web = web, 
                                                    model = model, 
                                                    arg0 = arg0, stat.obs = stat.obs )
              
              results_list   =  mclapply( X = 1:nrow( web$network ), 
                                          FUN = function( x ){
                                              update_n( x )
                                              fun_par_run( x )
                                              }, 
                                          mc.cores = cores )
              n_simulations  <-  n_simulations  +  nrow( web$network )
              results  =  as.data.frame( do.call( rbind, results_list ) )
            ##############
        }
        
        return( list( results      =  results, 
                      n_simulations =  n_simulations ) )
    }
    
    make_web_rings  =  function( web, model, arg0, number_of_nodes  =  1, cores, n_simulations ){
        
        results  =  NULL
        
        res_net  =  web$network
        dst    =  as.matrix( dist( x = res_net ) )
        
        model_paral  =  function( j, pnts, model, arg0, stat.obs ){
            res_1    =   pnts[ j, ]
            
            new.sim  =  do.call( what = model, args = c( arg0, list( x =  res_1 ) ) )
            # update_n( 1 )
            
            res_1[ , (ncol(res_1) + 1) : (ncol(res_1) + ncol(new.sim) ) ]  =  new.sim
            
            mse_1   =  MSE_sim( stat.obs = stat.obs, stat.sim = new.sim )
            mse_1   =  as.data.frame( mse_1 )
            names( mse_1)  = 'mse'
            
            res_1$mse    =  as.numeric( mse_1$mse )
            res_1$comm     =     'Cobweb_Ring'
            
            res_1$sim      =     web$sim_network[ i ]
            return( res_1 )
        }
        model_paral_run  =  function( x ) model_paral( j = x, pnts = pnts, 
                                               model = model, arg0 = arg0, 
                                               stat.obs = stat.obs ) 
        
        if ( nrow( web$network ) < 2 ) return( NULL )
        for( i in 1:( nrow( web$network ) - 2 ) ){
            
            pnt_1  =  web$network[ i, ]
            dst_1  =  dst[ -1,    ]
            dst    =  dst[ -1, -1 ]
            res_net  =  res_net[ -1, ]
            pnt_2  =  res_net[ as.integer( which.min( dst_1[ , 1 ] )[ 1 ] ), ]
            
            pnts  =  generate_points_between_two_points( pair = rbind( pnt_1, pnt_2 ), n = number_of_nodes + 2 )
            pnts  =  pnts[ -1, ]
            pnts  =  pnts[ -nrow( pnts ), ]
            
            ### Parallel implementation:
            results_list  =  mclapply( 
                X = 1:nrow( pnts ), 
                FUN =  function( x ){
                        update_n( x )
                        model_paral_run( x )
                        }, 
                mc.cores = cores 
                )
            n_simulations  <-  n_simulations  +  nrow( pnts )

            results  =  as.data.frame( do.call( rbind, results_list ) )
            ### 
        }
        return( list( results      =  results, 
                      n_simulations =  n_simulations ) )
    }
    
    for( i in 1:nmax ){
        
        # Progress BAR
        setTxtProgressBar(pb, i)
        
        results  =  data.frame( NULL )
        
        
        ### Get new parameters and corresponding results of simulations based on each set of psi and t
        # source( './lib_iKernel.R' )
        for ( j in 1:nrow( psi_t ) ){
            
            web = meta_sampling( param = par_sim, stat.sim = stat_sim, stat.obs = stat_obs, 
                                  psi = psi_t$psi[ j ], t = psi_t$t[ j ], talkative = FALSE )
            
            res_1  =  make_res( web = web, 
                                use_network  = include_top, 
                                model =  model, 
                                arg0  =  arg0, 
                                cores = cores, 
                                n_simulations = n_simulations )
            n_simulations  =  res_1$n_simulations
            res_1  =  res_1$results
            
            if ( include_top & include_web_rings ) {
                res_2  =  make_web_rings( web    =  web, 
                                          model  =  model, 
                                          arg0   =  arg0, 
                                          number_of_nodes = number_of_nodes_in_ring, 
                                          cores = cores, 
                                          n_simulations = n_simulations )
                
                n_simulations  =  res_2$n_simulations
                res_2  =  res_2$results
                
            } else res_2  =  NULL
            
            if ( FALSE ){
                res_1  =  web$par.best
                
                new_best.sim  =  do.call( what = model, args = c( arg0, list( x =  web$par.best ) ) )
                n_simulations  =  n_simulations  +  1
                
                res_1[ , (ncol(res_1) + 1) : (ncol(res_1) + ncol(new_best.sim) ) ]  =  new_best.sim
                
                mse_1   =  MSE_sim( stat.obs = stat.obs, stat.sim = new_best.sim )
                mse_1   =  as.data.frame( mse_1 )
                names( mse_1)  = 'mse'
                
                res_1$mse    =  as.numeric( mse_1$mse )
                res_1$comm     =     'Best'
                
                
                res_1$sim      =     web$sim.best
            }
            
            results  =  rbind( results, res_1, res_2 )
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
        if ( check_err & ( abs(err - err_previous) < epsilon ) ) {
            if ( rep_check ){
                break 
            } else { rep_check  =  TRUE }
        } else rep_check  =  FALSE
        
        err_previous  =  err
        if ( !is.na( n_simulation_stop ) & 
             ( n_simulation_stop <= n_simulations ) ) break
        
    }
    
    best  =  results_ALL[ which.min(results_ALL$mse ) , ]
    end_time = Sys.time()
    
    #Stop progress BAR:
    close(pb)
    
    return( list( results = results_ALL, 
                  best = best, 
                  MSE_min = err, 
                  number_of_iterations = i, 
                  time = as.numeric( difftime( end_time, start_time, units = "secs")[[1]] ),
                  n_simulations  =  n_simulations
    ) )
}



