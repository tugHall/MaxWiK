# This is a sampler for ABC method that based on maxima weighted isolation kernel method

#' @describeIn Function to generate parameters and simulate a model based on MaxWiK algorithm 
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
                              nmax = 30, G,
                              model_name, dimension, stochastic_term, par.truth ){
    
    stat.sim_itt =  stat.sim
    par.sim_itt  =  par.sim
    
    res = NULL
    for( itt in 1:nmax ){
        
        input  =  restrict_data( par.sim = par.sim_itt, 
                                  stat.sim = stat.sim_itt, 
                                  stat.obs = stat.obs, 
                                  size = size ) 
        stat.sim_itt =  input$stat.sim
        par.sim_itt  =  input$par.sim

        new_par = Get_parameter( method_name = method_name, 
                            kernel_name = kernel_name, 
                            stat.obs   =  stat.obs, 
                            stat.sim   =  stat.sim_itt, 
                            par.sim    =  par.sim_itt, 
                            G = G, 
                            par.truth = par.truth )
        new_par$stat$iteration = itt
        res = rbind( res, new_par$stat )
        
        par.est  =  new_par$par.est
        names( par.est )  =  names( par.sim_itt )
        if ( is.data.frame(  new_par$par.est ) ){
            add_arg  =  list( parameter =  as.data.frame( new_par$par.est ) )
        } else {
            add_arg  =  list( parameter =  data.frame( matrix( new_par$par.est, nrow = 1 ) ) )
        }
        new_sim  =  do.call( what = model, args = c( arg0, add_arg ) )

        names( new_sim )  =  names( stat.sim_itt )
        par.sim_itt  =  rbind( par.sim_itt, par.est )
        stat.sim_itt =  rbind( stat.sim_itt, new_sim )
    }
    
    return( list( results = res, 
                  par.sim_itt = par.sim_itt,
                  stat.sim_itt = stat.sim_itt ) )
}


#' @describeIn Function to call all the methods to get estimation of parameter and MSE
#'
#' @return \code{sampler_all_methods()} returns data.frame with MSE for all defined methods
#' 
#' @param cores Number of cores for parallel calculation with iterations
#' 
#' @export
#'
#' @examples
#' NULL 
sampler_all_methods  <-  function( model_name, dimension, stochastic_term, 
                                    stat.obs, stat.sim, par.sim, G, par.truth,
                                    cores = 4, nmax ){
    
    DF  =  NULL
    Meth_Kern  =  data.frame( Method = c('K2-ABC', 'K2-ABC', 'K2-ABC', 'Rejection', 
                                         'Loclinear', 'Neuralnet', 'Ridge',
                                         'MaxWiK_MAP', 'MaxWiK' ), 
                              Kernel = c('Gaussian', 'Laplacian', 'iKernel', '',
                                         '',          '',          'epanechnikov',
                                         'iKernel', 'iKernel') 
    )
    
    # Define model function:
    if ( model_name == 'Gaussian' ) {
        model  =  model
        arg0 = list(  name = c( 'Gaussian', 'Linear' )[1],
                      x0 = par.truth, 
                      stat.obs = stat.obs, 
                      noise = stochastic_term )
    } else {
        model  =  model
        arg0 = list(  name = c( 'Gaussian', 'Linear' )[2],
                      x0 = par.truth, 
                      stat.obs = stat.obs, 
                      noise = stochastic_term )
    } 

    DF_1_S  =  mclapply( 1:nrow(Meth_Kern) , FUN = function( mk ){   
                        sampler_method( method_name =  as.character( Meth_Kern$Method[mk] ), 
                                        kernel_name =  as.character( Meth_Kern$Kernel[mk] ), 
                                        model_name  =  model_name, 
                                        dimension   =  dimension, 
                                        stochastic_term  =  stochastic_term, 
                                        nmax = nmax,
                                        stat.obs   =  stat.obs, 
                                        stat.sim   =  stat.sim, 
                                        par.sim    =  par.sim, 
                                        G          =  G, 
                                        par.truth  =  par.truth, 
                                        model      =  model, 
                                        arg0       =  arg0,  
                                        size       =  500  
                    )
                }, mc.cores  =  cores )

        # Check an error
        bad   =  sapply( DF_1_S, inherits, what = "try-error" )
        # If NO errors:
        if ( any( !bad ) ){
            for( mk in ( 1:nrow(Meth_Kern) )[ !bad ] ){
                DF_1  =  DF_1_S[[ mk ]]$results
                DF    =  rbind( DF, DF_1 )
            }
            # DF_2  =  do.call( rbind, DF_1_S[ !bad ]$results )
            # DF    =  rbind( DF, DF_2 )
        }
        # If error in some core(s):
        if ( any( bad ) ){
            its = which( bad )
            
            DF_3  =   data.frame(   method_name = as.character( Meth_Kern$Method[ its ] ),
                                    kernel_name = as.character( Meth_Kern$Kernel[ its ] ),
                                    MSE = NA, 
                                    running_time = NA ) 
            # Add circle for its:
            for( i in its ){
                DF_3[ 1, 'iteration' ]  =  i
                DF    =  rbind( DF, DF_3 )
            }
        }

    return( DF )
}


#' Function to generate parameters and simulate a model based on MaxWiK algorithm 
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
#' 
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
    smpl_1  =  sampler_MaxWiK( stat.obs, stat.sim, par.sim,  model = model, 
                                         arg0 = list(  name = c( 'Gauss', 'Linear' )[1],
                                                       x0 = x0, 
                                                       stat.obs = stat.obs, 
                                                       noise = 0 ), 
                                         size = 200, psi_t, epsilon = 1E-12, 
                                         nmax = 30, include_top = TRUE,
                                         slowly = TRUE, rate = rat[ r ] )
    if ( r == 1 ){ smpl = smpl_1 } else { smpl = Map( c, smpl, smpl_1 ) }
}





