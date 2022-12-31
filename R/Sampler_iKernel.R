# This is a sampler for ABC method that based on maxima weighted isolation kernel method

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
        
        par.est  =  Get_MAP( as.data.frame( lc ) )
        
        #if (nrow( unique.data.frame(lc) ) > 1 ){
        #    par.est  =  point_estimate( lc )$MAP
        #} else {
        #    par.est  =  unique.data.frame(lc)
        #}
    }
    
    if ( method_name == 'Neuralnet' ){
        
        res      =  adjust_ABC_tolerance( par.sim = par.sim, stat.sim = stat.sim, 
                                          stat.obs = stat.obs )
        tol   =   res$tolerance
        if ( tol * nrow(stat.sim) < n_min ) tol  =  n_min / nrow( stat.sim )
        
        nn   =  abc(   target = stat.obs, param = par.sim, sumstat = stat.sim, 
                       tol = tol, method  =  'neuralnet', hcorr   =  TRUE, 
                       transf=c("none","log"), lambda = 0.0001, trace = FALSE )
        
        # nn_adj = round( nn$adj.values + runif( n= length( nn$adj.values ) ) / 1E6, digits = 9 )
        # par.est  =  point_estimate( nn_adj )$MAP
        par.est  =  Get_MAP( as.data.frame( nn$adj.values ) )
    }
    
    if ( method_name == 'Ridge' ){
        
        res      =  adjust_ABC_tolerance( par.sim = par.sim, stat.sim = stat.sim, 
                                          stat.obs = stat.obs )
        tol   =   res$tolerance
        if ( tol * nrow(stat.sim) < n_min ) tol = n_min / nrow( stat.sim )
        
        rdg   =  abc(   target = stat.obs, param = par.sim, sumstat = stat.sim, 
                        tol = tol, method  =  'ridge', hcorr   =  FALSE, 
                        transf=c("none","log"), kernel = 'epanechnikov' )
        
        # rdg_adj = round( rdg$adj.values + runif( n= length( rdg$adj.values ) )/1E6, digits = 9 )
        
        # if ( nrow( unique.data.frame( rdg_adj ) ) > 1 ){
        #     par.est  =  point_estimate( rdg_adj )$MAP
        # } else {
        #     par.est  =  unique.data.frame(rdg_adj)
        # }
        
        par.est  =  Get_MAP( as.data.frame( rdg$adj.values ) )
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
    
    running_time  =  as.numeric( difftime(Sys.time(), time_start, units = "secs")[[ 1 ]] )
    
    return( list( stat = data.frame(    method_name = method_name,
                                        kernel_name = kernel_name,
                                        MSE = MSE, 
                                        running_time = running_time), 
                  par.est = par.est ) )
}

#' @describeIn Get_parameter Function to generate parameters and simulate a model based on MaxWiK algorithm 
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


#' @describeIn Get_parameter  Function to call all the methods to get estimation of parameter and MSE
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


#' @describeIn Get_parameter Make experiments with sampling for all the methods and one case of toy model and dimension
#'
#' @param file_name File to save results
#' @param rng Range of points for each dimension, by default \code{rng = range(0,10)}
#' @param restrict_points_number Number of points which will be used in each method and that are close to observation point
#'
#' @return \code{experiment_samplers()} returns data frame with results of experiment for sampling using all the methods
#' 
#' @export
#'
#' @examples
#' NULL
experiment_samplers  <-  function( file_name = './output.txt', 
                                   model_name = 'Gaussian',
                                   dimension = 6, 
                                   stochastic_term  =  5,
                                   rng  =  c( 0,10 ), 
                                   restrict_points_number = 500, 
                                   nmax = 30 ){
    ### Check installation of libraries:
    check_packages()
    
    # delete old file
    if ( file.exists( file_name) ) unlink( file_name )
    
    input  =  NULL
    x0  =  round( runif( n = dimension, min = rng[1], max = rng[2] ), digits = 4 )
    Number_of_points  =  max( c( 50 * dimension, restrict_points_number ) )
    
    A  =  ( ( 1:dimension ) + 5 ) * 100
    sigma  =  rep( ( rng[2] - rng[1] ) / 10, dimension )
    if ( model_name == 'Gaussian' ) {
        
        input = get_dataset_of_Gaussian_model( d = dimension, x0 = x0, probability = TRUE, 
                                n = Number_of_points, r = rng, A = A, sigma = sigma, 
                                noise = stochastic_term )
    }
    if ( model_name == 'Linear' ) {
        input  =  get_dataset_of_Linear_model( d = dimension, x0 = x0, probability = TRUE, 
                                n = Number_of_points, r = rng, A = A, 
                                noise = stochastic_term )
    }
    
    if ( is.null( input ) ) stop( 'Model name is incorrect' )
    stat.sim_origin  =  input$stat.sim
    stat.obs  =  input$stat.obs
    par.sim_origin  =  input$par.sim
    rm( input )
    
    # Apply restrict number of points:
    tol = restrict_points_number / nrow( stat.sim_origin )
    rej = abc::abc( target = stat.obs, param = par.sim_origin, sumstat = stat.sim_origin,
                    method = 'rejection', tol = tol )
    
    stat.sim  =  stat.sim_origin[ rej$region, ]
    par.sim   =  par.sim_origin[ rej$region, ] 
    par.truth =  data.frame( matrix( x0, ncol = dimension ) )
    rm(x0)
    psi_t  =  adjust_psi_t( par.sim = par.sim, stat.sim = stat.sim, 
                            stat.obs = stat.obs, talkative = FALSE, 
                            check_pos_def = FALSE, 
                            n_best = 10, cores = 4 )
    
    ikern  =  iKernelABC( psi = psi_t$psi[1], t = psi_t$t[1], 
                          param = par.sim, 
                          stat.sim = stat.sim, 
                          stat.obs = stat.obs, 
                          talkative = FALSE, 
                          check_pos_def = FALSE )
    
    G = matrix( data = ikern$similarity, ncol = 1 )
    
    RES = sampler_all_methods(  model_name = model_name, 
                                dimension = dimension, 
                                stochastic_term = stochastic_term, 
                                stat.obs = stat.obs, 
                                stat.sim = stat.sim, 
                                par.sim = par.sim, 
                                par.truth = par.truth, 
                                G = G,
                                nmax = nmax
    )
    # Sampler simulation:
    if( FALSE ){
        smpl_MaxWiK  =  sampler_MaxWiK( stat.obs = stat.obs, stat.sim = stat.sim, 
                                        par.sim  = par.sim,  model = model, 
                                        arg0 = arg0, 
                                        size = 500, psi_t, epsilon = 1E-12, 
                                        nmax = 30, include_top = FALSE,
                                        slowly = TRUE, rate = 0.1 )
    }
    
    
    return( RES )
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





