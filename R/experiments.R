
### File to organize toy experiments

#' Function to call a method to get MSE
#'
#' @param method_name Name of a method
#' @param kernel_name Name of kernel function
#' @param model_name Name of a model
#' @param dimension  Dimension of the model
#' @param stochastic_term A number (usually in the range \code{[0,1]}) of stochastic term in the dataset \code{stat.sim}
#' @param iteration Iteration number of trial with the same model and other parameters
#' @param stat.obs Data frame of statistics of observation point
#' @param stat.sim Data frame of statistics of simulations
#' @param par.sim Data frame of parameters
#' @param G Matrix of similarities for K2-ABC based on isolation kernel
#' @param par.truth Truth parameter value to check result of estimation
#'
#' @return \code{Get_call()} returns the list: \cr
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
Get_call  <-  function( method_name, kernel_name = '', model_name, dimension, stochastic_term, iteration, 
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
    
    return( data.frame(   method_name = method_name,
                          kernel_name = kernel_name,
                          model_name  = model_name,
                          dimension   = dimension, 
                          stochastic_term = stochastic_term,
                          MSE = MSE, 
                          running_time = running_time, 
                          iteration = iteration ) )
}


#' @describeIn Get_call  Function to call all the methods to get estimation of parameter and MSE
#'
#' @return \code{Get_call_all_methods()} returns data.frame with MSE for all defined methods
#' 
#' @param cores Number of cores for parallel calculation with iterations
#' 
#' @export
#'
#' @examples
#' NULL 
Get_call_all_methods  <-  function( model_name, dimension, stochastic_term, iterations, 
                                    stat.obs, stat.sim, par.sim, G, par.truth,
                                    cores = 4 ){
    
    DF  =  NULL
    Meth_Kern  =  data.frame( Method = c('K2-ABC', 'K2-ABC', 'K2-ABC', 'Rejection', 
                                         'Loclinear', 'Neuralnet', 'Ridge',
                                         'MaxWiK_MAP', 'MaxWiK' ), 
                              Kernel = c('Gaussian', 'Laplacian', 'iKernel', '',
                                         '',          '',          'epanechnikov',
                                         'iKernel', 'iKernel') 
                              )
    
    for( mk in 1:nrow(Meth_Kern) ){
         DF_1  =  mclapply( iterations , FUN = function( x ){     
                        Get_call(  method_name =  as.character( Meth_Kern$Method[mk] ), 
                                   kernel_name =  as.character( Meth_Kern$Kernel[mk] ), 
                                   model_name  =  model_name, 
                                   dimension   =  dimension, 
                                   stochastic_term  =  stochastic_term, 
                                   iteration  =  x,
                                   stat.obs   =  stat.obs, 
                                   stat.sim   =  stat.sim, 
                                   par.sim    =  par.sim, 
                                   G          =  G, 
                                   par.truth  =  par.truth 
                                )
        }, mc.cores  =  cores )
        # Check an error
        bad   =  sapply( DF_1, inherits, what = "try-error" )
        # If NO errors:
        if ( any( !bad ) ){
            DF_2  =  do.call( rbind, DF_1[ !bad ] )
            DF    =  rbind( DF, DF_2 )
        }
        # If error in some core(s):
        if ( any( bad ) ){
            its = which( bad )
            DF_3  =   data.frame(   method_name = as.character( Meth_Kern$Method[mk] ),
                                    kernel_name = as.character( Meth_Kern$Kernel[mk] ),
                                    model_name  = model_name,
                                    dimension   = dimension, 
                                    stochastic_term = stochastic_term,
                                    MSE = NA, 
                                    running_time = NA, 
                                    iteration = NA ) 
            # Add circle for its:
            for( i in its ){
                DF_3[ 1, 'iteration' ]  =  i
                DF    =  rbind( DF, DF_3 )
            }
        }
    }
    
    return( DF )
}

#' Function to prepare toy experiments
#' 
#' @param file_name Name of file to output results
#' @param models Names of models for simulation, by default
#' \code{models = c( 'Gaussian', 'Linear' )}
#' @param dimensions Dimensions of models, by default
#' \code{ dimensions = (1:20)*2 }
#' @param stochastic_terms Stochastic terms for each model, by default
#' \code{ stochastic_terms  =  c( 0, 0.1, 0.3, 0.7, 1 ) }
#' @param rng Range for each variable, by default
#' \code{rng  =  c( 0,10 ) }
#'
#' @return
#' @export
#'
#' @examples
#' NULL
experiment_models  <-  function( file_name = 'output.txt', 
                                 models = c( 'Gaussian', 'Linear' ),
                                 dimensions = (1:20)*2, 
                                 stochastic_terms  =  c( 0, 0.1, 0.3, 0.7, 1, 1.5 ),
                                 rng  =  c( 0,10 ), 
                                 restrict_points_number = 1000 ){
    
    ### Check installation of libraries:
    check_packages()
    
    # delete old file
    if ( file.exists( file_name) ) unlink( file_name )
    
    DF = NULL
    for( model in models ){
        for( dimension in dimensions ){
            for( stochastic_term in stochastic_terms ){
                
                input  =  NULL
                x0  =  runif( n = dimension, min = rng[1], max = rng[2] )
                Number_of_points  =  500 * dimension
                
                if ( model == 'Gaussian' ) {
                    input = Gaussian_model( d = dimension, x0 = x0, probability = TRUE, 
                                            n = Number_of_points, r = rng )
                }
                if ( model == 'Linear' ) {
                    input  =  linear_model( d = dimension, x0 = x0, probability = TRUE, 
                                            n = Number_of_points, r = rng,
                                            noise = stochastic_term )
                }
                
                if ( is.null( input ) ) stop( 'Model name is incorrect' )
                stat.sim_origin  =  input$stat.sim
                stat.obs  =  input$stat.obs
                par.sim_origin  =  input$par.sim
                rm( input )
                
                # Apply restict number of points:
                tol = restrict_points_number / nrow( stat.sim_origin )
                rej = abc::abc( target = stat.obs, param = par.sim_origin, sumstat = stat.sim_origin,
                                method = 'rejection', tol = tol )
                
                stat.sim  =  stat.sim_origin[ rej$region, ]
                par.sim   =   par.sim_origin[ rej$region, ] 
                
                psi_t  =  adjust_psi_t( par.sim = par.sim, stat.sim = stat.sim, 
                                        stat.obs = stat.obs, talkative = FALSE, 
                                        check_pos_def = FALSE, 
                                        n_best = 8, cores = 4 )
                
                ikern  =  iKernelABC( psi = psi_t$psi[1], t = psi_t$t[1], 
                                      param = par.sim, 
                                      stat.sim = stat.sim, 
                                      stat.obs = stat.obs, 
                                      talkative = FALSE, 
                                      check_pos_def = FALSE )
                
                G = matrix( data = ikern$similarity, ncol = 1 )
                DF_new  =  Get_call_all_methods(    
                                    model_name = model, 
                                    dimension  = dimension,
                                    stochastic_term = stochastic_term, 
                                    iterations  =  1:12,
                                    stat.obs = stat.obs, 
                                    stat.sim = stat.sim, 
                                    par.sim  = par.sim, 
                                    G        = G, 
                                    par.truth  =  x0 )
                DF  =  rbind( DF, DF_new )
                
                if ( file.exists( file_name ) ){
                    write.table(file = file_name, x = DF_new , append = TRUE, sep = '\t', 
                                row.names = FALSE, col.names = FALSE )
                } else {
                    write.table(file = file_name, x = DF_new , append = TRUE, sep = '\t', 
                                row.names = FALSE, col.names = TRUE )
                }
            }
        }
    }
    
    return( DF )
}