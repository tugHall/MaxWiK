
### File to organize toy experiments

#' Function to call a method to get MSE
#'
#' @param method_name Name of a method
#' @param kernel_name Name of kernel function
#' @param dimension  Dimension of the model
#' @param iteration Iteration number of trial with the same model and other parameters
#' @param stat.obs Data frame of statistics of observation point
#' @param stat.sim Data frame of statistics of simulations
#' @param par.sim Data frame of parameters
#' @param G Matrix of similarities for K2-ABC based on isolation kernel
#' @param par.truth Truth parameter value to check result of estimation
#' @param model_par List of parameters for a model function
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
Get_call  <-  function( method_name, kernel_name = '', dimension, iteration, 
                        stat.obs, stat.sim, par.sim, G = NULL, par.truth, 
                        model_par = list(name = 'Gaussian', 
                                         noise = 0, sigma = 1, A = 1 ) ){
    
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
        
        par.est  =  Get_MAP( as.data.frame( loclin$adj.values ) )
        # lc = round( loclin$adj.values + runif( n= length(loclin$adj.values) )/1E6, digits = 9 )
        
        # if (nrow( unique.data.frame(lc) ) > 1 ){
        #     par.est  =  point_estimate( lc )$MAP
        # } else {
        #     par.est  =  unique.data.frame(lc)
        # }
    }
    
    if ( method_name == 'Neuralnet' ){
        
        res      =  adjust_ABC_tolerance( par.sim = par.sim, stat.sim = stat.sim, 
                                          stat.obs = stat.obs )
        tol   =   res$tolerance
        if ( tol * nrow(stat.sim) < n_min ) tol  =  n_min / nrow( stat.sim )
        
        nn   =  abc(   target = stat.obs, param = par.sim, sumstat = stat.sim, 
                       tol = tol, method  =  'neuralnet', hcorr   =  TRUE, 
                       transf=c("none","log"), lambda = 0.0001, trace = FALSE )
        
        par.est  =  Get_MAP( as.data.frame( nn$adj.values ) )
        # nn_adj = round( nn$adj.values + runif( n= length( nn$adj.values ) ) / 1E6, digits = 9 )
        # par.est  =  point_estimate( nn_adj )$MAP
    }
    
    if ( method_name == 'Ridge' ){
        
        res      =  adjust_ABC_tolerance( par.sim = par.sim, stat.sim = stat.sim, 
                                          stat.obs = stat.obs )
        tol   =   res$tolerance
        if ( tol * nrow(stat.sim) < n_min ) tol = n_min / nrow( stat.sim )
        
        rdg   =  abc(   target = stat.obs, param = par.sim, sumstat = stat.sim, 
                        tol = tol, method  =  'ridge', hcorr   =  FALSE, 
                        transf=c("none","log"), kernel = 'epanechnikov' )
        
        #rdg_adj = round( rdg$adj.values + runif( n= length( rdg$adj.values ) )/1E6, digits = 9 )
        
        par.est  =  Get_MAP( as.data.frame( rdg$adj.values ) ) 
        
        # if ( nrow( unique.data.frame( rdg_adj ) ) > 1 ){
        #    par.est  =  point_estimate( rdg_adj )$MAP
        # } else {
        #    par.est  =  unique.data.frame(rdg_adj)
        # }
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
    if ( !is.na( par.est )[1] ){
        model_par_all  =  c( model_par, list( parameter = par.est, 
                                x0 = par.truth, stat.obs = stat.obs ) )
        sim_est  =  do.call( model, model_par_all )
        MSE  =  MSE_sim(stat.obs = stat.obs, stat.sim = sim_est )  # sum( ( par.truth - par.est ) ** 2  )
    }
    
    running_time  =  as.numeric( difftime(Sys.time(), time_start, units = "secs")[[1]] )
    
    return( data.frame(   method_name = method_name,
                          kernel_name = kernel_name,
                          model_name  = model_par$name,
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
Get_call_all_methods  <-  function( dimension, iterations, stat.obs, stat.sim, 
                                    par.sim, G, par.truth, cores = 4,
                                    model_par = list(name = 'Gaussian', 
                                                     noise = 0, sigma = 1, A = 1 ) ){
    
    DF  =  NULL
    Meth_Kern  =  data.frame( Method = c('K2-ABC', 'K2-ABC', 'K2-ABC', 'Rejection', 
                                         'Loclinear', 'Neuralnet', 'Ridge',
                                         'MaxWiK_MAP', 'MaxWiK' ), 
                              Kernel = c('Gaussian', 'Laplacian', 'iKernel', '',
                                         '',          '',          'epanechnikov',
                                         'iKernel', 'iKernel') 
                              )
    
    for( mk in 1:nrow(Meth_Kern) ){
         DF_1  =  mclapply( iterations, FUN = function( x ){     
                        Get_call(  method_name =  as.character( Meth_Kern$Method[mk] ), 
                                   kernel_name =  as.character( Meth_Kern$Kernel[mk] ),
                                   dimension   =  dimension, 
                                   iteration  =  x,
                                   stat.obs   =  stat.obs, 
                                   stat.sim   =  stat.sim, 
                                   par.sim    =  par.sim, 
                                   G          =  G, 
                                   par.truth  =  par.truth,
                                   model_par  =  model_par
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
                                    model_name  = model_par$name,
                                    dimension   = dimension, 
                                    stochastic_term = model_par$noise,
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
#' @return \code{experiment_models()} returns list of results of experiments
#' @export
#'
#' @examples
#' NULL
experiment_models  <-  function( file_name = 'output.txt', 
                                 models = c( 'Gaussian', 'Linear' ),
                                 dimensions = (1:20)*2, 
                                 stochastic_terms  =  c( 0, 0.1, 0.3, 0.7, 1, 1.5 ),
                                 rng  =  c( 0,10 ), 
                                 restrict_points_number = 300 ){
    
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
                Number_of_points  =  max( c( 50 * dimension, restrict_points_number ) )
                
                if ( model == 'Gaussian' ) {
                    input = get_dataset_of_Gaussian_model( d = dimension, x0 = x0, probability = TRUE, 
                                            n = Number_of_points, r = rng,
                                            noise = stochastic_term )
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


#' Function to get statistics from the results of experiments and find the best methods for each simulation case
#'
#' @param DF Data frame with results of experiments, output of \code{experiment_models()} function.
#'
#' @return List of statistical data with analysis of results: \cr
#' - ; \cr
#' - ; \cr
#' - ; \cr
#' - ; \cr
#' - ; \cr
#' - ; \cr
#' - .
#' 
#' @export
#'
#' @examples
#' NULL 
analyze_experiments  <-  function( DF, file_to_save = '../gplot.pdf'  ){
    
    check_packages()
    
    get_region  =  function( DF, model_name, dimension, stochastic_term ){
        
        region  =  NULL
        if ( dimension != '' ){
            if ( model_name == 'Linear'){ 
                region = which(DF$model_name       ==  'Linear'        & 
                               DF$stochastic_term  ==  stochastic_term & 
                               DF$dimension        ==  dimension ) 
            }
            
            if ( model_name == 'Gaussian'){ 
                region = which(DF$model_name       ==  'Gaussian'      & 
                               DF$dimension        ==  dimension ) 
            }
        } else {
            if ( model_name == 'Linear'){ 
                region = which(DF$model_name       ==  'Linear'        & 
                                   DF$stochastic_term  ==  stochastic_term ) 
            }
            
            if ( model_name == 'Gaussian'){ 
                region = which(DF$model_name       ==  'Gaussian'  ) 
            }
        }
        if ( is.null(region) ) print( 'Model name is incorrect' )
        
        check_NA  =  complete.cases( DF[ region, ] )
        region    =  region[ check_NA ]

        return( region )
    }
    get_region_method  =  function( DF, model_name, dimension, stochastic_term, method_name ){
        
        region  =  NULL
        
        if ( model_name == 'Linear'){ 
            region = which(DF$model_name       ==  'Linear'        & 
                               DF$stochastic_term  ==  stochastic_term & 
                               DF$dimension        ==  dimension  &
                               DF$method_name      ==  method_name  ) 
        }
        
        if ( model_name == 'Gaussian'){ 
            region = which(DF$model_name       ==  'Gaussian'      & 
                               DF$dimension        ==  dimension   &
                               DF$method_name      ==  method_name ) 
        }
        
        if ( is.null(region) ) print( 'Model name is incorrect' )
        
        check_NA  =  complete.cases( DF[ region, ] )
        region    =  region[ check_NA ]
        
        return( region )
    }
    get_region_method_kernel  =  function( DF, model_name, dimension, stochastic_term, method_name, krnl ){
        
        region  =  NULL
        
        if ( model_name == 'Linear'){ 
            region = which(DF$model_name       ==  'Linear'        & 
                               DF$stochastic_term  ==  stochastic_term & 
                               DF$dimension        ==  dimension  &
                               DF$method_name      ==  method_name &
                               DF$kernel_name      ==  krnl ) 
        }
        
        if ( model_name == 'Gaussian'){ 
            region = which(DF$model_name       ==  'Gaussian'      & 
                               DF$dimension        ==  dimension   &
                               DF$method_name      ==  method_name &
                               DF$kernel_name      ==  krnl ) 
        }
        
        if ( is.null(region) ) print( 'Model name is incorrect' )
        
        check_NA  =  complete.cases( DF[ region, ] )
        region    =  region[ check_NA ]
        
        return( region )
    }
    
    # Get ranges:
    models      =  as.character( unique( DF$model_name ) )
    dimensions  =  as.numeric( unique( DF$dimension  ) )
    stochastic_terms  = as.numeric( unique( DF$stochastic_term ) )
    
    # Get the best methods for each case:
    best_methods  =  NULL
    
    for( model_name in models ){
        # if ( model_name == 'Linear') stoch = stochastic_terms else stoch  =  0
        for( dimension in dimensions ){
            for( stochastic_term in stochastic_terms ){
                
                region  =  get_region(  DF, model_name, dimension, stochastic_term )
                bm = data.frame( model_name = model_name, dimension = dimension, 
                                 stochastic_term = stochastic_term, best_method = NA , kernel = NA )
                # min_MSE  =  min( DF[ region, 'MSE'] )
                w = which.min( DF[ region, 'MSE'] )
                if ( length( w ) > 0 ){
                    bm[ 1, 'best_method']  =  as.character( DF[ region, 'method_name' ][ w ] )
                    bm[ 1, 'kernel']       =  as.character( DF[ region, 'kernel_name' ][ w ] )
                } else {
                    bm[ 1, 'best_method']  =  NA
                    bm[ 1, 'kernel']       =  NA
                }
                best_methods  =   rbind( best_methods, bm )
                
            }
        }
    }
    
    methods_all  =  unique( DF[, c(1,2)] )
    MvsD = NULL
    pdf( file_to_save )
    for( model_name in models ){
        for( stochastic_term in stochastic_terms ){
            for( j in 1:nrow( methods_all ) ){
                mthd  =  methods_all[ j, 'method_name']
                krnl  =  methods_all[ j, 'kernel_name']
                for( dimension in dimensions ){    
                    region  =  get_region_method_kernel(  DF, model_name, dimension, stochastic_term, mthd, krnl )
                    
                    if ( length( region ) >0 ){
                        # krnl  =  DF[ region, 'kernel_name'][1]
                        if ( krnl == '' ) { mwk = mthd } else { mwk = paste(mthd, krnl, sep = '_') }
                        md = data.frame( model_name = model_name, dimension = dimension, 
                                           stochastic_term = stochastic_term, 
                                           method = mthd , kernel = krnl,
                                           RMSE_min_per_dim  = sqrt(min(  DF[ region, 'MSE'] )) / dimension, 
                                           RMSE_max_per_dim  = sqrt(max(  DF[ region, 'MSE'] )) / dimension, 
                                           RMSE_mean_per_dim = sqrt(mean( DF[ region, 'MSE'] )) / dimension,
                                           method_with_kernel  =  mwk )
                    }
                MvsD  =  rbind(  MvsD, md )
                }
            }
        # Plot data frame:
        region   =  get_region(  MvsD, model_name, dimension = '', stochastic_term )
        df_plot  =  MvsD[ region, ]
        if ( model_name == 'Linear')  { 
            ttl  =  paste( 'Linear model with stochastic term', stochastic_term )
        } else {
            ttl  =  paste( 'Gaussian model with stochastic term', stochastic_term )
        }
        p = ggplot( df_plot, aes( x = dimension, y = RMSE_min_per_dim ) ) + 
                geom_line( aes( color = method_with_kernel, 
                                linetype = method_with_kernel ), size = 1.2 ) + 
            scale_color_manual( values = c( 'darkmagenta', 'blue3', 'darkgreen', 
                                            'coral4', 'deeppink', 'darkorange', 
                                            'blueviolet', 'cyan3', 'black' ) ) +
            ggtitle( ttl ) +
            ylab( 'Min of RMSE per dimension' ) + xlab( 'Dimension' )
        
        print( p )    
        
        }
    }
    dev.off()
    print( paste('All the plots are saved into file', file_to_save ) )
    return( list( best_methods = best_methods,
                  MvsD  =  MvsD ) )
}
