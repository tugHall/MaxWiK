
### File to organize toy experiments which code is presented in pipeline_ABC.R 


#' Function to call a method to get a parameter estimation and MSE
#'
#' @param method_name Name of a method
#' @param kernel_name Name of kernel function
#' @param dimension  Dimension of the model
#' @param iteration Iteration number of trial with the same model and other parameters
#' @param stat.obs Data frame of statistics of observation point
#' @param stat.sim Data frame of statistics of simulations
#' @param par.sim Data frame of parameters
#' @param par.truth Truth parameter value to check result of estimation
#' @param model_function Function that is used as function to calculate output for an estimated parameter
#' @param model_par List of parameters for a model function
#' @param hyper List of hyperparameters to use for a method
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
                        stat.obs, stat.sim, par.sim, par.truth, 
                        model_function  =  Gauss_function,
                        model_par = list(d = 1, x0 = 3, r = range(0,10), noise = 0, 
                                         A = 1, sigma = 1 ),
                        hyper ){

    time_start  =  Sys.time( )

    par.est  =  Get_parameter( method_name =  method_name, 
                               kernel_name =  kernel_name,
                               stat.obs    =  stat.obs, 
                               stat.sim    =  stat.sim, 
                               par.sim     =  par.sim, 
                               hyper       =  hyper )
    
    ### Get MSE 
    MSE = NULL
    if ( !is.na( par.est )[1] ){
        model_par_all  =  c( model_par, list( par.sim1 = par.est ) )
        model_par_all$noise  =  0 
        sim_est  =  do.call( model_function, model_par_all )
        MSE  =  MSE_sim(stat.obs = stat.obs, stat.sim = sim_est ) 
    }
    
    running_time  =  as.numeric( difftime(Sys.time(), time_start, units = "secs")[[1]] )
    
    return( data.frame(   method_name = method_name,
                          kernel_name = kernel_name,
                          # model_name  = model_par$name,
                          dimension   = dimension, 
                          stochastic_term = stochastic_term,
                          MSE = round( x = MSE, digits = 2 ),
                          running_time = round( x = running_time, digits = 2 ), 
                          iteration = iteration ) )
}



#' @describeIn Get_call Function to get a parameter estimation using a method
#'
#' @return \code{Get_parameter()} returns the parameter estimation
#' 
#' @export
#' 
#' @examples
#' NULL
Get_parameter  <-  function( method_name, kernel_name = '',
                             stat.obs, stat.sim, par.sim,
                             hyper ){
    
    n_min  =  50 
    par.est  =  NA 
    
    tol  =  hyper[[ 'tolerance' ]]  
    if ( tol * nrow(stat.sim) < n_min ) tol = n_min / nrow( stat.sim )
    
    if ( method_name == 'K2-ABC' & kernel_name != 'iKernel' ){
        if ( kernel_name == 'Gaussian'){
            kernel  =  rbfdot
        } else {
            kernel  =  laplacedot
        }
        par.est  =  adjust_K2_ABC( par.sim = par.sim, stat.sim = stat.sim, 
                                   stat.obs = stat.obs, kernel = kernel, 
                                   epsilon = hyper[[ kernel_name ]][[ 'epsilon' ]]
        )[[ 'par.est' ]]
    }
    
    if ( method_name == 'K2-ABC' & kernel_name == 'iKernel' ){
        
        ikern  =  iKernelABC( psi = hyper$iKernel$psi_t$psi[ 1 ], 
                              t   = hyper$iKernel$psi_t$t[ 1 ], 
                              param = par.sim, 
                              stat.sim = stat.sim, 
                              stat.obs = stat.obs, 
                              talkative = FALSE, 
                              check_pos_def = FALSE )
        
        G = matrix( data = ikern$similarity, ncol = 1 )
        
        par.est  =  adjust_K2_ABC_iKernel( par.sim = par.sim, stat.sim = stat.sim, 
                                           stat.obs = stat.obs, G = G )[[ 'par.est' ]]
    }
    
    if ( method_name == 'Rejection' ){
        
        res  =  abc( target = stat.obs, param = par.sim, sumstat = stat.sim, 
                     tol = tol, method  =  'rejection' )
        par.est  =  Get_MAP( as.data.frame( res$unadj.values ) )
    }
    
    if ( method_name == 'Loclinear' ){
        
        loclin   =  abc(   target = stat.obs, param = par.sim, sumstat = stat.sim, 
                           tol = tol, method  =  'loclinear', hcorr   =  FALSE, 
                           transf=c( "none" ) )
        
        par.est  =  Get_MAP( as.data.frame( loclin$adj.values ) )
    }
    
    if ( method_name == 'Neuralnet' ){
        
        nn   =  abc(   target = stat.obs, param = par.sim, sumstat = stat.sim, 
                       tol = tol, method  =  'neuralnet', hcorr   =  TRUE, 
                       transf=c("none","log"), lambda = 0.0001, trace = FALSE )
        
        par.est  =  Get_MAP( as.data.frame( nn$adj.values ) )
    }
    
    if ( method_name == 'Ridge' ){
        
        rdg   =  abc(   target = stat.obs, param = par.sim, sumstat = stat.sim, 
                        tol = tol, method  =  'ridge', hcorr   =  FALSE, 
                        transf=c("none","log"), kernel = 'epanechnikov' )
        
        par.est  =  Get_MAP( as.data.frame( rdg$adj.values ) ) 
    }
    
    if ( method_name == 'MaxWiK_MAP' ){
        par.est  =  get_Spider_MAP( par.sim = par.sim, stat.sim = stat.sim, stat.obs = stat.obs )
    }
    
    if ( method_name == 'MaxWiK' ){
        par.est  =  get_Spider_MAP( par.sim = par.sim, stat.sim = stat.sim, 
                                    stat.obs = stat.obs, n_best = 1 )
    }
    
    return( par.est )
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
                                    par.sim, par.truth, cores = 4,
                                    model_function  =  Gauss_function,
                                    model_par = list(d = 1, x0 = 3, r = range(0,10), noise = 0, 
                                                     A = 1, sigma = 1 ), 
                                    hyper ){
    
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
                                   par.truth  =  par.truth,
                                   model_function  =  model_function,
                                   model_par  =  model_par, 
                                   hyper      =  hyper
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
                                    # model_name  = model_par$name,
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

