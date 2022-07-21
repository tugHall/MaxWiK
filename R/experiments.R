
### file to organize toy experiments

#' Function to call a method to get MSE
#'
#' @param method_name Name of a method
#' @param kernel_name Name of kernel function
#' @param model_name Name of a model
#' @param stochastic_term A number (usually in the range \code{[0,1]}) of stochastic term in the dataset \code{stat.sim}
#' @param iteration Iteration number of trial with the same model and other parameters
#' @param args List of arguments for a function to get MSE
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
Get_call  <-  function( method_name, kernel_name = '', model_name, stochastic_term, iteration, 
                        stat.obs, stat.sim, par.sim, G = NULL, par.truth ){
    
    time_start  =  Sys.time()
    par.est  =  NA 
    
    if ( method_name == 'K2_ABC' & kernel_name != 'iKernel' ){
        if ( kernel_name == 'Gaussian'){
            kernel  =  rbfdot
        } else {
            kernel  =  laplacedot
        }
        par.est  =  adjust_K2_ABC( par.sim = par.sim, stat.sim = stat.sim, 
                                   stat.obs = stat.obs, kernel = kernel )
    }
    
    if ( method_name == 'K2_ABC' & kernel_name == 'iKernel' ){
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
        if ( tol * nrow(stat.sim) < 5 ) tol = 5 / nrow( stat.sim )
        
        loclin   =  abc(   target = stat.obs, param = par.sim, sumstat = stat.sim, 
                           tol = tol, method  =  'loclinear', hcorr   =  FALSE, 
                           transf=c("none","log") )
        lc = round( loclin$adj.values, digits = 9 )
        
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
        if ( tol * nrow(stat.sim) < 5 ) tol = 5 / nrow( stat.sim )
        
        nn   =  abc(   target = stat.obs, param = par.sim, sumstat = stat.sim, 
                       tol = tol, method  =  'neuralnet', hcorr   =  FALSE, 
                       transf=c("none","log") )
        
        par.est  =  point_estimate( nn$adj.values )$MAP
    }
    
    if ( method_name == 'Ridge' ){
        
        res      =  adjust_ABC_tolerance( par.sim = par.sim, stat.sim = stat.sim, 
                                          stat.obs = stat.obs )
        tol   =   res$tolerance
        if ( tol * nrow(stat.sim) < 5 ) tol = 5 / nrow( stat.sim )
        
        rdg   =  abc(   target = stat.obs, param = par.sim, sumstat = stat.sim, 
                        tol = tol, method  =  'ridge', hcorr   =  FALSE, 
                        transf=c("none","log"), kernel = 'epanechnikov' )
        
        rdg_adj = round( rdg$adj.values, digits = 9 )
        
        if ( nrow( unique.data.frame( rdg_adj ) ) > 1 ){
            par.est  =  point_estimate( rdg_adj )$MAP
        } else {
            par.est  =  unique.data.frame(rdg_adj)
        }
    }
    
    if ( method_name == 'MaxWiK' ){
        par.est  =  get_Spider_MAP( par.sim = par.sim, stat.sim = stat.sim, stat.obs = stat.obs )
    }
    
    ### Get MSE 
    MSE = NULL
    if ( !is.na( par.est ) ) MSE  =  sum( ( par.thruth - par.est ) ** 2  )
    
    running_time  =  as.numeric( difftime(Sys.time(), time_start, units = "secs")[[1]] )
    
    return( data.frame(   method_name = method_name,
                          kernel_name = kernel_name,
                          model_name  = model_name,
                          stochastic_term = stochastic_term,
                          MSE = MSE, 
                          running_time = running_time, 
                          iteration = iteration ) )
}
