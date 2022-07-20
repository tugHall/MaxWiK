
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
Get_call  <-  function( method_name, kernel_name = '', model_name, stochastic_term, iteration, args ){
    
    time_start  =  Sys.time()
    par.est  =  NA 
    
    if ( method_name == 'K2_ABC' & kernel_name != 'iKernel' ){
        if ( kernel_name == 'Gaussian'){
            args$kernel  =  rbfdot
        } else {
            args$kernel  =  laplacedot
        }
        par.est  =  do.call( what = adjust_K2_ABC, args = args )
    }
    if ( method_name == 'K2_ABC' & kernel_name == 'iKernel' ){
        par.est  =  do.call( what = adjust_K2_ABC_iKernel, args = args )
    }
    
    if ( method_name == 'Rejection' ){
        par.est  =  do.call( what = adjust_ABC_tolerance$par.est, args = args )
    }
    
    if ( method_name == '' ){
        par.est  =  do.call( what = ... , args = args )
    }
    
    if ( method_name == '' ){
        par.est  =  do.call( what = ... , args = args )
    }
    
    if ( method_name == '' ){
        par.est  =  do.call( what = ... , args = args )
    }
    
    if ( method_name == '' ){
        par.est  =  do.call( what = ... , args = args )
    }
    
    
    ### Get MSE 
    if ( !is.na(par.est) ) MSE  =  
    
    running_time  =  as.numeric( difftime(Sys.time(), time_start, units = "secs")[[1]] )
    
    return( data.frame(   method_name = method_name,
                          kernel_name = kernel_name,
                          model_name  = model_name,
                          stochastic_term = stochastic_term,
                          MSE = MSE, 
                          running_time = running_time, 
                          iteration = iteration ) )
}
