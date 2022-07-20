
### file to organize toy experiments

Get_call  <-  function( method_name, kernel_name, model_name, stochastics_term, itteration, args ){
    
    time_start  =  Sys.time()
    MSE = NA 
    
    if ( method_name == '' ){
        MSE  =  do.call( what = ... , args = args )
    }
    
    running_time  =  as.numeric( difftime(Sys.time(), time_start, units = "secs")[[1]] )
    
    return( data.frame(   method_name = method_name,
                          kernel_name = kernel_name,
                          model_name  = model_name,
                          stochastics_term = stochastics_term,
                          MSE = MSE, 
                          running_time = running_time, 
                          itteration = itteration ) )
}
