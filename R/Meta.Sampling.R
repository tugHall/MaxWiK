

# SPIDERWEB algorithm -----------------------------------------------------

#'
#' @describeIn get.MaxWiK The function to get the best value of parameter corresponding to 
#' Maxima Weighted Isolation Kernel mapping which is related to an observation point
#' 
#' @description The function \code{meta_sampling()} itteratively generates tracer based on the simple procedure: \cr
#' - making a reflection of the top points from the best point, \cr 
#' - and then generating the point tracers between them, \cr
#' - finally, the algorithm chooses again the top points and the best point (\code{sudoku()} function is used),
#' - repeat all the steps until condition to be \code{TRUE}: \cr
#' \code{abs( min( sim_tracers ) - sim_previous ) < epsilon }
#'
#'
#' @param n_bullets Number of generating points between two 
#' @param n_best Number of the best points to construct the next web net
#' @param halfwidth Parameter for the algorithm of deleting of generated points
#' @param epsilon Criterion to stop meta-sampling
#' @param rate Rate to renew points in the web net of generated points
#' @param max_iteration Maximum of iterations during meta-sampling
#' @param save_web Logical to save all the generated points (web net)
#' @param use.iKernelABC The iKernelABC object to use for meta-sampling. By default it is NULL and is generated.
#' 
#' @return The function \code{meta_sampling()} returns the list of the next objects:
#' - input.parameters that is the list of all the input parameters for Isolation Kernel ABC method;
#' - iteration that is iteration value when algorithm stopped;
#' - network that is network points when algorithm stopped;
#' - par.best that is data frame of one point that is the best from all the generated tracer points;
#' - sim.best that is numeric value of the similarity of the best tracer point;
#' - iKernelABC that is result of the function \code{get.MaxWiK()} given on \code{input parameters};
#' - spiderweb that is the list of all the networks during the meta-sampling.
#' 
#' @export
#' 
#' @examples
#' NULL
meta_sampling  <-  function(  psi = 4, t = 35, 
                              param, stat.sim, stat.obs, 
                              talkative = FALSE, check_pos_def = FALSE ,
                              n_bullets = 16, n_best = 10, halfwidth = 0.5, 
                              epsilon = 0.001, rate = 0.1, 
                              max_iteration = 15, save_web = TRUE, 
                              use.iKernelABC = NULL ){
    
    input.parameters  =  list( psi = psi, t = t, param = param, 
                               stat.sim = stat.sim, stat.obs = stat.obs, 
                               talkative = talkative, check_pos_def = check_pos_def,
                               n_bullets = n_bullets, n_best = n_best, 
                               halfwidth = halfwidth, epsilon = epsilon, rate = rate,
                               max_iteration = max_iteration, save_web = save_web )
    if ( is.null( use.iKernelABC ) ){
        data.iKernelABC  = get.MaxWiK( psi = psi, t = t, param = param, 
                                           stat.sim = stat.sim, stat.obs = stat.obs, 
                                            talkative = talkative, check_pos_def = check_pos_def )
    } else {
        data.iKernelABC  =  use.iKernelABC
    }
    ### The function to apply SUDOKU algorithm to get the best tracer bullets
    rslt  =  sudoku( DT = param , iKernelABC = data.iKernelABC, 
                     n_bullets = n_bullets, n_best = n_best, halfwidth = halfwidth )
    
    ### Get subset of parameters of Voronoi sites corresponding to an observation point
    network  =  get_subset_of_feature_map( dtst  =  param, 
                                           Matrix_Voronoi = data.iKernelABC$parameters_Matrix_Voronoi, 
                                           iFeature_point = data.iKernelABC$kernel_mean_embedding )
    network  =  unique.data.frame( network )
    
    # Number of Voronoi sites and number of rows for pool of points for network
    N_Voronoi  =  nrow( network )
    n_network  =   round( N_Voronoi * ( 1 + rate ) ) 
    
    # To add additional rows:
    n_add      =   n_network  -   nrow( network )
    
    ### Get the top:
    tracers = rslt$tracer_bullets
    ### Add the top of tracers:
    network[ N_Voronoi + ( 1 : n_add ),  ]  =  tracers[ order( rslt$similarity_to_mean , 
                                                               decreasing = TRUE)[ 1 : n_add ], ]
    
    par.best  =  tracers[ order( rslt$similarity_to_mean , decreasing = TRUE)[1], ]
    rm( tracers )
    sim_previous  =  0  
    
    sim.best  =  -1
    iteration  =  0 
    spiderweb  =  NULL
    
    ### START META-SAMPLING via iterations
    while( TRUE ){
        iteration  =  iteration  +  1
        
        if ( save_web ) spiderweb[[ iteration ]]  =  network 
        
        ### Reflect network through par.best 
        par.reflect  =  network  #  including par.best
        for( i in 1:n_network )  par.reflect[ i, ]  =  2 * par.best - network[ i ,  ] 
        
        ### Generate points between par.top and par.reflect:
        tracers  =  rbind( network, par.reflect )
        tracers  =  unique.data.frame( tracers )
        
        ### Algorithm to get diverse generation of points:
        
        gen_tr  =  get_tracer_bullets( DF = tracers, n_bullets = n_bullets )
        tracers  =  rbind( tracers, gen_tr )
        tracers  =  unique.data.frame( tracers )
        
        ### calculate the similarity for all the new points:
        feature_tracers  =  get_voronoi_feature_PART_dataset( data = rbind( param, tracers ), 
                                                              talkative = talkative, start_row = nrow( param ) + 1 ,  
                                                              Matrix_Voronoi = data.iKernelABC$parameters_Matrix_Voronoi )
        
        sim_tracers  =  iKernel_point_dataset( Matrix_iKernel = feature_tracers$M_iKernel, 
                                               t = data.iKernelABC$t, nr = nrow( feature_tracers$M_iKernel ), 
                                               iFeature_point = data.iKernelABC$kernel_mean_embedding )
        
        # new best point and top points (the cut tracers)
        rdr  =  order( sim_tracers, decreasing = TRUE )
        tracers   =   tracers[ rdr[ 1:n_network ],  ]  # Sort all the tracers and leave only the n_network top
        par.best  =   tracers[ 1, ]
        sim.best  =   sim_tracers[ rdr[1] ]
        
        ### Combine with network
        ## 1. Get similarities for network dataset
        feature_network  =  get_voronoi_feature_PART_dataset( data = rbind( param, network ), 
                                                              talkative = talkative, start_row = nrow( param ) + 1 ,  
                                                              Matrix_Voronoi = data.iKernelABC$parameters_Matrix_Voronoi )
        
        sim_network  =  iKernel_point_dataset( Matrix_iKernel = feature_network$M_iKernel, 
                                               t = data.iKernelABC$t, nr = nrow( feature_network$M_iKernel ), 
                                               iFeature_point = data.iKernelABC$kernel_mean_embedding )
        
        rdr  =  order( sim_network, decreasing = TRUE )[1:N_Voronoi]
        network   =    network[ rdr , ]  #  Cut and remove the n_add worst rows
        
        network   =    unique.data.frame( rbind( network, tracers ) )[ 1:n_network, ]  #  Renew data in the spiderweb
        rm( tracers )
        
        sim_network  =  sim_network[ rdr  ]
        sim.slow     =  sim_network[ N_Voronoi ]  # Get the worst value
        
        ### Check for break from iterations:
        if ( ( abs( sim.slow - sim_previous ) < epsilon ) | ( iteration >= max_iteration ) )   break
        sim_previous   =   sim.best
    }
    
    ### Calculate similarities of final network:
    feature_network  =  get_voronoi_feature_PART_dataset( data = rbind( param, network ), 
                                                          talkative = talkative, start_row = nrow( param ) + 1 ,  
                                                          Matrix_Voronoi = data.iKernelABC$parameters_Matrix_Voronoi )
    
    sim_network  =  iKernel_point_dataset( Matrix_iKernel = feature_network$M_iKernel, 
                                           t = data.iKernelABC$t, nr = nrow( feature_network$M_iKernel ), 
                                           iFeature_point = data.iKernelABC$kernel_mean_embedding )
    
    return( list( input.parameters = input.parameters, 
                  iteration  =  iteration, 
                  network  =  network, 
                  sim_network  =  sim_network,
                  par.best = par.best,
                  sim.best = sim.best, 
                  iKernelABC = data.iKernelABC, 
                  spiderweb = spiderweb ) )
}




# Predictor ---------------------------------------------------------------

#'
#' @describeIn get.MaxWiK The function to get the prediction of output based on a new parameter and MaxWiK
#' 
#' @description The function \code{MaxWiK.predictor()} uses the meta-sampling for a prediction
#'
#' @return The function \code{MaxWiK.predictor()} returns the list of the next objects:
#' - input.parameters that is the list of all the input parameters for Isolation Kernel ABC method;
#' - iteration that is iteration value when algorithm stopped;
#' - network that is network points when algorithm stopped;
#' - prediction.best that is data frame of one point that is the best from all the generated tracer points;
#' - sim.best that is numeric value of the similarity of the best tracer point;
#' - iKernelABC that is result of the function \code{get.MaxWiK()} given on \code{input parameters};
#' - spiderweb that is the list of all the networks during the meta-sampling.
#' 
#' @export
#' 
#' @examples
#' NULL
MaxWiK.predictor  <-  function( psi = 4, t = 35, 
                                param, stat.sim, new.param, 
                                talkative = FALSE, check_pos_def = FALSE ,
                                n_bullets = 16, n_best = 10, halfwidth = 0.5, 
                                epsilon = 0.001, rate = 0.1, 
                                max_iteration = 15, save_web = TRUE, 
                                use.iKernelABC = NULL ){

    input.parameters  =  list( psi = psi, t = t, param = param, 
                               stat.sim = stat.sim, new.param = new.param, 
                               talkative = talkative, check_pos_def = check_pos_def,
                               n_bullets = n_bullets, n_best = n_best, 
                               halfwidth = halfwidth, epsilon = epsilon, rate = rate,
                               max_iteration = max_iteration, save_web = save_web,
                               use.iKernelABC = use.iKernelABC )
    
    # !!! Replace the parameters and simulation output to use the same meta-sampling algorithm
    meta.sampling  =  meta_sampling( psi           = psi, 
                                     t             = t, 
                                     param         = stat.sim, 
                                     stat.sim      = param, 
                                     stat.obs      = new.param, 
                                     talkative     = talkative, 
                                     check_pos_def = check_pos_def,
                                     n_bullets     = n_bullets, 
                                     n_best        = n_best, 
                                     halfwidth     = halfwidth, 
                                     epsilon       = epsilon, 
                                     rate          = rate, 
                                     max_iteration = max_iteration, 
                                     save_web      = save_web, 
                                    use.iKernelABC = use.iKernelABC 
                                    )

    
    return( list( input.parameters = input.parameters, 
                  iteration    = meta.sampling$iteration, 
                  network      = meta.sampling$network, 
                  sim_network  = meta.sampling$sim_network,
                  prediction.best     = meta.sampling$par.best,
                  sim.best     = meta.sampling$sim.best, 
                  iKernelABC   = meta.sampling$iKernelABC, 
                  spiderweb    = meta.sampling$spiderweb ) )
}

