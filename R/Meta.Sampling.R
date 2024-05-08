

# SPIDERWEB algorithm -----------------------------------------------------

#' @describeIn get.iKernelABC The function to get the best value of parameter corresponding to 
#' Maxima Weighted Isolation Kernel mapping which is related to an observation point
#' 
#' @description The function \code{spiderweb()} itteratively generates tracer points gotten 
#' from \code{sudoku()} algorithm, based on the simple procedure: \cr
#' - making a reflection of the top points from the best point, \cr 
#' - and then generating the point tracers between them, \cr
#' - finally, the algorithm chooses again the top points and the best point (\code{sudoku()} function is used),
#' - repeat all the steps until condition to be \code{TRUE}: \cr
#' \code{abs( max( sim_tracers ) - sim_previous ) < epsilon }
#'
#' @param n_bullets Integer number of tracer bullets / additional points between the TWO most distant points
#' @param n_best Integer number of the best tracer bullets / points 
#' to consider them at the next algorithmic step
#' @param halfwidth Criterion to choose the best tracer points like: \cr
#' \code{if similarity_of_point >= halfwidth} then it is the point to be included to the poool of the best points
#' @param epsilon Criterion to stop algorithm \code{spiderweb()} that isused to check: \cr
#' \code{if ( abs( max( sim_tracers ) - sim_previous ) < epsilon ) break}
#'
#' @return The function \code{spiderweb()} returns the list of the next objects:
#' - input.parameters the list of all the input parameters for Isolation Kernel ABC method;
#' - par.best that is data frame of one point that is the best from all the generated tracer points;
#' - par.top that is data frame of n_best points that are the top from all the generated tracer points; 
#' - sim.top that is numeric vecor of similarities of the top points;
#' - sim.best that is numeric value of the similarity of the best tracer point;
#' tracers_all that is data frame of all the generated tracer points; 
#' - sim.tracers_all that is numeric vector of similarities of all the generated tacer points;
#' - iKernelABC that is result of the function \code{get.iKernelABC()} given on \code{input parameters}.
#' 
#' 
#' 
#' @examples
#' NULL
spiderweb_old  <-  function( psi = 4, t = 35, param = param, 
                             stat.sim = stat.sim, stat.obs = stat.obs, 
                             talkative = FALSE, check_pos_def = FALSE ,
                             n_bullets = 5, n_best = 10, halfwidth = 0.5, 
                             epsilon = 0.001 ){
    
    input.parameters  =  list( psi = psi, t = t, param = param, 
                               stat.sim = stat.sim, stat.obs = stat.obs, 
                               talkative = talkative, check_pos_def = check_pos_def,
                               n_bullets = n_bullets, n_best = n_best, 
                               halfwidth = halfwidth, epsilon = epsilon )
    
    iKernelABC  = get.iKernelABC( psi = psi, t = t, param = param, 
                                    stat.sim = stat.sim, stat.obs = stat.obs, 
                                    talkative = talkative, check_pos_def = check_pos_def )
    
    ### The function to apply SUDOKU algorithm to get the best tracer bullets
    rslt  =  sudoku( DT = param , iKernelABC = iKernelABC, 
                     n_bullets = n_bullets, n_best = n_best, halfwidth = halfwidth )
    
    ### Get the top:
    tracers = rslt$tracer_bullets
    
    tracers_all  =  tracers
    sim.tracers_all  =  rslt$similarity_to_mean
    
    par.top   =  tracers[ order( rslt$similarity_to_mean , decreasing = TRUE)[1:n_best], ]
    par.best  =  par.top[ 1, ]
    par.top   =  par.top[2:n_best, ]
    rm( tracers )
    sim_previous  =  0  # max( rslt$similarity_to_mean )
    
    sim.top   =  NULL
    sim.best  =  -1
    while( TRUE ){
        ### Reflect par.top through par.best 
        par.reflect  =  par.top
        for( i in 1:nrow( par.top ) )  par.reflect[ i, ]  =  2 * par.best - par.top[ i ,  ] 
        
        ### Generate points between par.top and par.reflect:
        tracers  =  rbind( par.best, par.top, par.reflect )
        for( i in 1:nrow( par.top  ) ){
            gen_tr  =  generate_points_between_two_points( pair = rbind( par.top[ i ,] , 
                                                                         par.reflect[ i, ] ), n = n_bullets )
            tracers  =  rbind( tracers, gen_tr )
        }
        
        ### calculate the similarity for new points:
        feature_tracers  =  get_voronoi_feature_PART_dataset( data = rbind( param, tracers ), 
                                                              talkative = talkative, start_row = nrow( param ) + 1 ,  
                                                              Matrix_Voronoi = iKernelABC$parameters_Matrix_Voronoi )
        
        sim_tracers  =  iKernel_point_dataset( Matrix_iKernel = feature_tracers$M_iKernel, 
                                               t = iKernelABC$t, nr = nrow( feature_tracers$M_iKernel ), 
                                               iFeature_point = iKernelABC$kernel_mean_embedding )
        
        tracers_all  =  rbind( tracers_all, tracers )
        sim.tracers_all  =  c( sim.tracers_all, sim_tracers )
        # new best point
        par.best     =  tracers[ which.max(sim_tracers ), ]
        
        if ( abs( max( sim_tracers ) - sim_previous ) < epsilon ) break
        sim_previous   =   max( sim_tracers )
        
        ### rename new tracers:
        sim.top  =  sort( sim_tracers , decreasing = TRUE)[1:n_best]
        sim.best =  sim.top[ 1 ]
        sim.top  =  sim.top[2:n_best]
        
        par.top  =  tracers[ order( sim_tracers , decreasing = TRUE)[1:n_best], ]
        par.best  =  par.top[ 1, ]
        par.top   = par.top[2:n_best, ]
        
        rm( tracers )
    }
    
    return( list( input.parameters = input.parameters,  par.best = par.best, par.top = par.top, 
                  sim.top = sim.top, sim.best = sim.best, tracers_all = tracers_all, 
                  sim.tracers_all = sim.tracers_all, iKernelABC = iKernelABC ) )
}

#' @describeIn get.iKernelABC The function to get the best value of parameter corresponding to 
#' Maxima Weighted Isolation Kernel mapping which is related to an observation point
#' 
#' @description The function \code{spiderweb()} itteratively generates tracer points gotten 
#' from \code{sudoku()} algorithm, based on the simple procedure: \cr
#' - making a reflection of the top points from the best point, \cr 
#' - and then generating the point tracers between them, \cr
#' - finally, the algorithm chooses again the top points and the best point (\code{sudoku()} function is used),
#' - repeat all the steps until condition to be \code{TRUE}: \cr
#' \code{abs( max( sim_tracers ) - sim_previous ) < epsilon }
#'
#' @param n_bullets Integer number of tracer bullets / additional points between the TWO most distant points
#' @param n_best Integer number of the best tracer bullets / points 
#' to consider them at the next algorithmic step
#' @param halfwidth Criterion to choose the best tracer points like: \cr
#' \code{if similarity_of_point >= halfwidth} then it is the point to be included to the poool of the best points
#' @param epsilon Criterion to stop algorithm \code{spiderweb()} that isused to check: \cr
#' \code{if ( abs( max( sim_tracers ) - sim_previous ) < epsilon ) break}
#' @param rate Numeric rate from 0 to 1 that gives rate of changing of surround points of proposed max of similarity
#' or part of changing of network during meta-sampling
#' @param max_iteration Maximal number of iteration in the function
#' @param save_web Logical to save or do not save network during meta-sampling
#'
#' @return The function \code{spiderweb()} returns the list of the next objects:
#' - input.parameters that is the list of all the input parameters for Isolation Kernel ABC method;
#' - iteration that is iteration value when algorithm stopped;
#' - network that is network points when algorithm stopped;
#' - par.best that is data frame of one point that is the best from all the generated tracer points;
#' - sim.best that is numeric value of the similarity of the best tracer point;
#' - iKernelABC that is result of the function \code{get.iKernelABC()} given on \code{input parameters};
#' - spiderweb that is the list of all the networks during the meta-sampling.
#' 
#' 
#' 
#' @examples
#' NULL
spiderweb  <-  function( psi = 4, t = 35, param = param, 
                         stat.sim = stat.sim, stat.obs = stat.obs, 
                         talkative = FALSE, check_pos_def = FALSE ,
                         n_bullets = 16, n_best = 10, halfwidth = 0.5, 
                         epsilon = 0.001, rate = 0.1, 
                         max_iteration = 5, save_web = TRUE ){
    
    input.parameters  =  list( psi = psi, t = t, param = param, 
                               stat.sim = stat.sim, stat.obs = stat.obs, 
                               talkative = talkative, check_pos_def = check_pos_def,
                               n_bullets = n_bullets, n_best = n_best, 
                               halfwidth = halfwidth, epsilon = epsilon, rate = rate,
                               max_iteration = max_iteration, save_web = save_web )
    
    iKernelABC  = get.iKernelABC( psi = psi, t = t, param = param, 
                                    stat.sim = stat.sim, stat.obs = stat.obs, 
                                    talkative = talkative, check_pos_def = check_pos_def )
    
    ### The function to apply SUDOKU algorithm to get the best tracer bullets
    rslt  =  sudoku( DT = param , iKernelABC = iKernelABC, 
                     n_bullets = n_bullets, n_best = n_best, halfwidth = halfwidth )
    
    ### Get subset of parameters of Voronoi sites corresponding to an observation point
    network  =  get_subset_of_feature_map( dtst  =  param, 
                                           Matrix_Voronoi = iKernelABC$parameters_Matrix_Voronoi, 
                                           iFeature_point = iKernelABC$kernel_mean_embedding )
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
        
        if ( FALSE ){    
            ### Change IT!!!! Incorrect!!! This algorithm seeks only across the best tracer not around,
            ### This algorithm restricts diversity of points in the high dimensional space
            for( i in 1:n_network ){
                gen_tr  =  generate_points_between_two_points( pair = rbind( network[ i ,] , 
                                                                             par.reflect[ i, ] ), n = n_bullets )
                tracers  =  rbind( tracers, gen_tr )
            }
            ### Change IT!!!! Incorrect!!!
        }
        
        ### Alternative algorithm to get diverse generation of points:
        if ( TRUE ){
            
            gen_tr  =  get_tracer_bullets( DF = tracers, n_bullets = n_bullets )
            
            tracers  =  rbind( tracers, gen_tr )
            
            tracers  =  unique.data.frame( tracers )
        }
        
        ### calculate the similarity for all the new points:
        feature_tracers  =  get_voronoi_feature_PART_dataset( data = rbind( param, tracers ), 
                                                              talkative = talkative, start_row = nrow( param ) + 1 ,  
                                                              Matrix_Voronoi = iKernelABC$parameters_Matrix_Voronoi )
        
        sim_tracers  =  iKernel_point_dataset( Matrix_iKernel = feature_tracers$M_iKernel, 
                                               t = iKernelABC$t, nr = nrow( feature_tracers$M_iKernel ), 
                                               iFeature_point = iKernelABC$kernel_mean_embedding )
        
        # tracers_all  =  rbind( tracers_all, tracers )  # Change later 
        # sim.tracers_all  =  c( sim.tracers_all, sim_tracers )  # Change later 
        
        # new best point and top points (the cut tracers)
        rdr  =  order( sim_tracers, decreasing = TRUE )
        tracers   =   tracers[ rdr[ 1:n_add ],  ]  # Remove all the tracer and leave only the n_add top
        par.best  =   tracers[ 1, ]
        sim.best  =   sim_tracers[ rdr[1] ]
        
        ### Combine with network
        ## 1. Get similarities for network dataset
        feature_network  =  get_voronoi_feature_PART_dataset( data = rbind( param, network ), 
                                                              talkative = talkative, start_row = nrow( param ) + 1 ,  
                                                              Matrix_Voronoi = iKernelABC$parameters_Matrix_Voronoi )
        
        sim_network  =  iKernel_point_dataset( Matrix_iKernel = feature_network$M_iKernel, 
                                               t = iKernelABC$t, nr = nrow( feature_network$M_iKernel ), 
                                               iFeature_point = iKernelABC$kernel_mean_embedding )
        
        network   =    network[ order( sim_network, decreasing = TRUE )[1:N_Voronoi] , ]  #  Cut and remove the n_add worst rows
        
        network   =    rbind( network, tracers )  #  Renew data in the spiderweb
        rm( tracers )
        
        ### Check for break from iterations:
        if ( ( abs( sim.best - sim_previous ) < epsilon ) | ( iteration >= max_iteration ) )   break
        sim_previous   =   sim.best
    }
    
    return( list( input.parameters = input.parameters, 
                  iteration  =  iteration, 
                  network  =  network, 
                  par.best = par.best,
                  sim.best = sim.best, 
                  iKernelABC = iKernelABC, 
                  spiderweb = spiderweb ) )
}

#' @describeIn get.iKernelABC The function to get the best value of parameter corresponding to 
#' Maxima Weighted Isolation Kernel mapping which is related to an observation point
#' 
#' @description The function \code{meta_sampling()} itteratively generates tracer points gotten 
#' from \code{sudoku()} algorithm, based on the simple procedure: \cr
#' - making a reflection of the top points from the best point, \cr 
#' - and then generating the point tracers between them, \cr
#' - finally, the algorithm chooses again the top points and the best point (\code{sudoku()} function is used),
#' - repeat all the steps until condition to be \code{TRUE}: \cr
#' \code{abs( min( sim_tracers ) - sim_previous ) < epsilon }
#'
#'
#' @return The function \code{spiderweb()} returns the list of the next objects:
#' - input.parameters that is the list of all the input parameters for Isolation Kernel ABC method;
#' - iteration that is iteration value when algorithm stopped;
#' - network that is network points when algorithm stopped;
#' - par.best that is data frame of one point that is the best from all the generated tracer points;
#' - sim.best that is numeric value of the similarity of the best tracer point;
#' - iKernelABC that is result of the function \code{get.iKernelABC()} given on \code{input parameters};
#' - spiderweb that is the list of all the networks during the meta-sampling.
#' 
#' @export
#' 
#' @examples
#' NULL
meta_sampling  <-  function( psi = 4, t = 35, param = param, 
                              stat.sim = stat.sim, stat.obs = stat.obs, 
                              talkative = FALSE, check_pos_def = FALSE ,
                              n_bullets = 16, n_best = 10, halfwidth = 0.5, 
                              epsilon = 0.001, rate = 0.1, 
                              max_iteration = 15, save_web = TRUE, 
                              use.iKernelABC = NA ){
    
    input.parameters  =  list( psi = psi, t = t, param = param, 
                               stat.sim = stat.sim, stat.obs = stat.obs, 
                               talkative = talkative, check_pos_def = check_pos_def,
                               n_bullets = n_bullets, n_best = n_best, 
                               halfwidth = halfwidth, epsilon = epsilon, rate = rate,
                               max_iteration = max_iteration, save_web = save_web )
    if ( is.na( use.iKernelABC ) ){
        data.iKernelABC  = get.iKernelABC( psi = psi, t = t, param = param, 
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






