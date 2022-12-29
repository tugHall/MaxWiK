#' Example of simulation for lazy start
#'
#' @param verbose Logical type to show or do not show messages during execution
#' @param to_plot Logical type to plot or do not plot graphical results of a simulation
#' @param seed Numeric type to set seed for a simulation, if seed = NA (by default) then it will be skipped
#' @param model Name of the toy model, can be 'Gaussian' or 'linear' only
#' @param d Integer number of dimension of the model
#' @param x0 Truth value of parameter corresponding to observation point
#' @param probability Logical to apply the probabilistic distribution to input data generation 
#' @param n Integer number of points to generate for input data frame 
#' @param r Range of parameters to plot, by default \code{r = range(0, 10)}
#' @param psi Integer number. Size of each Voronoi diagram or number of areas/points in the Voronoi diagrams
#' @param t Integer number of trees in the Isolation Forest
#' @param restrict_points_number Integer number of points which will be considered in the simulation, \cr 
#' other points will be skipped by rejection ABC
#'
#' @return List of results of simulation with default values for all the parameters
#' @export
#'
#' @examples
#' NULL
#' # it takes a time for a simulation and then it will demonstrates results, \cr
#' # so, please, wait for a while
#' sim = simulation_example( verbose = FALSE , to_plot = FALSE )
simulation_example  <-  function( verbose = TRUE , to_plot = TRUE, seed = NA, 
                                  model = c('Gaussian', 'linear')[2] ,
                                  d = 2, x0 = c(3,4), probability = TRUE, 
                                  n = 1000, r = range(0, 10),
                                  psi = 8, t = 8, 
                                  restrict_points_number = 300 ){

    if ( !is.na( seed ) ) set.seed( seed = seed )

    if ( verbose ) print('This code will be executed: ')
    if ( verbose ) print( simulation_example )

    # Attach packages from import list
    check_packages()

    if ( model == 'Gaussian' ) {
        input = get_dataset_of_Gaussian_model( d = d, x0 = x0, probability = probability, 
                                n = n, r = r)
    } else {
        input  =  linear_model( d = d, x0 = x0, probability = probability, 
                                 n = n, r = r)
    }
    stat.sim_origin  =  input$stat.sim
    stat.obs  =  input$stat.obs
    par.sim_origin  =  input$par.sim
    rm( input )
    
    tol = restrict_points_number / n 
    
    rej = abc::abc( target = stat.obs, param = par.sim_origin, sumstat = stat.sim_origin,
                    method = 'rejection', tol = tol )
    
    stat.sim  =  stat.sim_origin[ rej$region, ]
    par.sim   =   par.sim_origin[ rej$region, ]    
    
    web  =  spiderweb( psi = psi, t = t, param = par.sim_origin, 
                       stat.sim = stat.sim_origin, stat.obs = stat.obs, 
                       talkative = TRUE, check_pos_def = FALSE ,
                       n_bullets = 5, n_best = 20, halfwidth = 0.5, 
                       epsilon = 0.001 )
    
    plot_web_2D( stat.sim = stat.sim_origin , par.sim = par.sim_origin, par.truth = data.frame( t(x0) ), 
                 iKernelABC = web$iKernelABC, web = web, ind_X = 1, ind_Y = 2, 
                 names = c( 'Parameter 1', 'Parameter 2' ), xlim = c(0,10), ylim = c(0,10),
                 show_tracer = TRUE, show_obs = TRUE, show_network = TRUE, show_best = TRUE,
                 show_u_point = TRUE, show_legend = TRUE )
    title( main = ' Plot of results of spiderweb function' )
    
    readline( 'That is plot of spiderweb algorithm, press enter to show plot of sudoku algorithm' )
    
    ### The function to apply SUDOKU algorithm to get the best tracer bullets
    rslt  =  sudoku( DT = par.sim_origin , iKernelABC = web$iKernelABC, 
                      n_bullets = 5, n_best = 20, halfwidth = 0.5 )
    
    plot_sudoku_2D( stat.sim = stat.sim_origin , par.sim = par.sim_origin, par.truth = data.frame( t(x0) ), 
                    iKernelABC = web$iKernelABC, rslt = rslt, ind_X = 1, ind_Y = 2, 
                    names = c( 'Parameter 1', 'Parameter 2' ), 
                    xlim = c(0,10), ylim = c(0,10),
                    show_tracer = TRUE, show_obs = TRUE, show_appropriate = TRUE, 
                    show_best = TRUE, show_u_point = TRUE, show_legend = TRUE )
    title( main = ' Plot of results of sudoku function' )
    
    return( list( stat.sim_origin  =  stat.sim_origin, 
                  par.sim_origin  =  par.sim_origin,
                  x0 = x0,
                  stat.obs  =  stat.obs, 
                  stat.sim  =  stat.sim,
                  par.sim  =  par.sim,
                  spiderweb  = web, sudoku  =  rslt) )
}



#' @describeIn simulation_example  Example of simulation for lazy start and different psi / t hyperparameters
#'
#' @param psi_t Data.frame with different set of psi and t hyperparameters
#' @param cores Number of cores for parallel calculation
#'
#' @return List of results of simulation with default values for all the parameters
#' @export
#'
#' @examples
#' NULL
#' # it takes a time for a simulation and then it will demonstrates results, \cr
#' # so, please, wait for a while
#' sim = simulation_example_many_psi_t( verbose = FALSE , to_plot = FALSE )
simulation_example_many_psi_t  <-  function( verbose = TRUE , to_plot = TRUE, seed = NA, 
                                  model = c('Gaussian', 'linear')[2] ,
                                  d = 2, x0 = c(3,4), probability = TRUE, 
                                  n = 1000, r = range(0, 10),
                                  psi_t = data.frame( psi = c( 4,  4, 8, 8,  8,  8, 10, 10 ), 
                                                      t   = c( 8, 20, 6, 8, 16, 20,  6, 12 )  ), 
                                  restrict_points_number = 300, cores = 4 ){
    
    if ( !is.na( seed ) ) set.seed( seed = seed )
    
    if ( verbose ) print('This code will be executed: ')
    if ( verbose ) print( simulation_example )
    
    # Attach packages from import list
    check_packages()

    if ( model == 'Gaussian' ) {
        input = get_dataset_of_Gaussian_model( d = d, x0 = x0, probability = probability, 
                                n = n, r = r)
    } else {
        input  =  linear_model( d = d, x0 = x0, probability = probability, 
                                n = n, r = r)
    }
    stat.sim_origin  =  input$stat.sim
    stat.obs  =  input$stat.obs
    par.sim_origin  =  input$par.sim
    rm( input )
    
    tol = restrict_points_number / n 
    
    rej = abc::abc( target = stat.obs, param = par.sim_origin, sumstat = stat.sim_origin,
                    method = 'rejection', tol = tol )
    
    stat.sim  =  stat.sim_origin[ rej$region, ]
    par.sim   =   par.sim_origin[ rej$region, ]    
    
    webnet  =  list( )
    SIM  =  function( j ){
        psi =  psi_t$psi[ j ]
        t   =  psi_t$t[ j ]
        web  =  spiderweb_slow( psi = psi, t = t, param = par.sim_origin, 
                                stat.sim = stat.sim_origin, stat.obs = stat.obs, 
                                talkative = TRUE, check_pos_def = FALSE ,
                                n_bullets = 5, n_best = 20, halfwidth = 0.5, 
                                epsilon = 0.001 )
        return( web )
    } 
    # for( j in 1:nrow( psi_t ) ){
    webnet  =  mclapply( 1:nrow( psi_t ) , FUN = SIM, mc.cores = cores )

    
    simnet  =  list(  stat.sim_origin  =  stat.sim_origin, 
                      par.sim_origin   =  par.sim_origin,
                      x0 = x0,
                      stat.obs  =  stat.obs, 
                      stat.sim  =  stat.sim,
                      par.sim   =  par.sim,
                      psi_t     =  psi_t,
                      webnet  = webnet )
    
    simnet$networks  =  get_network_from_simnet( simnet = simnet )
    
    simnet$MAP  =  Get_MAP( simnet$networks )  # point_estimate( simnet$networks )$MAP
    
    return( simnet )
}



#' Function to extract all the best parameters estimation from results of the function 
#' \code{simulation_example_many_psi_t()}
#'
#' @param simnet results of the function \code{simulation_example_many_psi_t()}
#'
#' @return \code{get_par_best_from_simnet()} returns all the best parameters estimation 
#' 
#' 
#' @export
#'
#' @examples
#' NULL
get_par_best_from_simnet  <-  function( simnet ){
    best = NULL
    sim = NULL
    for( j in 1:length( simnet$webnet)  ){
        par.best  =  simnet$webnet[[ j ]]$par.best 
        
        best  =  rbind( best, par.best )
        sim   =  c( sim, simnet$webnet[[ j ]]$sim.best )
        
    }
    
    best$sim  =  sim
    return( best )
}


#' @describeIn get_par_best_from_simnet  Function to extract all the spiderwebs from results of the function 
#' \code{simulation_example_many_psi_t()}
#'
#'
#' @return \code{get_spiderweb_from_simnet()} returns all the spiderwebs from results of the function 
#' \code{simulation_example_many_psi_t()}
#' 
#' 
#' @export
#'
#' @examples
#' NULL
get_spiderweb_from_simnet  <-  function( simnet ){
    
    webnet  =  simnet$webnet
    param   =  simnet$webnet[[ 1 ]]$par.best   #  Get the format of parameters
    
    network =  list()
    for ( j in 1:length( webnet )){
        res = NULL
        for( p in 1:length(webnet[[ j ]]$spiderweb ) ){
            web  =  webnet[[ j ]]$spiderweb[[ p ]]
            web$iteration  =  p
            res  =  rbind( res, web )
        }
        network[[ j ]]  =  res 
    }
    
    return( network )
}



#' @describeIn get_par_best_from_simnet  Function to extract all the networks from results of the function 
#' \code{simulation_example_many_psi_t()}
#'
#'
#' @return \code{get_network_from_simnet()} returns all the networks from results of the function 
#' \code{simulation_example_many_psi_t()}
#' 
#' 
#' @export
#'
#' @examples
#' NULL
get_network_from_simnet  <-  function( simnet ){
    
    webnet  =  simnet$webnet
    param   =  simnet$webnet[[ 1 ]]$par.best   #  Get the format of parameters
    
    network =  data.frame()
    # iteration  =  NULL
    for ( j in 1:length( webnet )){

        network  =  rbind( network, webnet[[ j ]]$network )
        # iteration  =  c( iteration, rep( j, nrow( simnet$webnet[[ j ]]$network ) ) )
        
    }
    # network$iteration  =  iteration
    return( network )
}


