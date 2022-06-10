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
#' sim = simulation_example( verbose = FALSE , to_plot = TRUE )
simulation_example  <-  function( verbose = TRUE , to_plot = TRUE, seed = NA, 
                                  model = c('Gaussian', 'linear')[1] ,
                                  d = 2, x0 = c(3,4), probability = TRUE, 
                                  n = 1000, r = range(0, 10),
                                  psi = 8, t = 12, 
                                  restrict_points_number = 300 ){

    if ( !is.na( seed ) ) set.seed( seed = seed )

    if ( verbose ) print('This code will be executed: ')
    if ( verbose ) print( simulation_example )

    # Attach packages from import list
    packages  =  list(  randomcoloR = 'randomColor',
                        scales = 'alpha',
                        methods = 'new',
                        stats = c('aggregate', 'rbinom', 'rexp', 'rnorm', 'runif', 'dist' ),
                        utils = c('read.delim', 'read.table', 'write.table', 'globalVariables' ),
                        grDevices = c('dev.off', 'pdf', 'rgb'),
                        graphics = c('axis', 'legend', 'lines', 'par', 'plot', 'text', 'title', 'points' ),
                        abc = 'abc' )

    for( pck in names( packages ) ){
        library( package = pck, character.only = TRUE, include.only = packages[[ pck ]])
    }

    if ( model == 'Gaussian' ) {
        input = Gaussian_model( d = d, x0 = x0, probability = probability, 
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
    
    stat.sim  =  stat.sim_origin[ abc$region, ]
    par.sim   =   par.sim_origin[ abc$region, ]    
    
    web  =  spiderweb( psi = psi, t = t, param = par.sim, 
                       stat.sim = stat.sim, stat.obs = stat.obs, 
                       talkative = TRUE, check_pos_def = FALSE ,
                       n_bullets = 5, n_best = 20, halfwidth = 0.5, 
                       epsilon = 0.001 )
    
    plot_web_2D( stat.sim = stat.sim , par.sim = par.sim, par.truth = data.frame( t(x0) ), 
                 iKernelABC = web$iKernelABC, web = web, ind_X = 1, ind_Y = 2, 
                 names = c( 'Parameter 1', 'Parameter 2' ), xlim = c(0,10), ylim = c(0,10),
                 show_tracer = TRUE, show_obs = TRUE, show_top = TRUE, show_best = TRUE,
                 show_u_point = TRUE, show_legend = TRUE )
    
    readline( 'That is plot of spiderweb algorithm, press enter to show plot of sudoku algorithm' )
    
    ### The function to apply SUDOKU algorithm to get the best tracer bullets
    rslt  =  sudoku( DT = par.sim , iKernelABC = web$iKernelABC, 
                      n_bullets = 5, n_best = 20, halfwidth = 0.5 )
    
    plot_sudoku_2D( stat.sim = stat.sim , par.sim = par.sim, par.truth = data.frame( t(x0) ), 
                    iKernelABC = web$iKernelABC, rslt = rslt, ind_X = 1, ind_Y = 2, 
                    names = c( 'Parameter 1', 'Parameter 2' ), 
                    xlim = c(0,10), ylim = c(0,10),
                    show_tracer = TRUE, show_obs = TRUE, show_appropriate = TRUE, 
                    show_best = TRUE, show_u_point = TRUE, show_legend = TRUE )
    
    return( list( stat.sim  =  stat.sim, stat.obs  =  stat.obs, par.sim  =  par.sim,
                  spiderweb  = web, sudoku  =  rslt) )
}
