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




#' @describeIn simulation_example Function to plot whole space of parameters 
#' related to similarity of the observation point:
#'
#' @param sim Results of the function \code{simulation_example}
#' @param n Integer number of points for each axes to get net for plotly 3D graph
#'
#' @return Plot of 3D to show the similarity of each point of the parameter space
#' @export
#'
#' @examples
plot_3d_similarities  <-  function( sim, r = range(0, 10), n = 12 ){
    
    gen_x  =  ( 0 : (n-1) ) * ( r[2] - r[1] ) / ( n - 1 )
    gen_y  =  gen_x
    
    stat.obs = sim$stat.obs
    par.sim = sim$par.sim_origin
    stat.sim = sim$stat.sim_origin
    x0 = sim$x0
    iKernelABC  =  sim$spiderweb$iKernelABC

    
    
    
    # Prepare 3D plotly
    library(plotly)
    generate_net  <-  function( st = (0:20)/2 ){
        # st  =  (0:20)/2
        net  =  NULL
        for( x in st ){
            net_1  =  generate_points_between_two_points(n = length( st ),
                                                         pair = data.frame( x1  =  c( x,x) , x2  = c(0,10) ) )
            net  =  rbind( net, net_1 )
        }
        return( net )
    }
    
    net = generate_net()
    web = sim$spiderweb
    iKernelABC  =  web$iKernelABC
    
    plot_web_2D( stat.sim = stat.sim, par.sim = par.sim, par.truth = data.frame( t(x0) ), 
                 iKernelABC = web$iKernelABC, web = web, ind_X = 1, ind_Y = 2, 
                 names = c( 'Parameter 1', 'Parameter 2' ), xlim = c(0,10), ylim = c(0,10),
                 show_tracer = TRUE, show_obs = TRUE, show_network = TRUE, show_best = TRUE,
                 show_u_point = TRUE, show_legend = TRUE )
    
    points( net$x1, net$x2, pch = 20, col = alpha('grey', 0.2 ), cex = 0.8 )
    
    readline( 'Press enter to continue ')
    feature_tracers  =  get_voronoi_feature_PART_dataset( data = rbind( par.sim, net ), 
                                                          talkative = TRUE, start_row = nrow( par.sim ) + 1 ,  
                                                          Matrix_Voronoi = web$iKernelABC$parameters_Matrix_Voronoi) 
    # iKernelABC$parameters_Matrix_Voronoi )
    
    sim_tracers  =  iKernel_point_dataset( Matrix_iKernel = feature_tracers$M_iKernel, 
                                           t = iKernelABC$t, nr = nrow( feature_tracers$M_iKernel ), 
                                           iFeature_point = iKernelABC$kernel_mean_embedding )
    
    # plot_ly net with similarity:
    net$sim  =  sim_tracers
    x1 = unique( net$x1 )
    x2 = unique( net$x2 )
    m = matrix( data = net$sim, ncol = length(x1), nrow = length(x2) )
    
    p2 <- plot_ly(x = x1, y = x2, z = m ) %>% add_surface(
        contours = list(
            z = list(
                show=TRUE,
                usecolormap=TRUE,
                highlightcolor="#ff0000",
                project = list( z = TRUE )
            )
        )
    ) %>% add_trace( type = "scatter3d", mode="line", # markers
                     x = c(x0[1],x0[1]), y = c(x0[2], x0[2]), z = c(0, max(net$sim)) ) 
    

    print( p2 )
    
    
    
    
    
    return( NULL )
}

