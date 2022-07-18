

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
#' NULL
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




#' @describeIn simulation_example Function to plot whole space of parameters 
#' related to similarity of the observation point for each set of hyperparameters psi and t 
#' 
#' @param simnet Results of the function \code{simulation_example_many_psi_t}
#' 
#' @return Many 3D plots to show the similarity of each point of the parameter space
#' for each set of hyperparameters psi and t 
#' 
#' @export
#'
#' @examples
#' NULL
plot_3d_net_similarities  <-  function( simnet, r = range(0, 10), n = 20 ){
    
    gen_x  =  ( 0 : (n-1) ) * ( r[2] - r[1] ) / ( n - 1 )
    # gen_y  =  gen_x
    
    stat.obs = simnet$stat.obs 
    par.sim = simnet$par.sim_origin
    stat.sim = simnet$stat.sim_origin
    x0 = simnet$x0
    # iKernelABC  =  sim$spiderweb$iKernelABC
    psi_t  =  simnet$psi_t
    
    
    
    # Prepare 3D plotly
    # library(plotly)
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
    
    for( j in 1:length( simnet$webnet)  ){
        
        net = generate_net( gen_x )
        web  =  simnet$webnet[[ j ]] 
        iKernelABC  =  simnet$webnet[[ j ]]$iKernelABC
        Xbest  =  as.numeric( simnet$webnet[[ j ]]$par.best )
        
        network  =  simnet$webnet[[ j ]]$network
        sim_network  =  simnet$webnet[[ j ]]$sim_network
        
        if ( FALSE ){
            iKernelABC  =  web$iKernelABC
            
            plot_web_2D( stat.sim = stat.sim, par.sim = par.sim, par.truth = data.frame( t(x0) ), 
                         iKernelABC = web$iKernelABC, web = web, ind_X = 1, ind_Y = 2, 
                         names = c( 'Parameter 1', 'Parameter 2' ), xlim = c(0,10), ylim = c(0,10),
                         show_tracer = TRUE, show_obs = TRUE, show_network = TRUE, show_best = TRUE,
                         show_u_point = TRUE, show_legend = TRUE )
        
            points( net$x1, net$x2, pch = 20, col = alpha('grey', 0.2 ), cex = 0.8 )
        
            readline( 'Press enter to continue ')
        }
        
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
        )  %>% add_trace( type = "scatter3d", mode="line",  # markers
                          name = 'Network',
                          x = network$x1, y = network$x2, z = sim_network 
        )  %>% add_trace( type = "scatter3d", mode="line", # markers
                           name = 'Estimation',
                           x = c( Xbest[1], Xbest[1] ), 
                           y = c( Xbest[2], Xbest[2] ), 
                           z = c(0, max(net$sim) )
        ) %>% layout( title = paste('Psi = ', psi_t$psi[j], ',  t = ', psi_t$t[j] ) 
        ) %>% add_trace( type = "scatter3d", mode="line", # markers
                         name = 'Truth',
                        x = c(x0[1],x0[1]), y = c(x0[2], x0[2]), z = c(0, max(net$sim))  
        )
        
        
        print( p2 )
        
        readline( 'Press enter to continue ')
    }

    return( NULL )
}



### Example of density plot of meta-sampling distribution:
if ( FALSE){
    ggplot( data = simnet$networks, aes( x = x1 ) ) + 
        geom_histogram( aes( y = ..density.. ), colour = 'white', fill = 4, alpha = 0.1) +
        geom_density(lwd = 2, colour = 4, fill = 4, alpha = 0.25)
}




