

# Libraries ---------------------------------------------------------------
# library( randomcoloR )
# library(stringr)


# Plot 2D --------------------------------------------------------------------

#' Function to plot 2D figure of points y = y(x)
#'
#' @param x Input data for axes X
#' @param y Input data for axes Y
#' @param names Vector of two characters with names for X and Y axes
#' @param pch Parameter pch for plot function
#' @param col Colors of points
#' @param cex Parameter cex for plot function
#' @param xr Range for X
#' @param yr Range for Y
#' @param safe_pdf Indicator to save plot to a file or not
#' @param filename Name of file to save plot if safe_pdf == TRUE
#' @param change_par Indicator to change par() or not for a plot. By default change_par = TRUE, after plot it will be returned to initial values
#'
#' @return NULL, making 2D plot using points
#' 
#'
#' @examples
#' plot_2D( x=-5:5, y=-3:7 )
plot_2D   <-  function( x, y, names = c( 'X', 'Y' ), pch = 18, col = 'blue', cex = 1.2,
                        xr = c(-10,10), yr = c(-10,10),
                        safe_pdf = FALSE, filename = './plot.pdf', change_par = TRUE ){
    sapply( X = c('grDevices', 'graphics' ), FUN = function( x ) check_pkg( x ) )
    # define_par_for_plot( change_par_back = change_par )
    rp = 1
    if ( safe_pdf )    {
        pdf( filename, width = 8, height = 8 )
        rp = 2
    }
    for( i in 1:rp ){
        par( mgp = c(2.2, 0.5, 0), font.axis = 2, font.lab = 2 )
        plot( x, y, pch = pch, xlab = names[1], xlim = xr, ylim = yr,
              ylab = names[2], axes = FALSE, cex.lab = 1.5, col = col,
              cex = cex )

        axis( 1, font = 2, tck = 0.03, cex.axis = 1.4 )
        axis( 2, font = 2, tck = 0.03, cex.axis = 1.4)
        axis( 3, font = 2, tck = 0.03, cex.axis = 1.4, labels = FALSE )
        axis( 4, font = 2, tck = 0.03, cex.axis = 1.4, labels = FALSE )

        if ( safe_pdf && i == 1 )      dev.off( )
    }
}


#' Function to plot 2D figure of lines  \code{ y_i = DF[, nl[ i ] ] ) }, i - index
#'
#' @param x Input data for axes X
#' @param DF data.frame with data to plot
#' @param nl indexes of columns in DF to plot
#' @param names Vector of two characters with names for X and Y axes
#' @param legend_names Name of legend
#' @param col Vector of colors for lines
#' @param cex Parameter cex for plot function
#' @param lwd Vector of width of lines
#' @param lt Vector of types of lines
#' @param xr Range for X
#' @param yr Range for Y
#' @param safe_pdf Indicator to save plot to a file or not
#' @param filename Name of file to save plot if safe_pdf == TRUE
#' @param type Parameter type in plot function
#' @param logscale Parameter logscale in plot function
#' @param draw_key Indicator to draw key or not
#' @param change_par Indicator to change par() or not for a plot. By default change_par = TRUE, after plot it will be returned to initial values
#'
#' @return NULL, making 2D plot using lines
#' 
#'
#' @examples
#' NULL
#' # plot_2D_lines( x = DF[, 1 ], DF, nl = 8:12 , xr = c(1,max(DF$Time) ), yr = c(0,1) )
#' # xr = c(1,max(DF$Time) )
#' # yr = c(0,max(DF[,14],DF[,16],DF[,17] ))
#' # plot_2D_lines( x = DF[, 1 ], DF, nl = c(14,16,17) , xr =xr, yr = yr )
#' # plot_2D_lines( x = DF[, 1 ], DF, nl = 18:22 , xr = c(1,max(DF$Time) ), yr = c(0,1) )
plot_2D_lines   <-  function( x, DF, nl = 1:2, names = c( 'X', 'Y' ),
                              legend_names = '',
                               col = c( 'blue3', 'darkmagenta', 'red', 'green4',
                                        'darkorange', 'steelblue1' ),
                              cex = 1.2, lwd = 2.0,
                              lt = c(1:6), xr = c(-10,10), yr = c(-10,10),
                        safe_pdf = FALSE, filename = './plot.pdf',
                        type = 'l' , logscale = '' , draw_key  =  TRUE, change_par = TRUE ){

    sapply( X = c( 'grDevices', 'graphics' ), FUN = function( x ) check_pkg( x ) )
    # define_par_for_plot( change_par_back = change_par )
    rp = 1
    if ( safe_pdf )    {
        pdf( filename, width = 8, height = 8 )
        rp = 2
    }
    for( i in 1:rp ){
        par( mgp = c(2.2, 0.5, 0), font.axis = 2, font.lab = 2 )
        ### Plot the first line:
        y = DF[, nl[1] ]
        plot( x, y, xlab = names[1], xlim = xr, ylim = yr,
              ylab = names[2], axes = FALSE, cex.lab = 1.5, col = col[1],
              lwd = lwd, lty = lt[1],
              type = type, log = logscale )

        axis( 1, font = 2, tck = 0.03, cex.axis = 1.4 )
        axis( 2, font = 2, tck = 0.03, cex.axis = 1.4)
        axis( 3, font = 2, tck = 0.03, cex.axis = 1.4, labels = FALSE )
        axis( 4, font = 2, tck = 0.03, cex.axis = 1.4, labels = FALSE )

        ### More than 1 line:
        if ( length( nl ) > 1 ){
            for( j in 2:length(nl) ){
                st = nl[j]
                lines( x, DF[, st ], lwd = lwd, col = col[j], lty = lt[j] )
            }
        }
        if ( draw_key ){
            key = names( DF[ nl ] )
            if( legend_names[1] != '') key = legend_names
            legend( x = 'bottom', legend = key,
                    horiz = TRUE, xpd = TRUE,  inset = c(0, 1.03), box.col = "white",
                    lty = lt[ 1:length(nl) ], col = col[ 1:length(nl) ]  )
        }
        if ( safe_pdf & i == 1 )      dev.off( )
    }
}




# PLOT SUDOKU AND SPIDERWEB--------------------------------------------------------------------

#' Function to plot all the results of the function \code{sudoku()}
#' 
#' @description The function \code{plot_sudoku_2D} draws all the results of the function \code{sudoku()}
#' 
#' @param stat.sim Summary statistics of the simulations (model output)
#' @param par.sim Data frame of parameters of the model
#' @param par.truth Truth value of the parameter corresponding to observation point (if known)
#' @param iKernelABC Result of calculations based on Isolation Kernel ABC 
#' that can be gotten by the function \code{iKernelABC()}
#' @param rslt Results of function \code{sudoku()}
#' @param ind_X Column index of the par.sim data frame to plot as X axes 
#' @param ind_Y Column index of the par.sim data frame to plot as Y axes 
#' @param names Vector of axes names, by default \code{names = c( 'Parameter_1', 'Parameter_2' )}
#' @param xlim Numeric vector of the range for X axes
#' @param ylim Numeric vector of the range for Y axes
#' @param show_tracer Logical to show or do not show all the tracer points
#' @param show_obs Logical to show or do not show the stat.obs or parameter of observation point (if known)
#' @param show_appropriate Logical to show or do not show all the points of 
#' the object \code{rslt$surroundings_best_points} from results of \code{sudoku()} function
#' @param show_best Logical to show or do not show the best point of the results of \code{sudoku()} function
#' @param show_u_point Logical to show or do not show the estimated point of parameter for an observation 
#' @param show_legend Logical to show or do not show the legend 
#'
#' @return Plot all the results of \code{sudoku()} function 
#' 
#' 
#' 
#' @examples
#' NULL
plot_sudoku_2D   <-  function( stat.sim, par.sim, par.truth, iKernelABC, rslt, 
                               ind_X, ind_Y, names = c( 'Parameter_1', 'Parameter_2' ), 
                               xlim, ylim,
                               show_tracer = TRUE, show_obs = TRUE, show_appropriate = TRUE, show_best = TRUE,
                               show_u_point = TRUE, show_legend = FALSE ){
    
    sbst_feature_Param  =  get_subset_of_feature_map( dtst  =  par.sim,
                                                      Matrix_Voronoi = iKernelABC$parameters_Matrix_Voronoi,
                                                      iFeature_point = iKernelABC$kernel_mean_embedding )
    # if ( draw_Y ){
    #    sbst_feature_Y  =  get_subset_of_feature_map( dtst  =  stat.sim,
    #                                              Matrix_Voronoi = iKernelABC$Matrix_Voronoi,
    #                                              iFeature_point = iKernelABC$iFeature_point )
    # }
    # pr = par()
    
    l = par.sim[ , c( ind_X, ind_Y ) ]
    par( mgp = c(2.2, 0.5, 0), font.axis = 2, font.lab = 2 )
    plot( l[,1], l[,2], pch = 18, xlab = names[1], ylab = names[2], 
          axes = FALSE, cex.lab = 1.5, xlim = xlim, ylim = ylim )
    
    # cex.lab = 1.3 , cex.axis = 1.4 , mar = c(4,4,2,2),
    # tck = 0.04) ,
    # box()
    
    axis( 1, font = 2, tck = 0.03, cex.axis = 1.4 )
    axis( 2, font = 2, tck = 0.03, cex.axis = 1.4)    
    axis( 3, font = 2, tck = 0.03, cex.axis = 1.4, labels = FALSE )
    axis( 4, font = 2, tck = 0.03, cex.axis = 1.4, labels = FALSE ) 
    
    
    
    # if ( draw_Y )    points( sbst_feature_Y[,ind_X], sbst_feature_Y[,ind_Y], col = 'red', cex = 1.6 )
    if ( show_u_point ) points( sbst_feature_Param[,ind_X], sbst_feature_Param[,ind_Y], col = 'blue', cex = 1.4, pch = 0 )
    if ( show_tracer ) points( rslt$tracer_bullets[,ind_X], rslt$tracer_bullets[,ind_Y], col = alpha('grey', 0.2 ), cex = 0.8, pch = 20 )
    if ( show_appropriate ) points( rslt$surroundings_best_points[,ind_X], rslt$surroundings_best_points[,ind_Y], 
                                    col = 'blue', cex = 0.5, pch = 20 )
    if ( show_best ) points( rslt$best_tracer_bullets[,ind_X], rslt$best_tracer_bullets[,ind_Y], col = 'red', cex = 0.8, pch = 5 )
    if (show_obs ) points( par.truth[,ind_X], par.truth[,ind_Y], col = 'red', cex = 2.4, pch = 4, lwd = 3 )
    
    if ( show_legend ) legend( x = 7.7, y = 10.3, bg = '#E0FFFF', box.col = 'grey',
                               col = c('black', 'red', 'blue', 'grey', 'blue', 'red' ), 
                               pch = c( 18,      4,       0 ,  1,    20,   5  ), bty = "o", 
                               legend = c('Simulation','Truth observation', 
                                          expression( paste('Voronoi sites ', upsilon[j]^'*' ) ),
                                          'Tracers', 'iKernel > 0.5', 
                                          'Top 20 points' )  
    )
    
}



#' @describeIn plot_sudoku_2D The function to plot all the results of the function \code{spiderweb()}
#'
#' @description The function \code{plot_web_2D()} draws all the results of the function \code{spiderweb()}
#'
#' @param web Results of the function \code{spiderweb()} to draw on the plot
#' @param show_network Logical to show or do not show network points
#'
#' @return The function \code{plot_web_2D()} draws all the results of the function \code{spiderweb()}
#' 
#' 
#'
#' @examples
#' NULL
plot_web_2D   <-  function( stat.sim, par.sim, par.truth, iKernelABC, web, ind_X, ind_Y, 
                            names = c( 'P1', 'P2' ),
                            xlim, ylim,
                            show_tracer = TRUE, show_obs = TRUE, show_network = TRUE, show_best = TRUE,
                            show_u_point = TRUE, show_legend = FALSE ){
    
    sbst_feature_Param  =  get_subset_of_feature_map( dtst  =  par.sim,
                                                      Matrix_Voronoi = iKernelABC$parameters_Matrix_Voronoi,
                                                      iFeature_point = iKernelABC$kernel_mean_embedding )
    # if ( draw_Y ){
    #     sbst_feature_Y  =  get_subset_of_feature_map( dtst  =  stat.sim,
    #                                                   Matrix_Voronoi = iKernelABC$Matrix_Voronoi,
    #                                                   iFeature_point = iKernelABC$iFeature_point )
    # }
    # pr = par()
    
    l = par.sim[ , c( ind_X, ind_Y ) ]
    par( mgp = c(2.2, 0.5, 0), font.axis = 2, font.lab = 2 )
    plot( l[,1], l[,2], pch = 18, xlab = names[1], ylab = names[2], 
          axes = FALSE, cex.lab = 1.5, xlim = xlim, ylim = ylim )
    
    # cex.lab = 1.3 , cex.axis = 1.4 , mar = c(4,4,2,2),
    # tck = 0.04) ,
    # box()
    
    axis( 1, font = 2, tck = 0.03, cex.axis = 1.4 )
    axis( 2, font = 2, tck = 0.03, cex.axis = 1.4)    
    axis( 3, font = 2, tck = 0.03, cex.axis = 1.4, labels = FALSE )
    axis( 4, font = 2, tck = 0.03, cex.axis = 1.4, labels = FALSE ) 
    
    
    
    # if ( draw_Y )    points( sbst_feature_Y[,ind_X], sbst_feature_Y[,ind_Y], col = 'red', cex = 1.6 )
    if ( show_u_point ) points( sbst_feature_Param[,ind_X], sbst_feature_Param[,ind_Y], col = 'blue', cex = 1.4, pch = 0 )
    # if ( show_tracer ) points( web$tracers_all[,ind_X], web$tracers_all[,ind_Y], col = alpha('grey', 0.2 ), cex = 0.8, pch = 20 )
    if ( show_network ) points( web$network[,ind_X], web$network[,ind_Y], col = 'blue', cex = 0.5, pch = 20 )
    if ( show_best ) points( web$par.best[,ind_X], web$par.best[,ind_Y], col = 'red', cex = 0.8, pch = 5 )
    if (show_obs ) points( par.truth[,ind_X], par.truth[,ind_Y], col = 'red', cex = 2.4, pch = 4, lwd = 3 )
    
    if ( show_legend ) legend( x = 7.7, y = 10.3, bg = '#E0FFFF', box.col = 'grey',
                               col = c('black', 'red', 'blue', 'grey', 'blue', 'red' ), 
                               pch = c( 18,      4,       0 ,  1,    20,   5  ), bty = "o", 
                               legend = c('Simulation','Truth observation', 
                                          expression( paste('Voronoi sites ', upsilon[j]^'*' ) ),
                                          'Tracers', 'Top 20 points', 'The best point' )  
    )
    
}

