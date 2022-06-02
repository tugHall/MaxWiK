

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
#' @export
#'
#' @examples
#' plot_2D( x=-5:5, y=-3:7 )
plot_2D   <-  function( x, y, names = c( 'X', 'Y' ), pch = 18, col = 'blue', cex = 1.2,
                        xr = c(-10,10), yr = c(-10,10),
                        safe_pdf = FALSE, filename = './plot.pdf', change_par = TRUE ){
    sapply( X = c('grDevices', 'graphics' ), FUN = function( x ) check_pkg( x ) )
    define_par_for_plot( change_par_back = change_par )
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
#' @export
#'
#' @examples
#' DF = tugHall_dataset$data_avg
#' plot_2D_lines( x = DF[, 1 ], DF, nl = 8:12 , xr = c(1,max(DF$Time) ), yr = c(0,1) )
#' xr = c(1,max(DF$Time) )
#' yr = c(0,max(DF[,14],DF[,16],DF[,17] ))
#' plot_2D_lines( x = DF[, 1 ], DF, nl = c(14,16,17) , xr =xr, yr = yr )
#' plot_2D_lines( x = DF[, 1 ], DF, nl = 18:22 , xr = c(1,max(DF$Time) ), yr = c(0,1) )
plot_2D_lines   <-  function( x, DF, nl = 1:2, names = c( 'X', 'Y' ),
                              legend_names = '',
                               col = c( 'blue3', 'darkmagenta', 'red', 'green4',
                                        'darkorange', 'steelblue1' ),
                              cex = 1.2, lwd = 2.0,
                              lt = c(1:6), xr = c(-10,10), yr = c(-10,10),
                        safe_pdf = FALSE, filename = './plot.pdf',
                        type = 'l' , logscale = '' , draw_key  =  TRUE, change_par = TRUE ){

    sapply( X = c( 'grDevices', 'graphics' ), FUN = function( x ) check_pkg( x ) )
    define_par_for_plot( change_par_back = change_par )
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


#' Function to plot order of genes dysfunction as a step function with number of cells related to each order
#'
#' @param rdr_dysf Order of genes dysfunction as a data.frame
#' @param pos Coordinates of list of order of genes dysfunction
#' @param logscale Parameter logscale for plot function
#' @param cex Parameter cex for plot function
#' @param change_par Indicator to change par() or not for a plot. By default change_par = TRUE, after plot it will be returned to initial values
#'
#' @return NULL, making plot with step function of order of genes' dysfunction
#' @export
#'
#' @examples
#' rdr_dysf = tugHall_dataset$rdr_dysf
#' plot_order_dysfunction( rdr_dysf , logscale = '', pos = c(8, 5000), cex = 1.4)
#' plot_order_dysfunction( rdr_dysf , logscale = 'y', pos = c(10, 100), cex = 1.2)
plot_order_dysfunction  <-  function( rdr_dysf , pos = c(0,100),
                                      logscale = 'y', cex = 1, change_par = TRUE ){
    sapply( X = c( 'grDevices', 'graphics' ), FUN = function( x ) check_pkg( x ) )
    define_par_for_plot( change_par_back = change_par )
    tbl_rdr_dysf  =  aggregate( N_cells ~ order, data = rdr_dysf, FUN = sum )
    tbl_rdr_dysf  =  tbl_rdr_dysf[ order( tbl_rdr_dysf$N_cells, decreasing = T ), ]

    if ( logscale == '' ) cfcnt = 1.05 else cfcnt = 10.5

    plot_2D_lines( x = 1:nrow(tbl_rdr_dysf), DF = tbl_rdr_dysf, nl = 2 ,
                   names = c( 'Index', 'Number of cells' ),
                   yr = c(1, max( tbl_rdr_dysf$N_cells ) * cfcnt ),
                   xr = c(0.1, round( nrow( tbl_rdr_dysf )+5, digits = -1) ),
                   type = 's', logscale = logscale,  draw_key  =  FALSE )
    txt = NULL
    for( i in 1:nrow( tbl_rdr_dysf) ) {
        txt  = paste( txt, paste( i, tbl_rdr_dysf$order[ i ] ) , '\n', collapse = '   ')
    }
    text( x = pos[1], y = pos[2],
          labels = txt , pos = 4, cex = cex , col = 'red')
}


# Plot clones evolution ---------------------------------------------------


#' Function to plot clone evolution
#'
#' @param data_flow data.frame with results of simulation at each time step
#' @param threshold Vector two numbers from 0 to 1 to show clones with relative final numbers of cells in the range of threshold
#' @param lwd Line width in the plot function
#' @param hue Parameter hue in the function randomColor from library randomcoloR
#' @param luminosity Parameter luminosity in the function randomColor from library randomcoloR
#' @param yr Range for Y axes
#' @param add_initial Indicator to add or do not add initial clones to plot
#' @param log_scale Indicator to use log_scale or not for Y axes
#' @param change_par Indicator to change par() or not for a plot. By default change_par = TRUE, after plot it will be returned to initial values
#'
#' @return NULL, making plot with clones evolution
#' @export
#'
#' @examples
#' data_flow = tugHall_dataset$data_flow
#' plot_clone_evolution( data_flow, threshold = c(0.01, 1 ), add_initial = TRUE, log_scale = FALSE )
#' plot_clone_evolution( data_flow, threshold = c(0, 0.01 ), add_initial = FALSE, log_scale = TRUE )
plot_clone_evolution  <-  function( data_flow, threshold = c(0.05,1.0), lwd = 2.0,
                                    hue = c(" ", "random", "red", "orange", "yellow",
                                            "green", "blue", "purple", "pink", "monochrome")[1],
                                    luminosity = c(" ", "random", "light", "bright", "dark")[5] ,
                                    yr = NA , add_initial = TRUE, log_scale = FALSE,
                                    change_par = TRUE ){
    sapply( X = c( 'grDevices', 'graphics', 'randomcoloR' ), FUN = function( x ) check_pkg( x ) )
    define_par_for_plot( change_par_back = change_par )
    clones_flow  =  data_flow[ ,c('Time', 'ID', 'ParentID', 'Birth_time', 'N_cells' ) ]
    Nmax  =  max( clones_flow$N_cells )
    time_max  =  max( clones_flow$Time )
    Nthreshold  =  round( Nmax * threshold[1] )
    N_max  =  round( Nmax * threshold[2] )

    # delete clones with number of cells less than Nthreshold
    w = sapply( X = ( max(clones_flow$ID[ which( clones_flow$Time == 0 ) ]) + 1 ):( max(clones_flow$ID ) ),
                FUN = function( x ) {
                    if ( max( clones_flow$N_cells[ which( clones_flow$ID == x ) ] ) > Nthreshold &
                         max( clones_flow$N_cells[ which( clones_flow$ID == x ) ] ) < N_max )
                        return( x )
                    else
                        return(NULL)
                    }
                )
    w  =  unlist( w )
    if ( add_initial )  w  =  c( clones_flow$ID[ which( clones_flow$Time == 0 ) ], w )

    # clrs  =  gen_colors( nm = length( w ) )

    if (TRUE)    clrs  =  randomColor(count = length( w ),
                         hue = hue,
                         luminosity = luminosity ) else {
                             clrs = gen_colors( length( w ) )
                             clrs = clrs$color
                         }
    DF  =  list( )
    for ( i in 1:length( w ) ){
        ss  =  which( clones_flow$ID == w[ i ] )
        DF[[ i ]]  =  data.frame( x = clones_flow$Time[ ss ], y = clones_flow$N_cells[ ss ] )
    }
    if ( is.na(yr) ) yr = c( 1, N_max )
    df_pl  =  as.data.frame( DF[[ 1 ]] )
    plot_2D_lines( x = df_pl$x, DF = df_pl, nl = 2, names = c( 'Time step', 'Number of cells'),
                    xr = c( 1, round( time_max+4, digits = -1 ) ), yr = yr, draw_key = FALSE,
                   logscale = ifelse( log_scale, 'y', '' ) )

    if ( length( w ) > 1 ){
        for( i in 2:length( w ) ){
            df_pl  =  as.data.frame( DF[[ i ]] )
            lines( x = df_pl$x, y = df_pl$y, lwd = lwd, col = clrs[ i ] ) # clrs[ i, 'color']  )
        }
    }
    title( main = paste0('Number of shown clones is ', length( w ), ' from ',  max( clones_flow$ID ), ' clones' ) )
}





