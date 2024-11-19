

#' Density plot
#'
#' @param title Title of the plot
#' @param datafr1 data frame 1
#' @param datafr2 data frame 2
#' @param var.df Variables to show
#' @param obs.true True observation if so, NULL by default
#' @param best.sim The best point from a simulation if so, NULL by default
#' @param clrs Colors to plot, by default it is c( "#a9b322", "#f9b3a2", 'red', 'blue' )
#' @param alpha Transparency values for density plots
#' @param lw Line widths
#' @param lt Line types
#'
#' @returns Make and return the ggplot object of the densities of the data frames 
#' @export
#'
#' @examples
#' MaxWiK::MaxWiK_templates(dir = tempdir()) # See the templates and vignettes for usage. 
#' # Function \code{MaxWiK.ggplot.density()} is used in the MaxWiK.ABC.R and 
#' # MaxWiK.Predictor.R templates.
MaxWiK.ggplot.density  <-  function( title = '', datafr1, datafr2, var.df, 
                                 obs.true = NULL, 
                                 best.sim = NULL,  
                                 clrs = c( "#a9b322", "#f9b3a2", 'red', 'blue' ), 
                                 alpha = c(0.1, 0.4),
                                 lw = c( 0.7, 0.7 ),
                                 lt = c( 'dashed', 'dotted' ) ){
    
    pl.base = ggplot(  ) 
    
    if ( is.null(datafr1)) {
        geom.1 = NULL
    } else {
        geom.1  =    geom_density( data = datafr1, 
                                   aes_string( x = var.df ), 
                                   fill=clrs[ 1 ], color= clrs[ 1 ], 
                                   alpha=alpha[ 1 ] ) 
    }
    
    if ( is.null(datafr2)) {
        geom.2 = NULL
    } else {
        geom.2  =    geom_density( data = datafr2, 
                                   aes_string( x = var.df ), 
                                   fill=clrs[ 2 ], color= clrs[ 2 ], 
                                   alpha=alpha[ 2 ] ) 
    }
    
    if ( is.null(obs.true)) {
        l1 = NULL
    } else {
        l1 =     geom_vline( aes(xintercept= obs.true),   
                             color=clrs[ 3 ], linetype=lt[ 1 ], linewidth=lw[ 1 ])
    }
    
    if ( is.null(best.sim)) {
        l2 = NULL
    } else {
        l2 =     geom_vline( aes(xintercept= best.sim),   
                             color=clrs[ 4 ], linetype=lt[ 2 ], linewidth=lw[ 2 ])
    }
    
    return( pl.base + geom.1 + geom.2 + l1 + l2 + ggtitle( title ) )
}
