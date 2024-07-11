###############  LIBRARY and FUNCTION  ####################
library('MaxWiK')
library('ggplot2')

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
                                   aes( x =  .data[[ var.df ]] ), 
                                   fill=clrs[ 1 ], color= clrs[ 1 ], 
                                    alpha=alpha[ 1 ] ) 
    }
    
    if ( is.null(datafr2)) {
        geom.2 = NULL
    } else {
        geom.2  =    geom_density( data = datafr2, 
                                   aes( x =  .data[[ var.df ]] ), 
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
    
    print( pl.base + geom.1 + geom.2 + l1 + l2 + ggtitle( title ) )
}
# LOAD the data: --------------------------------------------------------

hyper = read_hyperparameters( input = './HYPERPARAMETERS/Hyperparameters.Dim-2.Noise-0.txt' )

Matrix.Voronoi  =  readRDS( file = hyper$files$Matrix.Voronoi )

file.nm  =  paste0( './Input.data/', hyper$dm, 'D/input.data.noise.', hyper$noise, '.txt' )
DF = read.table(file = file.nm, header = TRUE, 
                sep = '\t', stringsAsFactors = FALSE)

names(DF)
par.sim   =  DF[ , (hyper$dm+1) : (2*hyper$dm) ]

stat.sim  =  DF[ , 1 : hyper$dm ]

rm(DF)

obs = read.table(file = paste0( './Input.data/', hyper$dm , 'D/Gaussian.txt' ),
                 header = TRUE, sep = '\t', stringsAsFactors = FALSE)


# Approximate Bayesian Computation ---------------------------------------------------------------------




res.MaxWik  =  get.MaxWiK( 
                            psi = hyper$psi, 
                            t   = hyper$t, 
                            param = par.sim, 
                            stat.sim = stat.sim, 
                            stat.obs = as.data.frame( t( obs$A ) ), 
                            talkative = TRUE, 
                            check_pos_def = TRUE, 
                            Matrix_Voronoi = Matrix.Voronoi
                            )

w.sim  =  which(res.MaxWik$similarity > 0 )
posteriori.MaxWiK  =  par.sim[ w.sim, ]
posteriori.weights =  res.MaxWik$similarity[ w.sim ] / sum( res.MaxWik$similarity[ w.sim ] )


MaxWiK.ggplot.density( title = ' Posteriori distribution of X2 parameter (Red line is a true value) ', 
                   datafr1 = posteriori.MaxWiK, 
                   datafr2 = NULL, 
                   var.df  = 'par.sim.X2', 
                   obs.true = obs$x0[ 2 ], 
                   best.sim = NULL
                )

MaxWiK.ggplot.density( title = ' Posteriori distribution of X1 parameter (Red line is a true value) ', 
                   datafr1 = posteriori.MaxWiK, 
                   datafr2 = NULL, 
                   var.df  = 'par.sim.X1', 
                   obs.true = obs$x0[ 1 ], 
                   best.sim = NULL
)

##########################  META-SAMPLING  #######################################

meta.sampling  =  meta_sampling( psi = hyper$psi, t = hyper$t, param = par.sim, 
                              stat.sim = stat.sim, stat.obs = as.data.frame( t( obs$A ) ), 
                              talkative = FALSE, check_pos_def = FALSE ,
                              n_bullets = 42, n_best = 12, halfwidth = 0.5, 
                              epsilon = 0.001, rate = 0.2, 
                              max_iteration = 10, save_web = TRUE, 
                              use.iKernelABC = res.MaxWik )


network  = unique.data.frame( do.call(rbind.data.frame, meta.sampling$spiderweb ) )

MaxWiK.ggplot.density( title = ' Posteriori distribution of X1 parameter \n (Red line is a true value, dotted is the best simulation) \n red area is metasampling', 
                   datafr1 = posteriori.MaxWiK, 
                   datafr2 = network, 
                   var.df  = 'par.sim.X1', 
                   obs.true = obs$x0[ 1 ], 
                   best.sim = as.numeric( meta.sampling$par.best[ 1 ] )
        )


MaxWiK.ggplot.density( title = ' Posteriori distribution of X2 parameter \n (Red line is a true value, dotted is the best simulation) \n red area is metasampling', 
                   datafr1 = posteriori.MaxWiK, 
                   datafr2 = network, 
                   var.df  = 'par.sim.X2', 
                   obs.true = obs$x0[ 2 ], 
                   best.sim = as.numeric( meta.sampling$par.best[ 2 ] )
)


##########################  PREDICTOR  #######################################


pred.MaxWiK  =  get.MaxWiK( 
                        psi = hyper$psi, 
                        t   = hyper$t, 
                        param = stat.sim,
                        stat.sim =  par.sim,
                        stat.obs = as.data.frame( t( obs$x0 ) ), 
                        talkative = TRUE, 
                        check_pos_def = TRUE, 
                        Matrix_Voronoi = Matrix.Voronoi
                    )

w.sim  =  which(pred.MaxWiK$similarity > 0 )
posteriori.pred.MaxWiK  =  stat.sim[ w.sim, ]
posteriori.pred.weights =  pred.MaxWiK$similarity[ w.sim ] / sum( pred.MaxWiK$similarity[ w.sim ] )


MaxWiK.ggplot.density( title = ' Posteriori distribution of Y1 parameter \n (Red line is a true value)', 
                   datafr1 = posteriori.pred.MaxWiK, 
                   datafr2 = NULL, 
                   var.df  = 'stat.sim.Y1', 
                   obs.true = obs$A[ 1 ], 
                   best.sim = NULL
)


predictor  =  MaxWiK.predictor( psi = hyper$psi, t = hyper$t, param = par.sim, 
                                stat.sim = stat.sim, 
                                new.param = as.data.frame( t( obs$x0 ) ),
                                talkative = FALSE, check_pos_def = FALSE ,
                                n_bullets = 42, n_best = 12, halfwidth = 0.5, 
                                epsilon = 0.001, rate = 0.2, 
                                max_iteration = 10, save_web = TRUE, 
                                use.iKernelABC = res.MaxWik 
                                )



pred.network  = unique.data.frame( do.call(rbind.data.frame, predictor$spiderweb ) )

MaxWiK.ggplot.density( title = ' Posteriori distribution of Y1 parameter \n (Red line is a true value, blue one is prediction), \n red area is metasampling', 
                   datafr1 = posteriori.pred.MaxWiK, 
                   datafr2 = pred.network, 
                   var.df  = 'stat.sim.Y1', 
                   obs.true = obs$A[ 1 ], 
                   best.sim = as.numeric( predictor$prediction.best[ 1 ] )
)


MaxWiK.ggplot.density( title = ' Posteriori distribution of Y2 parameter at the last step of metasampling \n (Red line is a true value, blue one is prediction), \n red area is metasampling', 
                   datafr1 = posteriori.pred.MaxWiK, 
                   datafr2 = predictor$spiderweb[[predictor$iteration]], 
                   var.df  = 'stat.sim.Y2', 
                   obs.true = obs$A[ 2 ], 
                   best.sim = as.numeric( predictor$prediction.best[ 2 ] )
)




########################### SAMPLING FOR ABC --------------------------------------------------------

#' @describeIn Function to get output of Gaussian function for d dimension case
#'
#' @description The function \code{Gauss_function()} is used as 
#' a model based on Gaussian function for multidimensional case.
#' 
#' @return \code{Gauss_function()} returns output of Gaussian function for d dimensions 
#' 
#' @export
#'
#' @examples
#' NULL
Gauss_function  <-  function( d = 1, x0 = 3, r = range(0,10), noise = 0, 
                              A = 1, sigma = 1, x, nms = 'Y1' ){
    
    par.sim  =   as.numeric( x )
    sim1  =  data.frame( matrix( NA, nrow = 1, ncol = d ) )
    names( sim1 )  =  nms
    for( i in 1:d ){
        sim1[ 1, i ]  = A[ i ] * ( exp( x = - ( par.sim[ i ] - x0[ i ] ) ** 2 / 2 / sigma[ i ] / sigma[ i ] ) + 
                                       runif( n = 1, min = -0.5, max = 0.5 ) * noise )
    }
    
    return( sim1 )
}

model_par = list( d = obs$dimension[ nrow(obs) ], x0 = obs$x0, r = obs$range, 
                 A = obs$A, sigma = obs$sigma, 
                 noise = 0, nms = names( stat.sim)  ) 

model_function  =  Gauss_function

psi_t = data.frame( psi = hyper$psi, t = hyper$t) 

stat.obs = stat.sim[1,]
stat.obs[1,] = obs$A

# Check the call of the model:
do.call( what = model_function, args = c( model_par, list( x =  obs$x0 ) ) )


smpl_1  =  sampler_MaxWiK( stat.obs =  stat.obs, 
                           stat.sim =  stat.sim, 
                           par.sim  =  par.sim,  
                           model    =  model_function, 
                           arg0     =  model_par, 
                           size     =  1600, 
                           psi_t    =  psi_t, 
                           epsilon  =  1E-10, 
                           check_err  =  FALSE, 
                           nmax     =  60, 
                           include_top  =  TRUE,
                           slowly       =  TRUE, 
                           rate         =  0.05, 
                           n_simulation_stop = 1000,
                           include_web_rings = TRUE,
                           number_of_nodes_in_ring = 4  )


MSE  =  data.frame( x = 1:length( smpl_1$results$mse ), y = smpl_1$results$mse )
X12 = data.frame( i  = 1:length( smpl_1$results$par.sim.X1 ), 
                  x1 = smpl_1$results$par.sim.X1,  
                  x2 = smpl_1$results$par.sim.X2 )


library( ggplot2 )
th =     theme(
    plot.title = element_text(color="red", size=14, face="bold.italic"),
    axis.title.x = element_text(color="black", size=16, face="bold"),
    axis.title.y = element_text(color="black", size=16, face="bold"), 
    axis.text    = element_text(color="black", size=13             )
)

ggplot(data = MSE, aes( x, y ) ) + 
    geom_line( linewidth = 0.7 )  +  scale_y_log10() + 
    ylab("Mean Squared Error") + xlab("Number of simulations") + th

ggplot(data = X12, aes( i ) ) + 
    geom_line(aes(y = x1 ), linewidth = 0.5 ) +
    geom_line(aes(y = x2 ), linewidth = 0.5 ) + 
    geom_hline( aes(yintercept= obs$x0[1]),   
                color='red', linetype=2, linewidth=0.4) +
    geom_hline( aes(yintercept= obs$x0[2]),   
                color='red', linetype=2, linewidth=0.4) +
    ylab("Parameters X1 and X2") + xlab("Number of simulations") + th

    
