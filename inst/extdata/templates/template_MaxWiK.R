###############  LIBRARY and FUNCTION  ####################
library('MaxWiK')
library('ggplot2')

my.ggplot.density  <-  function( title = '', datafr1, datafr2, var.df, 
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


my.ggplot.density( title = ' Posteriori distribution of X2 parameter (Red line is a true value) ', 
                   datafr1 = posteriori.MaxWiK, 
                   datafr2 = NULL, 
                   var.df  = 'par.sim.X2', 
                   obs.true = obs$x0[ 2 ], 
                   best.sim = NULL
                )

my.ggplot.density( title = ' Posteriori distribution of X1 parameter (Red line is a true value) ', 
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

my.ggplot.density( title = ' Posteriori distribution of X1 parameter \n (Red line is a true value, dotted is the best simulation) \n red area is metasampling', 
                   datafr1 = posteriori.MaxWiK, 
                   datafr2 = network, 
                   var.df  = 'par.sim.X1', 
                   obs.true = obs$x0[ 1 ], 
                   best.sim = as.numeric( meta.sampling$par.best[ 1 ] )
        )


my.ggplot.density( title = ' Posteriori distribution of X2 parameter \n (Red line is a true value, dotted is the best simulation) \n red area is metasampling', 
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


my.ggplot.density( title = ' Posteriori distribution of Y1 parameter \n (Red line is a true value)', 
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

my.ggplot.density( title = ' Posteriori distribution of Y1 parameter \n (Red line is a true value, blue one is prediction), \n red area is metasampling', 
                   datafr1 = posteriori.pred.MaxWiK, 
                   datafr2 = pred.network, 
                   var.df  = 'stat.sim.Y1', 
                   obs.true = obs$A[ 1 ], 
                   best.sim = as.numeric( predictor$prediction.best[ 1 ] )
)


my.ggplot.density( title = ' Posteriori distribution of Y2 parameter at the last step of metasampling \n (Red line is a true value, blue one is prediction), \n red area is metasampling', 
                   datafr1 = posteriori.pred.MaxWiK, 
                   datafr2 = predictor$spiderweb[[predictor$iteration]], 
                   var.df  = 'stat.sim.Y2', 
                   obs.true = obs$A[ 2 ], 
                   best.sim = as.numeric( predictor$prediction.best[ 2 ] )
)


