
library('MaxWiK')
library('ggplot2')

# Load the data:

hyper     =  MaxWiK::Data.2D$ABC$hyperparameters
par.sim   =  MaxWiK::Data.2D$X
stat.sim  =  MaxWiK::Data.2D$Y
obs       =  MaxWiK::Data.2D$observation
Matrix.Voronoi  =  MaxWiK::Data.2D$ABC$Matrix.Voronoi



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
                                use.iKernelABC = pred.MaxWiK 
)



pred.network  = unique.data.frame( do.call(rbind.data.frame, predictor$spiderweb ) )

# Check the diapason
pred.network  = apply_range( diapason = c(0,1000), input.data = pred.network )
predictor$prediction.best = apply_range( diapason = c(0,1000), input.data = predictor$prediction.best )

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



