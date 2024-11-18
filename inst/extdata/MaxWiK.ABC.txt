
library('MaxWiK')
library('ggplot2')

# Load the data:

hyper     =  MaxWiK::Data.2D$ABC$hyperparameters
par.sim   =  MaxWiK::Data.2D$X
stat.sim  =  MaxWiK::Data.2D$Y
obs       =  MaxWiK::Data.2D$observation
Matrix.Voronoi  =  MaxWiK::Data.2D$ABC$Matrix.Voronoi

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

# Get the Posteriori where similarity more than 0
w.sim  =  which(res.MaxWik$similarity > 0 )
posteriori.MaxWiK  =  par.sim[ w.sim, ]
posteriori.weights =  res.MaxWik$similarity[ w.sim ] / sum( res.MaxWik$similarity[ w.sim ] )


MaxWiK.ggplot.density( title = ' Posteriori distribution of X1 parameter (Red line is a true value) ', 
                   datafr1 = posteriori.MaxWiK, 
                   datafr2 = NULL, 
                   var.df  = 'par.sim.X1', 
                   obs.true = obs$x0[ 1 ], 
                   best.sim = NULL
)

MaxWiK.ggplot.density( title = ' Posteriori distribution of X2 parameter (Red line is a true value) ', 
                   datafr1 = posteriori.MaxWiK, 
                   datafr2 = NULL, 
                   var.df  = 'par.sim.X2', 
                   obs.true = obs$x0[ 2 ], 
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


