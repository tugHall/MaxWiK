
library('MaxWiK')

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

library('ggplot2')
ggplot(  ) + 
    geom_density( data = posteriori.MaxWiK, 
                  aes( x =  par.sim.X2, weight = posteriori.weights ), 
                  fill="#69b3a2", color="#69b3a2", alpha=0.4 ) + 
    geom_density( data = posteriori.MaxWiK, 
                  aes( x =  par.sim.X2 ), 
                  fill="#f9b3a2", color="#f9b3a2", alpha=0.4 ) + 
    geom_vline( aes(xintercept= obs$x0[ 2 ]),   color="red", linetype="dashed", linewidth=0.7) +
    ggtitle( ' Posteriori distribution of X2 parameter (Red line is a line of interest) ' ) 


meta.sampling  =  meta_sampling( psi = hyper$psi, t = hyper$t, param = par.sim, 
                              stat.sim = stat.sim, stat.obs = as.data.frame( t( obs$A ) ), 
                              talkative = FALSE, check_pos_def = FALSE ,
                              n_bullets = 42, n_best = 12, halfwidth = 0.5, 
                              epsilon = 0.001, rate = 0.2, 
                              max_iteration = 10, save_web = TRUE, 
                              use.iKernelABC = res.MaxWik )


network  = unique.data.frame( do.call(rbind.data.frame, meta.sampling$spiderweb ) )

ggplot(  ) + 
    geom_density( data = posteriori.MaxWiK, 
                  aes( x =  par.sim.X1 ), 
                  fill="#f9b3a2", color="#f9b3a2", alpha=0.4 ) + 
    geom_vline( aes(xintercept= obs$x0[ 1 ]),   color="red", linetype="dashed", linewidth=0.7) + 
    geom_density( data = network, 
                  aes( x =  par.sim.X1 ), 
                  fill="#69b3a2", color="#69b3a2", alpha=0.4 ) +
    ggtitle( ' Posteriori distribution of X1 parameter (Red line is a line of interest) ' ) 

ggplot(  ) + 
    geom_density( data = posteriori.MaxWiK, 
                  aes( x =  par.sim.X2 ), 
                  fill="#f9b3a2", color="#f9b3a2", alpha=0.4 ) + 
    geom_vline( aes(xintercept= obs$x0[ 2 ]),   color="red", linetype="dashed", linewidth=0.7) + 
    geom_density( data = network, 
                  aes( x =  par.sim.X2 ), 
                  fill="#69b3a2", color="#69b3a2", alpha=0.4 ) +
    ggtitle( ' Posteriori distribution of X2 parameter (Red line is a line of interest) ' ) 
