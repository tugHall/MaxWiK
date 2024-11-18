
library('MaxWiK')
library('ggplot2')

# Load the data:

hyper     =  MaxWiK::Data.2D$ABC$hyperparameters
par.sim   =  MaxWiK::Data.2D$X
stat.sim  =  MaxWiK::Data.2D$Y
obs       =  MaxWiK::Data.2D$observation
Matrix.Voronoi  =  MaxWiK::Data.2D$ABC$Matrix.Voronoi


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
                           n_simulation_stop = 400,
                           include_web_rings = F,
                           number_of_nodes_in_ring = 1  )


MSE  =  data.frame( x = smpl_1$results$sim_ID, y = smpl_1$results$mse )
X12 = data.frame( i  = smpl_1$results$sim_ID, 
                  x1 = smpl_1$results$par.sim.X1,  
                  x2 = smpl_1$results$par.sim.X2 )


th =     theme(
    plot.title = element_text(color="red", size=14, face="bold.italic"),
    axis.title.x = element_text(color="black", size=16, face="bold"),
    axis.title.y = element_text(color="black", size=16, face="bold"), 
    axis.text    = element_text(color="black", size=13             )
)

ggplot(data = MSE, aes( x, y ) ) + 
    geom_line( linewidth = 0.7 )  +  scale_y_log10() + 
    geom_smooth(method = "lm", alpha = .5, formula= y~x ) + 
    ylab("Mean Squared Error") + xlab("Number of simulations") + th

ggplot(data = X12, aes( i ) ) + 
    geom_line(aes(y = x1 ), linewidth = 0.5 ) +
    geom_line(aes(y = x2 ), linewidth = 0.5 ) + 
    geom_hline( aes(yintercept= obs$x0[1]),   
                color='red', linetype=2, linewidth=0.4) +
    geom_hline( aes(yintercept= obs$x0[2]),   
                color='red', linetype=2, linewidth=0.4) +
    ylab("Parameters X1 and X2") + xlab("Number of simulations") + th

cat('The best value after sampling procedure is \n') ; print( smpl_1$best)
cat('The observation or true values of X are \n ') ; dg = obs$x0; names(dg) = names(smpl_1$best[1:2]) ; print( dg )
cat('The observation or true values of Y are \n ') ; dg = obs$A;  names(dg) = names(smpl_1$best[3:4]) ; print( dg )
