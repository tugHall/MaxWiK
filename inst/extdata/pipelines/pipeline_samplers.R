### 
###  Pipeline to get parameter estimation 

library('MaxWiK')

# Please, pay attention running of pipeline will take X-X hours for 4 cores. 

# DEfine the working folder:
wd  =  '../Simulation'

# And create it:
if ( dir.exists( wd ) ) {
    setwd( wd )
} else {
    dir.create( wd )
    setwd( wd )
}

### Check installation of libraries:
check_packages()


### This is the function to execute the pipeline automatically, 
###    so, if you need just to repeat the results of the dataset, 
###    please, use this short cut 

if ( FALSE ){
    setwd( dir = wd )
    RES  =  experiment_samplers(   file_name = './output.txt', 
                                   model_name = 'Gaussian',
                                   dimension = 2, 
                                   stochastic_term  =  2,
                                   rng  =  c( 0,10 ), 
                                   restrict_points_number = 100, 
                                   nmax = 30 )
}


### In order to understand procedure of execution all the method, please, see
### the functions:
###             *   sampler_all_methods    - to call all the methods, and
###             *   sampler_method         - to call particular method with 
###                                 fitting hyper parameters for each method

### This pipeline is the content of the experiment_samplers() function,
###      here you can change hyper parameters as well as the models' data like
###      name, dimension and stochastic term.

# We use these input parameters for experiment_samplers() function:

file_name = './output.txt'
model_name = 'Gaussian'
dimension = 2
stochastic_term  =  1
rng  =  c( 0,10 )
restrict_points_number = 100
nmax = 30 


# delete the file from the previous simulation
if ( file.exists( file_name) ) unlink( file_name )

input  =  NULL
x0  =  round( runif( n = dimension, min = rng[1], max = rng[2] ), digits = 4 )
Number_of_points  =  max( c( 50 * dimension, restrict_points_number ) )

if ( model_name == 'Gaussian' ) {
    input = Gaussian_model( d = dimension, x0 = x0, probability = TRUE, 
                            n = Number_of_points, r = rng,
                            noise = stochastic_term )
}
if ( model_name == 'Linear' ) {
    input  =  linear_model( d = dimension, x0 = x0, probability = TRUE, 
                            n = Number_of_points, r = rng,
                            noise = stochastic_term )
}

if ( is.null( input ) ) stop( 'Model name is incorrect' )
stat.sim_origin  =  input$stat.sim
stat.obs  =  input$stat.obs
par.sim_origin  =  input$par.sim
rm( input )

# Apply restrict number of points:
tol = restrict_points_number / nrow( stat.sim_origin )
rej = abc::abc( target = stat.obs, param = par.sim_origin, sumstat = stat.sim_origin,
                method = 'rejection', tol = tol )

stat.sim  =  stat.sim_origin[ rej$region, ]
par.sim   =  par.sim_origin[ rej$region, ] 
par.truth =  data.frame( matrix( x0, ncol = dimension ) )
rm(x0)
psi_t  =  adjust_psi_t( par.sim = par.sim, stat.sim = stat.sim, 
                        stat.obs = stat.obs, talkative = FALSE, 
                        check_pos_def = FALSE, 
                        n_best = 10, cores = 4 )

ikern  =  iKernelABC( psi = psi_t$psi[1], t = psi_t$t[1], 
                      param = par.sim, 
                      stat.sim = stat.sim, 
                      stat.obs = stat.obs, 
                      talkative = FALSE, 
                      check_pos_def = FALSE )

G = matrix( data = ikern$similarity, ncol = 1 )

### Next block of the code is the content of the sampler_all_methods() function:
if ( FALSE ){
    RES = sampler_all_methods(  model_name = model_name, 
                                dimension = dimension, 
                                stochastic_term = stochastic_term, 
                                stat.obs = stat.obs, 
                                stat.sim = stat.sim, 
                                par.sim = par.sim, 
                                par.truth = par.truth, 
                                G = G,
                                nmax = nmax
                            )
}

cores = 4

DF  =  NULL
Meth_Kern  =  data.frame( Method = c('K2-ABC', 'K2-ABC', 'K2-ABC', 'Rejection', 
                                     'Loclinear', 'Neuralnet', 'Ridge',
                                     'MaxWiK_MAP', 'MaxWiK' ), 
                          Kernel = c('Gaussian', 'Laplacian', 'iKernel', '',
                                     '',          '',          'epanechnikov',
                                     'iKernel', 'iKernel') 
)

# Define model function:
if ( model_name == 'Gaussian' ) {
    model  =  model
    arg0 = list(  name = c( 'Gaussian', 'Linear' )[1],
                  x0 = par.truth, 
                  stat.obs = stat.obs, 
                  noise = stochastic_term )
} else {
    model  =  model
    arg0 = list(  name = c( 'Gaussian', 'Linear' )[2],
                  x0 = par.truth, 
                  stat.obs = stat.obs, 
                  noise = stochastic_term )
} 

### Plot for 2D example:
plot_2D( x = par.sim$x1, y = par.sim$x2, xr = c(0,10), yr = c( 0,10) )
points( x = par.truth$X1, y = par.truth$X2, col = 'red', pch = 18, cex = 2.5 )



DF_1_S  =  mclapply( 1:nrow(Meth_Kern) , FUN = function( mk ){   
    sampler_method( method_name =  as.character( Meth_Kern$Method[mk] ), 
                    kernel_name =  as.character( Meth_Kern$Kernel[mk] ), 
                    model_name  =  model_name, 
                    dimension   =  dimension, 
                    stochastic_term  =  stochastic_term, 
                    nmax = nmax,
                    stat.obs   =  stat.obs, 
                    stat.sim   =  stat.sim, 
                    par.sim    =  par.sim, 
                    G          =  G, 
                    par.truth  =  par.truth, 
                    model      =  model, 
                    arg0       =  arg0,  
                    size       =  ifelse( Number_of_points > 500, 500 , Number_of_points )  
    )
}, mc.cores  =  cores )




# Check an error
bad   =  sapply( DF_1_S, inherits, what = "try-error" )
# If NO errors:
if ( any( !bad ) ){
    for( mk in ( 1:nrow(Meth_Kern) )[ !bad ] ){
        DF_1  =  DF_1_S[[ mk ]]$results
        DF    =  rbind( DF, DF_1 )
    }
    # DF_2  =  do.call( rbind, DF_1_S[ !bad ]$results )
    # DF    =  rbind( DF, DF_2 )
}
# If error in some core(s):
if ( any( bad ) ){
    its = which( bad )
    
    DF_3  =   data.frame(   method_name = as.character( Meth_Kern$Method[ its ] ),
                            kernel_name = as.character( Meth_Kern$Kernel[ its ] ),
                            MSE = NA, 
                            running_time = NA ) 
    # Add circle for its:
    for( i in its ){
        DF_3[ 1, 'iteration' ]  =  i
        DF    =  rbind( DF, DF_3 )
    }
}






