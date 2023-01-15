### 
###  Pipeline to get parameter estimation based on sampling schemes

library('MaxWiK')

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

### In order to understand procedure of execution all the method, please, see
### the functions:
###             *   sampler_all_methods    - to call all the methods, and
###             *   sampler_method         - to call particular method with 
###                                 fitting hyper parameters for each method

### This pipeline is the content of the experiment_samplers() function,
###      here you can change hyper parameters as well as the models' data like
###      name, dimension and stochastic term.


file_name   =  'output.txt'
model       =  c( 'Gaussian', 'Linear' )[ 1 ]
dimension   =  2
stochastic_term   =   c( 0, 0.1, 0.5, 1, 2 )[ 1 ]
rng  =  c( 0, 1000 )   # range of parameters
restrict_points_number  =  500
d = max( dimension )
A      =  ( ( 1:d ) + 12 ) * 100  # Amplitude for Gauss function / Linear function
sigma  =  rep( rng[2]/ 5, d )     # Sigma for Gauss function 

if ( model == 'Gaussian' ) {
    model_function  =  Gauss_function
}
if ( model == 'Linear' ) {
    model_function  =  Linear_function
}

#  By default number of processors in parallel calculations
#             cores = 4 in the function Get_call_all_methods
#   YOU can change this parameter to accelerate calculations ( see the code below )
cores = 4

# delete old file
if ( file.exists( file_name) ) unlink( file_name )


Get_data  =  function( dimension, rng, restrict_points_number, 
                       Number_of_points, model,
                       A, sigma, stochastic_term ){
    
    input  =  NULL
    x0  =  runif( n = dimension, min = rng[1], max = rng[2] )
    Number_of_points  =  max( c( 50 * dimension, restrict_points_number ) )
    
    if ( model == 'Gaussian' ) {
        input = get_dataset_of_Gaussian_model( d = dimension, x0 = x0, probability = TRUE, 
                                               n = Number_of_points, r = rng, 
                                               A = A[1:dimension], sigma = sigma[1:dimension], 
                                               noise = stochastic_term )
        
        model_par = list(d = dimension, x0 = x0, r = rng, 
                         A = A[1:dimension], sigma = sigma[1:dimension], 
                         noise = stochastic_term  ) 
    }
    if ( model == 'Linear' ) {
        input  =  get_dataset_of_Linear_model( d = dimension, x0 = x0, probability = TRUE, 
                                               n = Number_of_points, r = rng, A = A[1:dimension],
                                               noise = stochastic_term )
        model_par = list( d = dimension, x0 = x0, r = rng, 
                          A = A[1:dimension],
                          noise = stochastic_term  ) 
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
    par.sim   =   par.sim_origin[ rej$region, ] 
    
    return( list( stat.obs  =  stat.obs,
                  stat.sim  =  stat.sim, 
                  par.sim   =  par.sim,
                  model_par =  model_par ) )
}


input  =  Get_data( dimension = dimension, rng = rng, 
                    restrict_points_number = restrict_points_number, 
                    Number_of_points = Number_of_points,
                    model = 'Gaussian',
                    A = A, sigma = sigma, stochastic_term = 0 )
par.sim   =  input$par.sim 
stat.sim  =  input$stat.sim
stat.obs  =  input$stat.obs
par.truth =  input$model_par$x0
model_par =  input$model_par

stat.sim_init =  stat.sim
par.sim_init  =  par.sim

hyper  =  Get_hyperparameters(stat.obs = stat.obs, stat.sim = stat.sim, 
                              par.sim = par.sim, par.truth = par.truth )
utils::capture.output( hyper, file = 'HyperParameters.txt', append = FALSE )

psi_t  =  hyper$iKernel$psi_t


Meth_Kern  =  data.frame( Method = c('K2-ABC', 'K2-ABC', 'K2-ABC', 'Rejection', 
                                     'Loclinear', 'Neuralnet', 'Ridge',
                                     'MaxWiK_MAP', 'MaxWiK' ), 
                          Kernel = c('Gaussian', 'Laplacian', 'iKernel', '',
                                     '',          '',          'epanechnikov',
                                     'iKernel', 'iKernel') 
)


method_name  =  Meth_Kern$Method[ 4 ]
kernel_name  =  Meth_Kern$Kernel[ 4 ]

stat.sim_itt =  stat.sim
par.sim_itt  =  par.sim

# Maximal number of itteration in sampling for each method
nmax  =  100   

res = NULL
for( itt in 1:nmax ){
    
    input  =  restrict_data( par.sim = par.sim_itt, 
                             stat.sim = stat.sim_itt, 
                             stat.obs = stat.obs, 
                             size = restrict_points_number ) 
    stat.sim_itt =  input$stat.sim
    par.sim_itt  =  input$par.sim
    
    new_par = Get_parameter( method_name = method_name, 
                             kernel_name = kernel_name, 
                             stat.obs   =  stat.obs, 
                             stat.sim   =  stat.sim_itt, 
                             par.sim    =  par.sim_itt, 
                             G          =  NULL, 
                             hyper      =  hyper )
    
    res = rbind( res, new_par )
    
    add_arg  =  list( par.sim1 =  as.data.frame( new_par ) )

    new_sim  =  do.call( what = model_function, 
                         args = c( model_par, add_arg ) )
    
    par.sim_itt  =  rbind( par.sim_itt, new_par )
    stat.sim_itt =  rbind( stat.sim_itt, new_sim )
}




















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






