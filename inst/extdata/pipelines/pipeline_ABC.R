library('MaxWiK')

# Please, pay attention running of pipeline will take 4-6 hours for 4 cores. 

# Define the working folder:
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
###             *   Get_call_all_methods   - to call all the methods, and
###             *   Get_call               - to call particular method with 
###                                 fitting hyper parameters for each method

### This pipeline has two models with different parameters,
###      here you can change hyper parameters as well as the models' data like
###      dimension and stochastic term.

file_name   =  'output.txt'
models      =  c( 'Gaussian', 'Linear' )
dimensions  =  (1:20)*2
stochastic_terms   =   c( 0, 0.1, 0.5, 1, 2 )
rng  =  c( 0, 1000 )   # range of parameters
restrict_points_number  =  500
d = max( dimensions )
A      =  ( ( 1:d ) + 12 ) * 100  # Amplitude for Gauss function / Linear function
sigma  =  rep( rng[2]/ 5, d )     # Sigma for Gauss function 

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


input  =  Get_data( dimension = dimensions[5], rng = rng, 
                    restrict_points_number = restrict_points_number, 
                    Number_of_points = Number_of_points,
                    model = 'Gaussian',
                    A = A, sigma = sigma, stochastic_term = 0 )
par.sim   =  input$par.sim 
stat.sim  =  input$stat.sim
stat.obs  =  input$stat.obs
par.truth =  input$model_par$x0

hyper  =  Get_hyperparameters(stat.obs = stat.obs, stat.sim = stat.sim, 
                              par.sim = par.sim, par.truth = par.truth )
utils::capture.output( hyper, file = 'HyperParameters.txt', append = FALSE )

psi_t  =  hyper$iKernel$psi_t

DF = NULL   # Data frame to collect results of all the simulations
for( model in models ){
    
    if ( model == 'Gaussian' ) {
        model_function  =  Gauss_function
    }
    if ( model == 'Linear' ) {
        model_function  =  Linear_function
    }
    
    for( dimension in dimensions ){
        for( stochastic_term in stochastic_terms ){
            
            input  =  Get_data( dimension = dimension, rng = rng, 
                                restrict_points_number = restrict_points_number, 
                                Number_of_points = Number_of_points, 
                                model = model,
                                A = A, sigma = sigma, stochastic_term = stochastic_term )
            par.sim   =  input$par.sim 
            stat.sim  =  input$stat.sim
            stat.obs  =  input$stat.obs
            model_par =  input$model_par
            
            DF_new  =  Get_call_all_methods(    
                dimension  = dimension,
                iterations  =  1:12,
                stat.obs = stat.obs, 
                stat.sim = stat.sim, 
                par.sim  = par.sim, 
                par.truth  =  model_par$x0, 
                cores = cores,
                model_function = model_function, 
                model_par = model_par, 
                hyper  =  hyper )
            DF_new$model_name  =  model
            
            DF  =  rbind( DF, DF_new )
            
            if ( file.exists( file_name ) ){
                write.table(file = file_name, x = DF_new , append = TRUE, sep = '\t', 
                            row.names = FALSE, col.names = FALSE )
            } else {
                write.table(file = file_name, x = DF_new , append = TRUE, sep = '\t', 
                            row.names = FALSE, col.names = TRUE )
            }
        }
    }
}

###  Results of simulation is collected in dataset DF
print( 'Results of simulation is collected in dataset DF' )



