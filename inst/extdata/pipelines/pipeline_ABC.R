library('MaxWiK')

# Please, pay attention running of pipeline will take 4-6 hours for 4 cores. 

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
    DF  =  experiment_models( file_name = '../Results_ALL.txt', 
                              models = c( 'Gaussian', 'Linear' ),
                              dimensions = (1:10)*2, 
                              stochastic_terms  =  c( 0, 1, 5, 10, 20, 30 ),
                              rng  =  c( 0,10 ), 
                              restrict_points_number = 300 )
}

### In order to understand procedure of execution all the method, please, see
### the functions:
###             *   Get_call_all_methods   - to call all the methods, and
###             *   Get_call               - to call particular method with 
###                                 fitting hyper parameters for each method

### This pipeline is the content of the experiment_models() function,
###      here you can change hyper parameters as well as the models' data like
###      dimension and stochastic term.

file_name   =  'output.txt'
models      =  c( 'Gaussian', 'Linear' )
dimensions  =  (1:20)*2
stochastic_terms   =   c( 0, 0.1, 0.3, 0.7, 1, 1.5 )
rng  =  c( 0, 10 )   # range of parameters
restrict_points_number  =  300

#  By default number of processors in parallel calculations
#             cores = 4 in the function Get_call_all_methods
#   YOU can change this parameter to accelerate calculations ( see the code below )
cores = 4

# delete old file
if ( file.exists( file_name) ) unlink( file_name )


x0     =  c(10, 50, 90, 130, 180, 280, 390, 430, 520, 630, 1010, 1050, 1090, 1130, 1180, 1280, 1390, 1430, 1520, 1630)
sigma  =  x0 / 5
x = 1:20000 / 10 
A = 250
f <- function( A, x, x0, sigma ){
    Gaus  =  function( A, x, x0, sigma ) A * exp( -( x - x0 ) ** 2 / sigma )
    y = sum( unlist( lapply( 1:length( x0 ), FUN = function( i ) Gaus( A = A, x = x,
                                        x0 = x0[ i ], sigma  =  sigma[ i ]  )  )  )  )
    return( y )
}


y = sapply( 1:length( x ), FUN = function( i ) f( A = A, x = x[ i ], x0 = x0, sigma = sigma ) )

plot( x, y, type = 'l', xlim = c( 0, 200 ) )

#### ___________________________________________________________________

DF = NULL   # Data frame to collect results of all the simulations
for( model in models ){
    for( dimension in dimensions ){
        for( stochastic_term in stochastic_terms ){
            
            input  =  NULL
            x0  =  runif( n = dimension, min = rng[1], max = rng[2] )
            Number_of_points  =  max( c( 50 * dimension, restrict_points_number ) )
            
            if ( model == 'Gaussian' ) {
                input = get_dataset_of_Gaussian_model( d = dimension, x0 = x0, probability = TRUE, 
                                        n = Number_of_points, r = rng,
                                        noise = stochastic_term )
            }
            if ( model == 'Linear' ) {
                input  =  get_dataset_of_Linear_model( d = dimension, x0 = x0, probability = TRUE, 
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
            par.sim   =   par.sim_origin[ rej$region, ] 
            
            psi_t  =  adjust_psi_t( par.sim = par.sim, stat.sim = stat.sim, 
                                    stat.obs = stat.obs, talkative = FALSE, 
                                    check_pos_def = FALSE, 
                                    n_best = 8, cores = 4 )
            
            ikern  =  iKernelABC( psi = psi_t$psi[1], t = psi_t$t[1], 
                                  param = par.sim, 
                                  stat.sim = stat.sim, 
                                  stat.obs = stat.obs, 
                                  talkative = FALSE, 
                                  check_pos_def = FALSE )
            
            G = matrix( data = ikern$similarity, ncol = 1 )
            DF_new  =  Get_call_all_methods(    
                model_name = model, 
                dimension  = dimension,
                stochastic_term = stochastic_term, 
                iterations  =  1:12,
                stat.obs = stat.obs, 
                stat.sim = stat.sim, 
                par.sim  = par.sim, 
                G        = G, 
                par.truth  =  x0, 
                cores = cores )
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



