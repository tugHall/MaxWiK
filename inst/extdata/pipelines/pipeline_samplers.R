### 
###  Pipeline to get parameter estimation based on sampling schemes

library('MaxWiK')
library('parallel')

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
restrict_points_number  =  Number_of_points  =  500
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


Get_data  =  function( dimension, rng, 
                       Number_of_points, model,
                       A, sigma, stochastic_term, 
                       probability  =  TRUE ){
    
    input  =  NULL
    x0  =  runif( n = dimension, min = rng[1], max = rng[2] )
    Number_of_points  =  max( c( 50 * dimension, Number_of_points ) )
    
    if ( model == 'Gaussian' ) {
        input = get_dataset_of_Gaussian_model( d = dimension, x0 = x0, probability = probability, 
                                               n = Number_of_points, r = rng, 
                                               A = A[1:dimension], sigma = sigma[1:dimension], 
                                               noise = stochastic_term )
        
        model_par = list(d = dimension, x0 = x0, r = rng, 
                         A = A[1:dimension], sigma = sigma[1:dimension], 
                         noise = stochastic_term  ) 
    }
    if ( model == 'Linear' ) {
        input  =  get_dataset_of_Linear_model( d = dimension, x0 = x0, probability = probability, 
                                               n = Number_of_points, r = rng, A = A[1:dimension],
                                               noise = stochastic_term )
        model_par = list( d = dimension, x0 = x0, r = rng, 
                          A = A[1:dimension],
                          noise = stochastic_term  ) 
    }
    
    if ( is.null( input ) ) stop( 'Model name is incorrect' )
    stat.sim  =  input$stat.sim
    stat.obs  =  input$stat.obs
    par.sim   =  input$par.sim
    rm( input )
    
    return( list( stat.obs  =  stat.obs,
                  stat.sim  =  stat.sim, 
                  par.sim   =  par.sim,
                  model_par =  model_par ) )
}


input  =  Get_data( dimension = dimension, rng = rng, 
                    Number_of_points = Number_of_points,
                    model = 'Gaussian',
                    A = A, sigma = sigma, 
                    stochastic_term = stochastic_term, 
                    probability = TRUE )
par.sim   =  input$par.sim 
stat.sim  =  input$stat.sim
stat.obs  =  input$stat.obs
par.truth =  input$model_par$x0
model_par =  input$model_par

hyper    =  Get_hyperparameters(stat.obs = stat.obs, stat.sim = stat.sim, 
                              par.sim = par.sim, par.truth = par.truth )

utils::capture.output( hyper, file = 'HyperParameters.txt', append = FALSE )

psi_t  =  hyper$iKernel$psi_t

# Maximal number of iteration in sampling for each method
nmax  =  1000

if ( model  == 'Gaussian' ){
    
    input  =  get_dataset_of_Gaussian_model( d = dimension, 
                                             x0 = model_par$x0, 
                                             r = model_par$r, 
                                             noise = model_par$noise, 
                                             A = model_par$A, 
                                             sigma = model_par$sigma, 
                                             probability = FALSE, 
                                             n = nmax )
}

par.sim   =  rbind( par.sim,  input$par.sim )
stat.sim  =  rbind( stat.sim, input$stat.sim )
rm( input )


rej      =  abc( target = stat.obs, param = par.sim, 
                 sumstat = stat.sim, tol = hyper$tolerance, method = 'rejection' ) 

w        =  which( rej$region )
w        =  w[ which( w > ( nrow( par.sim ) - nmax ) ) ]

rm( rej )

if ( length( w ) == 0 ) stop( 'There is NO points closed to observation one!' )


data_par_est  =  list( init  =  par.sim[ 1 : ( nrow( par.sim) - nmax ), ], w = w )
data_MSE      =  data.frame( w = w )

############### Start of a method

Meth_Kern  =  data.frame( Method = c('K2-ABC', 'K2-ABC', 'K2-ABC', 'Rejection', 
                                     'Loclinear', 'Neuralnet', 'Ridge',
                                     'MaxWiK_MAP', 'MaxWiK' ), 
                          Kernel = c('Gaussian', 'Laplacian', 'iKernel', '',
                                     '',          '',          'epanechnikov',
                                     'iKernel', 'iKernel') 
)

# Start of looop for all the methods:
for( j in 1:length( Meth_Kern$Method ) ){
    
    method_name  =  Meth_Kern$Method[ j ]
    kernel_name  =  Meth_Kern$Kernel[ j ]
    
    MK  =  paste0( method_name, '_', kernel_name )
    
    print( MK )
    pb <- txtProgressBar(min = 0, max = max( w ), style = 3)
    # DF  =  NULL
    # time_start  =  Sys.time( )
    res  =  lapply( X = w, FUN = function( x ){
    
        G = NULL
        if ( kernel_name  ==  'iKernel' ){
            ikern  =  iKernelABC( psi = psi_t$psi[1], t = psi_t$t[1], 
                                  param = par.sim[1:x, ], 
                                  stat.sim = stat.sim[1:x, ], 
                                  stat.obs = stat.obs, 
                                  talkative = FALSE, 
                                  check_pos_def = FALSE )
            
            G = matrix( data = ikern$similarity, ncol = 1 )
        }
        new_par = Get_parameter( method_name = method_name, 
                                 kernel_name = kernel_name, 
                                 stat.obs   =  stat.obs, 
                                 stat.sim   =  stat.sim[1:x, ], 
                                 par.sim    =  par.sim[1:x, ], 
                                 G          =  NULL, 
                                 hyper      =  hyper )
        setTxtProgressBar(pb, x)
        return( new_par )
        
    } )
    # running_time  =  as.numeric( difftime(Sys.time(), time_start, units = "secs")[[1]] )
    
    data_par_est[[ MK ]]  =  as.data.frame( do.call( rbind, res ) ) 
    
    for( i in 1:nrow( data_par_est[[ MK ]] ) ){      
        
        new_par  =  data_par_est[[ MK ]][ i, ]
        model_par_all  =  c( model_par, list( par.sim1 = new_par ) )
        model_par_all$noise  =  0 
        sim_est  =  do.call( model_function, model_par_all )
        
        MSE  =  MSE_sim(stat.obs = stat.obs, stat.sim = sim_est ) 
    
        data_MSE[ i, MK ]  = MSE
    }
    
}   # End of loop for all the methods




plot(x = w, y = log( res_2 ), type = 'p' )
points( x = w, y = log( res_1 ), col = 'blue' )    
    






















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






