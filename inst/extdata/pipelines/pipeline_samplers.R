### 
###  Pipeline to get parameter estimation based on sampling schemes


# Input data --------------------------------------------------------------


library('MaxWiK')
library('randomcoloR')
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

file_name   =  'output.txt'
model       =  c( 'Gaussian', 'Linear' )[ 1 ]
dimension   =  20
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

Get_MSE  <-  function( new_par, model_par, model_function, stat.obs ){
    
    model_par_all  =  c( model_par, list( x = new_par ) )
    model_par_all$noise  =  0 
    sim_est  =  do.call( model_function, model_par_all )
    
    MSE  =  MSE_sim(stat.obs = stat.obs, stat.sim = sim_est ) 
    
    return( MSE )
}

input  =  Get_data( dimension = dimension, rng = rng, 
                    Number_of_points = Number_of_points,
                    model = 'Gaussian',
                    A = A, sigma = sigma, 
                    stochastic_term = stochastic_term, 
                    probability = FALSE )
par.sim   =  input$par.sim 
stat.sim  =  input$stat.sim
stat.obs  =  input$stat.obs
par.truth =  input$model_par$x0
model_par =  input$model_par

hyper    =  Get_hyperparameters(stat.obs = stat.obs, stat.sim = stat.sim, 
                              par.sim = par.sim, par.truth = par.truth )

utils::capture.output( hyper, file = 'HyperParameters.txt', append = FALSE )

psi_t  =  hyper$iKernel$psi_t



# Input extra data for sampling based on uniform prior --------------------



# Maximal number of iteration in sampling for each method
nmax  =  4500

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



# Sampling on kernels -----------------------------------------------------


############### Start of a method

Meth_Kern  =  data.frame( Method = c('K2-ABC', 'K2-ABC', 'K2-ABC', 
                                     'Rejection', 'Loclinear',
                                     'MaxWiK_MAP', 'MaxWiK' ), 
                          Kernel = c('Gaussian', 'Laplacian', 'iKernel',
                                     '', '',
                                     'iKernel', 'iKernel') 
)

# Start of loop for all the methods:
for( j in 1:length( Meth_Kern$Method ) ){
    
    method_name  =  Meth_Kern$Method[ j ]
    kernel_name  =  Meth_Kern$Kernel[ j ]
    
    MK  =  paste0( method_name, '_', kernel_name )
    
    print( MK )
    pb <- txtProgressBar(min = 0, max = max( w ), style = 3)
    # DF  =  NULL
    # time_start  =  Sys.time( )
    res  =  lapply( X = w, FUN = function( x ){
    
        # G = NULL
        if ( FALSE & method_name == 'K2-ABC' & kernel_name == 'iKernel' ){
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
                                 hyper      =  hyper )
        setTxtProgressBar(pb, x)
        return( new_par )
        
    } )
    # running_time  =  as.numeric( difftime(Sys.time(), time_start, units = "secs")[[1]] )
    
    data_par_est[[ MK ]]  =  as.data.frame( do.call( rbind, res ) ) 
    
    for( i in 1:nrow( data_par_est[[ MK ]] ) ){      
        
        new_par  =  data_par_est[[ MK ]][ i, ]
        
        MSE  =  Get_MSE( new_par    =  new_par, 
                         model_par  =  model_par, 
                         model_function = model_function, 
                         stat.obs = stat.obs )
    
        data_MSE[ i, MK ]  = MSE
    }
    
}   # End of loop for all the methods

nl  =  c(2:8 ) # 8:10)
l   =  length( nl )

hue = c(" ", "random", "red", "orange", "yellow",
        "green", "blue", "purple", "pink", "monochrome")[1]
luminosity = c(" ", "random", "light", "bright", "dark")[5] 

clrs  =  randomColor(count = l,
                     hue = hue,
                     luminosity = luminosity )

# palette.colors( n = l, palette = 'ggplot2' )  # as.vector( gen_colors( nm = l ) )

plot_2D_lines( x = data_MSE$w, DF = data_MSE, nl = nl, 
               names = c( 'Iterations', 'log of MSE'), 
               xr = c(1000, 5500), 
               yr = c(5E6, 6E7), 
               logscale = '', 
               col = clrs, 
               lwd = 2, lt = 1:l, cex = 1.5, 
               draw_key = TRUE )

# Save data to a file
data_par_est$MSE  =  data_MSE
saveRDS( object = data_par_est, 
         file = paste0('par_est_', dimension, '_', stochastic_term, '.RDS' ))




# Sampling from EasyABC package and methods from it ------------------------------------

library( 'EasyABC' )

toy_model  =  function( x ){
    
    cntr  <<-  cntr + 1
    
    new_par  =  data.frame( matrix( as.numeric(x), nrow = 1 ) )

    model_par_all  =  c( model_par, list( x = new_par ) )
    
    model_par_all$noise  =  0 
    
    y     =  as.numeric( do.call( model_function, model_par_all ) )
    
    return( y )
}

### Results:
MSE_samplings  =  list(  )

# toy_prior  =  list( c( "unif", model_par$r ), c( "unif", model_par$r ) )
el  =  c( "unif", model_par$r )
toy_prior  =  list( el )[ rep( 1, dimension ) ]
sum_stat_obs  = as.numeric( stat.obs )  #  c( 100, 100 )
set.seed(1234)


############# REJECTION ABC
n=2500
cntr = 0
ABC_rej  =  ABC_rejection( model = toy_model, prior = toy_prior,
                           nb_simul = n,
                           summary_stat_target = sum_stat_obs,
                           tol = 0.02,
                           progress_bar = TRUE )
print( paste0( 'The number of simulations is ', ABC_rej$nsim ) )

ABC_rej$param
ABC_rej$stats

hist( ABC_rej$param[ , 1] )
hist( ABC_rej$param[ , 2] )

MaxWiK:: Get_MAP( DF = as.data.frame( ABC_rej$param ) )

############# Adaptive ABC or sequential ABC scheme

### Ref:
# Beaumont, M. A., Cornuet, J., Marin, J., and Robert, C. P. (2009)
# Adaptive approximate Bayesian computation. Biometrika, 96, 983–990.

tolerance  =  ( 20 : 18 )*0.03  # c( 4E-1, 1E-1, 4E-2 )
n = 100
cntr = 0
ABC_Beaumont  =  ABC_sequential( method = "Beaumont",
                                 model  = toy_model,
                                 prior  = toy_prior,
                                 nb_simul = n,
                                 summary_stat_target = sum_stat_obs,
                                 tolerance_tab = tolerance,
                                 verbose = TRUE )

print( paste0( 'The number of simulations is ', cntr ) )

ABC_Beaumont$weights
ABC_Beaumont$param

hist( ABC_Beaumont$param[, 1])
hist( ABC_Beaumont$param[, 2])

# Weighted estimations:
sum( ABC_Beaumont$weights * ABC_Beaumont$param[ , 1 ] )
sum( ABC_Beaumont$weights * ABC_Beaumont$param[ , 2 ] )

MaxWiK:: Get_MAP( DF = as.data.frame( ABC_Beaumont$param ) )


MSE_samplings$ABC_Beaumont  =  data.frame()
for (i in 1 : length( ABC_Beaumont$intermediary ) ){
    
    new_par  =  sapply( X = 2 : ( 1 + ncol( par.sim ) ), FUN = function( x ) {
                        sum( ABC_Beaumont$intermediary[[ i ]]$posterior[ , x ] * 
                        ABC_Beaumont$intermediary[[ i ]]$posterior[ , 1 ] ) 
                    } )
    MSE_samplings$ABC_Beaumont[ i, 'MSE' ]  =  
        Get_MSE(new_par = new_par, model_par = model_par, 
                    model_function = model_function, stat.obs = stat.obs )
    MSE_samplings$ABC_Beaumont[ i, 'n_simul_tot' ]  =  
        ABC_Beaumont$intermediary[[ i ]]$n_simul_tot
}

plot(x = MSE_samplings$ABC_Beaumont$n_simul_tot, 
     y = MSE_samplings$ABC_Beaumont$MSE, log = 'y', type = 'l' )

##### Sequential Monte-Carlo
### See:
### Del Moral, P., Doucet, A., and Jasra, A. (2012) 
###     An adaptive sequential Monte Carlo method for approximate Bayesian computation. 
###     Statistics and Computing, 22, 1009–1020.

tolerance  =  0.1  # c( 4E-1, 1E-1, 4E-2 )
alpha_delmo  =  0.5
n = 300
cntr = 0
ABC_Delmoral_2  =  ABC_sequential( method = "Delmoral",
                                 model  = toy_model,
                                 prior  = toy_prior,
                                 nb_simul = n,
                                 summary_stat_target = sum_stat_obs,
                                 tolerance_target = tolerance, 
                                 alpha  =  alpha_delmo, 
                                 verbose = TRUE )

print( paste0( 'The number of simulations is ', cntr ) )

MSE_samplings$Delmoral  =  data.frame( n_simul_tot = NA, MSE = NA )
for (i in 1 : length( ABC_Delmoral$intermediary ) ){
    
    new_par  =  sapply( X = 2 : ( 1 + ncol( par.sim ) ), FUN = function( x ) {
        sum( ABC_Delmoral$intermediary[[ i ]]$posterior[ , x ] * 
                 ABC_Delmoral$intermediary[[ i ]]$posterior[ , 1 ] ) 
    } )
    MSE_samplings$Delmoral[ i, 'MSE' ]  =  
        Get_MSE(new_par = new_par, model_par = model_par, 
                model_function = model_function, stat.obs = stat.obs )
    MSE_samplings$Delmoral[ i, 'n_simul_tot' ]  =  
        ABC_Delmoral$intermediary[[ i ]]$n_simul_tot
}

plot(x = MSE_samplings$Delmoral$n_simul_tot, 
     y = MSE_samplings$Delmoral$MSE, log = '', type = 'l' , 
     xlim = c( 0, 7000 ))
points(x = data_MSE$w, data_MSE$`K2-ABC_Laplacian`, pch = 16 )






#### Performing a ABC-MCMC scheme

# n  =  20
n_between_sampling  =  10
MSE_samplings$ABC_Marjoram_original  =  data.frame( n = NA, MSE = NA )
for( i in 1:30 ){

    n  =  20 * i
    ABC_Marjoram_original  =  ABC_mcmc( method="Marjoram_original",
                                    model=toy_model,
                                    prior=toy_prior,
                                    summary_stat_target=sum_stat_obs,
                                    n_rec=n, 
                                    n_between_sampling = n_between_sampling )

    # print( paste0( 'The number of simulations is ', ABC_Marjoram_original$nsim ) )

    # ABC_Marjoram_original$param

    # hist( ABC_Marjoram_original$param[ , 1 ] )
    # hist( ABC_Marjoram_original$param[ , 2 ] )

    # ABC_Marjoram_original$stats_normalization

    new_par  =  MaxWiK:: Get_MAP( DF = as.data.frame( ABC_Marjoram_original$param ) )
    
    MSE  =  Get_MSE( new_par    =  new_par, 
                         model_par  =  model_par, 
                         model_function  =  model_function, 
                         stat.obs = stat.obs )
    MSE_samplings$ABC_Marjoram_original[ i, 'n']    =  ABC_Marjoram_original$nsim
    MSE_samplings$ABC_Marjoram_original[ i, 'MSE']  =  MSE
}

plot( x = MSE_samplings$ABC_Marjoram_original$n, 
      y = MSE_samplings$ABC_Marjoram_original$MSE, 
      log  =  '', ylim = c(1E4, 3E7) )






# MaxWiK sampling ---------------------------------------------------------

# Restrict number of initial simulations
stat.sim  =  stat.sim[ 1:1000, ]
par.sim   =  par.sim[ 1:1000, ]

smpl_1  =  sampler_MaxWiK( stat.obs =  stat.obs, 
                           stat.sim =  stat.sim, 
                           par.sim  =  par.sim,  
                           model    =  model_function, 
                           arg0     =  model_par, 
                           size     =  1000, 
                           psi_t    =  psi_t, 
                           epsilon  =  1E-8, 
                           nmax     =  30, 
                           include_top  =  TRUE,
                           slowly       =  TRUE, 
                           rate         =  0.2, 
                           n_simulation_stop = 3000  )
# Get correct MSE with noise = 0
smpl_1$results$mse  =  sapply(  X = 1:nrow(smpl_1$results), 
                                FUN = function( x ) Get_MSE(new_par = smpl_1$results[ x, 1:dimension ], 
                                                            model_par = model_par, 
                                                            model_function = model_function, 
                                                            stat.obs = stat.obs ) )


MSE_samplings$kernels  =  data_MSE 
MSE_samplings$MaxWiK   =  data.frame( n_sim_total  =  ( nrow(stat.sim) + 1 ) : ( nrow(stat.sim) + nrow( smpl_1$results ) ), 
                                      MSE          =  smpl_1$results$mse )

plot( MSE_samplings$MaxWiK$n_sim_total, 
      MSE_samplings$MaxWiK$MSE, type = 'l', log = 'y' )


### Save all the data:

saveRDS( object = MSE_samplings, 
         file = paste0('RESULTS_dim_', dimension, '_noise_', stochastic_term, '.RDS' ))

