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

Get_MSE  <-  function( new_par, model_par, model_function, stat.obs ){
    
    model_par_all  =  c( model_par, list( par.sim1 = new_par ) )
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
               names = c( 'Iterations', 'log of MSE'), xr = c(500, 1500), 
               yr = c(1E-5, 1E6), logscale = 'y', 
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
    
    new_par  =  data.frame( matrix( x, nrow = 1 ) )

    model_par_all  =  c( model_par, list( par.sim1 = new_par ) )
    
    model_par_all$noise  =  0 
    
    y     =  as.numeric( do.call( model_function, model_par_all ) )
    
    return( y )
}

### Results:
MSE_samplings  =  list(  )

toy_prior  =  list( c( "unif", model_par$r ), c( "unif", model_par$r ) )
sum_stat_obs  = as.numeric( stat.obs )  #  c( 100, 100 )
set.seed(1)


############# REJECTION ABC
n=1000
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

tolerance  =  ( 20 : 1 )*0.02  # c( 4E-1, 1E-1, 4E-2 )
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

tolerance  =  0.5  # c( 4E-1, 1E-1, 4E-2 )
alpha_delmo  =  0.5
n = 200
cntr = 0
ABC_Delmoral  =  ABC_sequential( method = "Delmoral",
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
     y = MSE_samplings$Delmoral$MSE, log = 'y', type = 'p' , 
     xlim = c( 0, 1000 ))
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
      log  =  'y', ylim = c(0.1, 1E6) )






# MaxWiK sampling ---------------------------------------------------------

smpl_1  =  sampler_MaxWiK( stat.obs =  stat.obs, 
                           stat.sim =  stat.sim, 
                           par.sim  =  par.sim,  
                           model    =  model_function, 
                           arg0     =  model_par, 
                           size     =  2000, 
                           psi_t    =  psi_t, 
                           epsilon  =  1E-12, 
                           nmax     =  30, 
                           include_top  =  TRUE,
                           slowly       =  TRUE, 
                           rate         =  0.2  )



# Sequential ABC - does not work ------------------------------------------



### Performing a A Simulated Annealing Approach to Approximate Bayes Computations scheme

### Ref:
# Albert C., Kunsch HR., Scheidegger A. (2014)
# A Simulated Annealing Approach to Approximate Bayes Computations.
# Stat. Comput., 1–16, arXiv:1208.2157.

# Sampler:
r.prior  =  function()   c( runif( 1, 1, 1000 ), runif( 1, 1, 1000 ) )

# Density:
d.prior  =  function(x)  dunif( x[1], 1, 1000 ) * dunif( x[2], 1, 1000 )


n.sample  =  300

iter.max  =  n.sample * 20

eps.init  =  2
cntr  =  0
ABC_Albert  =  SABC(   r.model  =  toy_model,
                       r.prior  =  r.prior,
                       d.prior  =  d.prior,
                       n.sample =  n.sample,
                       eps.init =  eps.init,
                       iter.max =  iter.max,
                       method   =  "informative", 
                       resample = 100, 
                       y        =  sum_stat_obs )

print( paste0( 'The number of simulations is ', cntr ) )

hist( ABC_Albert$E[ , 1 ], breaks = 25 )
hist( ABC_Albert$E[ , 2 ], breaks = 25 )


#   











# OLD SECTION -------------------------------------------------------------




### OLD CODE

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






