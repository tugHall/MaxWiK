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
    experiment_samplers(   file_name = './output.txt', 
                           model_name = 'Gaussian',
                           dimension = 6, 
                           stochastic_term  =  5,
                           rng  =  c( 0,10 ), 
                           restrict_points_number = 500, 
                           nmax = 300 )
}


### In order to understand procedure of execution all the method, please, see
### the functions:
###             *   Get_call_all_methods   - to call all the methods, and
###             *   Get_call               - to call particular method with 
###                                 fitting hyper parameters for each method

### This pipeline is the content of the experiment_models() function,
###      here you can change hyper parameters as well as the models' data like
###      dimension and stochastic term.

# We use these input parameters for experiment_samplers() function:

file_name = './output.txt'
model_name = 'Gaussian'
dimension = 2
stochastic_term  =  1
rng  =  c( 0,10 )
restrict_points_number = 500
nmax = 300 





