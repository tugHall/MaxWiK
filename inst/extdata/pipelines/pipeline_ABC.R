

# Please, pay attention running of pipeline will take 4-6 hours for 4 cores. 

library('MaxWiK')

check_packages()

DF  =  experiment_models( file_name = '../Results_ALL.txt', 
                          models = c( 'Gaussian', 'Linear' ),
                          dimensions = (1:10)*2, 
                          stochastic_terms  =  c( 0, 1, 5, 10, 20, 30 ),
                          rng  =  c( 0,10 ), 
                          restrict_points_number = 300 )