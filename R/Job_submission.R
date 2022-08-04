### Job submission
if ( FALSE ){
    library('MaxWiK')
    DF_EXP  =  experiment_models( file_name = '../Results.txt', 
                                  models = c( 'Gaussian', 'Linear' ),
                                  dimensions = (1:20)*2, 
                                  stochastic_terms  =  c( 0, 0.1, 0.3, 0.7, 1, 1.5 ),
                                  rng  =  c( 0,10 ), 
                                  restrict_points_number = 1000 )
}