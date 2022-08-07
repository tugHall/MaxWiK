### Job submission
if( FALSE ){
    library('MaxWiK')
    DF_lin  =  experiment_models( file_name = '../Results_lin.txt', 
                                      models = c( 'Linear' ),
                                      dimensions = (1:20)*2, 
                                      stochastic_terms  =  c( 2, 3, 4, 5, 7, 10 ),
                                      rng  =  c( 0,10 ), 
                                      restrict_points_number = 1000 )
}    
    