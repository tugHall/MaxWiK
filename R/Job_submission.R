### Job submission
if( FALSE ){
    library('MaxWiK')
    DF_Gauss  =  experiment_models( file_name = '../Results_Gauss.txt', 
                                      models = c( 'Gaussian' ),
                                      dimensions = (1:20)*2, 
                                      stochastic_terms  =  c( 0.5, 1, 5 ),
                                      rng  =  c( 0,10 ), 
                                      restrict_points_number = 1000 )
}    





