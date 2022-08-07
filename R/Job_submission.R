### Job submission
if( TRUE ){
    library('MaxWiK')
    DF_lin  =  experiment_models( file_name = '../Results_lin.txt', 
                                      models = c( 'Linear' ),
                                      dimensions = (1:20)*2, 
                                      stochastic_terms  =  c( 0, 0.1, 0.3, 0.5, 1, 2 ),
                                      rng  =  c( 0,10 ), 
                                      restrict_points_number = 1000 )
}    
    