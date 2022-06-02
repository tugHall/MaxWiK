
# Initialization of data and norm definition ------------------------------


### The kernel functions

### The norm function for vector
norm_vec     <-  function(x) sqrt(sum(x^2))
norm_vec_sq  <-  function(x)     (sum(x^2))

### The function is to check DATA.FRAME 
### that it has numeric format for ALL columns and it has NO 'NA' values
check_numeric_format  <-  function( l ) {

    if ( all(!is.na( l ) )  & all( sapply(l, is.numeric) ) & is.data.frame( l ) ) return( TRUE )
    msg  =  NULL
    if ( !all(!is.na( l )) )     msg = paste0( msg, ' - input data has NA value(s);\n')
    if ( !is.data.frame( l ) )   msg = paste0( msg, ' - input data should be TYPE of data frame;\n')
    if ( !all( sapply(l, is.numeric) ) )  msg = paste0( msg, ' - input data should has ONLY NUMERIC TYPE for ALL columns;\n')
    msg    =    paste0( substr(msg,1,nchar(msg)-2), '.') 
    stop( msg )
}


# Gaussian and Laplacian kernel -------------------------------------------
### 1. Radial basis function (RDF) or Gaussian kernel, and Laplacian kernel
###     https://en.wikipedia.org/wiki/Kernel_smoother 


### Let define parameters of kernels such that k(x,y) = 1 if x = y

Gaussian_kernel  <-  function( x1, x2 , sgm = 1)  exp( -1 * norm_vec_sq(x1 - x2) / (2 * sgm**2 ) )

Laplacian_kernel  <-  function( x1, x2 , b = 1 )  exp( -1 * norm_vec(x1 - x2)    / (2 * b ) )


# Isolation kernel --------------------------------------------------------
### ISOLATION KERNEL BASED ON VORONOI DIAGRAM 

### The common coefficients:
### t is a number of trees in iKernel or dimension of RKHS
### psi is a size of each Voronoi diagram
### nr is a size of dataset (number of rows in Dataset)

### The functions

### The function to get subset with size psi for Voronoi diagram 
GET_SUBSET  <-  function( data_set, pnts ){
    voronoi_subset  <-  data_set[  , pnts]
    return( voronoi_subset )
}

### The function to calculate a similarity between two points
### in the Matrix of real numbers.

### It returns the value of similarity or Isolation KERNEL for TWO points  
iKernel  <-  function( Matrix_iKernel, pnt_1, pnt_2, t ){
    ### Input data:
    ### t is a number of columns of Matrix_iKernel
    ### pnt_1 and pnt_2 are IDs of points 1 and 2 in the Matrix_iKernel
    ### 
    
    # t  <-  ncol( Matrix_iKernel )
    smlrt  <- sum( Matrix_iKernel[ pnt_1, ]  == Matrix_iKernel[ pnt_2, ] ) / t
    return( smlrt )
}

### The function to get Isolation Kernel between a new point and dataset
iKernel_point_dataset  <-  function( Matrix_iKernel, t, nr, iFeature_point ){
    ### Input data:
    ### t  is a number of columns of Matrix_iKernel
    ### nr is a number of rows of Matrix_iKernel
    ### iFeature_point is feature map in RKHS for a new point
    ### iFeature_point is a result of the add_new_point_iKernel()  function
    ### 
    
    smlrt_p_data  <- sapply( X = 1:nr, FUN = function( x ) sum( iFeature_point  == Matrix_iKernel[ x, ] ) / t  )
    
    return( smlrt_p_data )
}

### The function to get feature representation in RKHS based on Voronoi diagram for WHOLE dataset
get_voronoi_feature  <-  function( psi = 40, t = 350, data, talkative = FALSE, 
                                   new = TRUE,  Matrix_Voronoi = NULL ){
    ### Input parameters
    ### psi is number of points in a subset / dimension of RKHS 
    ### psi is also number of areas in any Voronoi diagram
    ### t is number of trees or number of Voronoi diagrams
    ### new is the logical parameter - Is Matrix_Voronoi new ?
    ### if new = FALSE then you should define Matrix_Voronoi
    
    ### Check the data format:
    if (talkative ) print( 'Check the data format')
    if ( talkative & check_numeric_format( data ) ) print( 'OK' )
    if ( new & talkative ) print( 'This is a new Voronoi diagram and Matrix') 
    if ( !new ) {
        if (talkative ) print( 'This is a calculation based on a GIVEN Voronoi diagram and Matrix')
        psi = ncol( Matrix_Voronoi )
        t   = nrow( Matrix_Voronoi )
    }
    
    if (talkative ) print( 'Finding matrix of distances between all points in data')
    dissim  <-  dist( data, method = 'euclidean', diag = TRUE, upper = TRUE)
    dissim  <-  as.matrix( dissim )
    nr  =  nrow(dissim)  ## number of columns and rows in matrix of distances
    if (talkative ) print( 'Done' )
    
    if (talkative ) print( 'Transform data to the Hilbert space related to Isolation Kernel' )
    ### The matrix 'Matrix_iKernel' keeps IDs of Voronoi area for each point and each tree
    Matrix_iKernel  =  matrix( data = NA, nrow = nr, ncol = t )
    
    ### The matrix 'Matrix_Voronoi' keeps the subsets for each tree (psi points IDs for each tree)
    if ( new ) Matrix_Voronoi  =  matrix( data = NA, nrow = t, ncol = psi )
    
    ### Get the psi points for EACH Voronoi diagram:
    for ( j in 1:t ) {
        if ( new ) {
            pnts  =  sort( sample( 1:nr, psi ) )
            Matrix_Voronoi[ j, ]  =  pnts 
        } else {
            pnts  =  Matrix_Voronoi[ j, ]
        }
        sub_data  =  dissim[, pnts] 
        
        for (i in 1:nr) {
            Matrix_iKernel[i, j]   =  which.min( sub_data[ i, ] )[1]  ### which( sub_data[ i, ] == min( sub_data[ i, ] ) )[1]
        }
    }
    if (talkative ) print( 'Done' )
    
    return( list( M_Voronoi = Matrix_Voronoi, 
                  M_iKernel = Matrix_iKernel,
                  M_dissim  = dissim ) )
}


### The function to get feature representation in RKHS based on Voronoi diagram for PART of dataset
### The NEW PART of data set show be at the end of PREVIOUS dataset
### The Matrix_Voronoi is based on the PREVIOUS dataset
get_voronoi_feature_PART_dataset  <-  function( data, talkative = FALSE, start_row,  Matrix_Voronoi ){
    ### Input parameters
    ### psi is a number of points in a subset / dimension of RKHS 
    ### psi is also number of areas in any Voronoi diagram
    ### t is number of trees or number of Voronoi diagrams
    ### start_row is a row number from which a new data have added
    ### Matrix_Voronoi is a matrix of Voronoi diagrams based on the PREVIOUS dataset
    
    ### Check the data format:
    if (talkative ) print( 'Check the data format')
    if ( talkative & check_numeric_format( data ) ) print( 'OK' )
    if (talkative ) print( 'This is a calculation based on a GIVEN Voronoi diagram and Matrix')
    psi = ncol( Matrix_Voronoi )
    t   = nrow( Matrix_Voronoi )
    
    if (talkative ) print( 'Finding matrix of distances between all points in data')
    dissim  <-  as.matrix( dist( data, method = 'euclidean', diag = TRUE, upper = TRUE) )
    nnw  =  nrow(dissim) - start_row  + 1  ### number of new points
    nr   =  nrow(dissim)                   ### number of rows 
    if (talkative ) print( 'Done' )
    
    if (talkative ) print( 'Transform a NEW data to the Hilbert space related to Isolation Kernel' )
    ### The matrix 'Matrix_iKernel' keeps IDs of Voronoi area for each point and each tree
    Matrix_iKernel  =  matrix( data = NA, nrow = nnw, ncol = t )
    
    ### The matrix 'Matrix_Voronoi' keeps the subsets for each tree (psi points IDs for each tree)
    
    ### Get the feature mapping for NEW points based on GIVEN Voronoi diagram:
    for ( j in 1:t ) {
        pnts  =  Matrix_Voronoi[ j, ]
        sub_data  =  dissim[ start_row:nr, pnts ]
        
        for (i in 1:nnw) {
            Matrix_iKernel[i, j]   =  which.min( sub_data[ i, ] )[1]  ### which( sub_data[ i, ] == min( sub_data[ i, ] ) )[1]
        }
    }
    if (talkative ) print( 'Done' )
    
    return( list( M_Voronoi = Matrix_Voronoi, 
                  M_iKernel = Matrix_iKernel,
                  M_dissim  = dissim ) )
}

### The function to get RKHS mapping for a new point (Isolation Kernel)
add_new_point_iKernel  <-  function( data, d1, Matrix_Voronoi, dissim, t, psi, nr ){
    ### Input data:
    ### d1 is a data point - usually it is an observation data point
    ### dissim is a matrix of dissimilarity or distances between all points
    ### nr is a number of rows in dissim matrix (matrix of distances)
    ### nr is also size of data
    
    dissim_new  <-  matrix( 0,  nr +1, nr +1 )
    dissim_new[1:nr, 1:nr]  <-  dissim 
    
    ### Get distances between new point and points from the dataset
    dlt  =  data
    dlt    <-    sapply( X = 1:ncol(dlt), FUN = function( x )  dlt[ , x ] - d1[ 1, x ] )
    dissim_new[ nr+1,1:nr ]   <-  sapply( X = 1:nr, FUN = function(x) norm_vec( dlt[ x, ] ) )
    dissim_new[ 1:nr, nr+1 ]  <-  dissim_new[ nr+1,1:nr ]
    
    
    ### Get the feature map of a new point for EACH Voronoi diagram:
    iFeature_point  <-  rep( 0, t)
    
    for ( j in 1:t ) {
        pnts      <-  Matrix_Voronoi[ j, ] 
        sub_data  <-  dissim_new[ , pnts] 
        iFeature_point[ j ]  <-  which.min( sub_data[ nr + 1, ] )[1]
    }
    
    if ( min(iFeature_point ) == 0) stop( 'Error in calculation of feature map for a new point.' )
    return( iFeature_point )
}

### The function to get weights from Feature mapping
get_weights_iKernel   <-   function( GI, Matrix_iKernel, t, nr, iFeature_point ){
    
    smlrt_p_data    =  iKernel_point_dataset( Matrix_iKernel, t, nr, iFeature_point )
    
    PHI_T_k_S       =  matrix( data = smlrt_p_data, ncol = 1 ) 
    
    WGHTS  =  GI  %*%  PHI_T_k_S
    
    return( list( weights_RKHS = WGHTS, weights_similarity = PHI_T_k_S ) )
}


### The function to calculate kernel mean embedding for Isolation Kernel of parameters
get_kernel_mean_embedding  <-  function( parameters_Matrix_iKernel, Hilbert_weights ){
    ### Input:
    ### parameters_Matrix_iKernel is a matrix of of all points of PARAMETERS in a Hilbert space
    ### (rows - points, columns - isolation trees)
    ### 
    ### Hilbert_weights is a weights in Hilbert space to get kernel mean embedding for parameters_Matrix_iKernel
    w  = Hilbert_weights[,1]                 #  weights
    t  =  ncol( parameters_Matrix_iKernel )  #  number of trees in the Isolation Forest
    mu = rep( 0, t)
    for( i in 1:t ){
        ### iteration for trees:
        vct   =  parameters_Matrix_iKernel[ , i ]
        pnts  =  sort( unique( vct ) )
        sms   =  rep( 0, length( pnts ) )
        ### iteration for points 
        for( j in 1:length( pnts ) ){
            ids       =  which( vct == pnts[ j ] )
            sms[ j ]  =  sum( w[ ids ] )   
        }
        mu[ i ]  =  pnts[ which.max( sms ) ]
    }
    
    return( mu )
}


### The function to get subset of points based on feature map:
get_subset_of_feature_map   <-   function(dtst, Matrix_Voronoi, iFeature_point ){
    ### Input
    ### dtst is a dataset of original points (all points)
    ### Matrix_Voronoi is a matrix of Voronoi diagram
    ### iFeature_point is a feature mapping of a point in the Hilbert space
    
    rws  =  NULL
    
    for( i in 1:length( iFeature_point ) ){
        id_map  =  iFeature_point[ i ] 
        rws  =  c( rws, Matrix_Voronoi[i, id_map ] )
    }
    rws  =  sort( unique( rws ) )
    if ( is.null( rws) ) stop( ' subset is NULL. ' )
    return( dtst[ rws, ] )
}




### Gram matrix -------------------------------------------------------------
### The function to calculate Gram matrix for Isolation Kernel method
GRAM_iKernel    <-  function( Matrix_iKernel, check_pos_def = FALSE ){
    t   <-  ncol( Matrix_iKernel )
    nr  <-  nrow( Matrix_iKernel )
    G   <-  matrix( data = 0, ncol = nr, nrow = nr ) 
    for (i in 1:nr){  
        for (j in 1:i) {
            G[i,j]  <-  sum( Matrix_iKernel[ i, ]  ==  Matrix_iKernel[ j, ] ) / t
            G[j,i]  <-  G[i,j]  ###  symmetric matrix
        }
    }
    
    if (check_pos_def ) {
        if ( !check_positive_definite( G, n = 50 ) ) stop( 'The Gram matrix is NOT positive definite !')
    }
    
    return( G )
}

### The function to get inverse Gram matrix
get_inverse_GRAM  <-  function( G, l = 1E-6, check_pos_def = FALSE ){
    ### l is a regularization constant
    n = nrow( G )
    I = diag( n ) * n * l
    G = G - I
    d = det( G )
    # if ( d == 0 & !check_pos_def ) print( 'The Gram matrix has determinant = 0 or close to zero' )    # stop( 'The Gram matrix has determinant = 0 \n' )
    
    GI  =  solve( G ) 
    if ( check_pos_def ) {
        if ( !check_positive_definite( GI, n = 50 ) ) stop( 'The Gram matrix is NOT positive definite !')
    }
    
    
    return( GI )
}


### The function to check the positive definite property of Gram matrix:
check_positive_definite  <-  function( G , n = 10 ){
    ### G - Gram matrix, n - number of iterations to check the positive definite property
    nr = nrow( G )
    check_pf  =  TRUE
    
    for( i in 1:n ){
        rand_vect  =  as.matrix( runif( nr, -1, 0 )  )
        x = t( rand_vect)  %*%  G  %*%  rand_vect
        if ( as.numeric( x ) < 0 )  check_pf = FALSE
        # print( paste( as.numeric( x ) , check_pf, sep = ' - ' ) )
    }
    
    return( check_pf )
}







# Kernel ABC --------------------------------------------------------------


iKernelABC  <-  function( psi = 40, t = 350, param, 
                          stat.sim, stat.obs, talkative = FALSE, check_pos_def = TRUE ){
    ### Input parameters
    ### psi is number of points in a subsets 
    ### t is number of trees or Voronoi diagrams 
    
    ### param is a data frame of parameters of the model
    ### stat.sim is a summary statistics of the simulations (model output)
    ### stat.obs is a summary statistics of the observation data
    par.sim  =  param
    ### Check the data format:
    if (talkative ) print( 'Check the data format of the summary statistics of observation data')
    if ( talkative & check_numeric_format(stat.obs) ) print( 'OK' )
    if (talkative ) print( 'Check the data format of the summary statistics of simulation data')
    if ( talkative & check_numeric_format(stat.sim) ) print( 'OK' )
    
    Voron  =  get_voronoi_feature( psi = psi, t = t, data = stat.sim, 
                                   talkative = talkative, new = TRUE, Matrix_Voronoi = NULL )
    
    Matrix_iKernel = Voron[['M_iKernel']]
    Matrix_Voronoi = Voron[['M_Voronoi']]
    dissim         = Voron[['M_dissim']]
    
    nr = nrow(stat.sim)
    
    ### The function to get RKHS mapping for a new point (Isolation Kernel)
    ### The new point is observation point of the summary statistics
    iFeature_point  =  add_new_point_iKernel( data = stat.sim, d1 = stat.obs, 
                                              Matrix_Voronoi, dissim, t, psi, nr )
    
    if (talkative ) print( 'Get the Gram matrix and inverse Gram matrix' )
    ### Get Gram matrix
    G = GRAM_iKernel( Matrix_iKernel, check_pos_def = check_pos_def )
    
    ### Get inverse Gram matrix, l - regularization coefficient
    GI  <-  get_inverse_GRAM( G , l = 1E-5 )
    if (talkative ) print( 'OK') 
    
    if (talkative ) print( 'Get the weights for Isolation Kernel' )
    weights_iKernel  =  get_weights_iKernel( GI, Matrix_iKernel, t, nr, iFeature_point )
    if (talkative ) print( 'OK') 
    
    if (talkative ) print( 'Transform the parameters data set to Hilbert space' )
    param_Voron  =  get_voronoi_feature( psi = psi, t = t, data = par.sim, talkative = talkative, 
                                         new = FALSE, Matrix_Voronoi = Matrix_Voronoi )
    parameters_Matrix_iKernel = param_Voron[['M_iKernel']]
    parameters_Matrix_Voronoi = param_Voron[['M_Voronoi']]
    
    kernel_mean_embedding  =  get_kernel_mean_embedding( parameters_Matrix_iKernel = parameters_Matrix_iKernel, 
                                                         Hilbert_weights = weights_iKernel[['weights_RKHS']] )  # t( parameters_Matrix_iKernel )  %*%  weights_iKernel[['weights_RKHS']]
    if (talkative ) print( 'Finish') 
    
    description = "The result of Kernel Approximate Bayesian Computation
    based on Isolation Kernel:
    - kernel_mean_embedding is a kernel mean embedding of observation point;
    - parameters_Matrix_Voronoi is a matrix of information about Voronoi trees (rows - trees, columns - Voronoi points/areas IDs) for parameters data set;
    - parameters_Matrix_iKernel is a matrix of of all points of PARAMETERS in a Hilbert space (rows - points, columns - isolation trees);
    - Hilbert_weights is a weights in Hilbert space to get kernel mean embedding for parameters_Matrix_iKernel;
    - Matrix_iKernel is a matrix of all points of simulations in a Hilbert space (rows - points, columns - isolation trees);
    - iFeature_point is a feature embedding map for OBSERVATION point;
    - similarity is a vector of similarity between the simulation points and observation point;
    - Matrix_Voronoi is a matrix of information about Voronoi trees (rows - trees, columns - Voronoi points/areas IDs);
    - t is a number of trees in the Isolation Forest; 
    - psi is a number of areas/points in the Voronoi diagrams."
    
    if (talkative ) message( description )
    
    return( list( kernel_mean_embedding = kernel_mean_embedding, 
                  parameters_Matrix_Voronoi  =  parameters_Matrix_Voronoi,
                  parameters_Matrix_iKernel  =  parameters_Matrix_iKernel, 
                  Hilbert_weights  =  weights_iKernel[['weights_RKHS']], 
                  similarity = weights_iKernel[['weights_similarity']], 
                  iFeature_point  =  iFeature_point, 
                  Matrix_iKernel = Matrix_iKernel, 
                  Matrix_Voronoi = Matrix_Voronoi, 
                  description    = description, 
                  t = t, 
                  psi = psi )  )
}


# SUDOKU ------------------------------------------------------------------

### This is a heuristic algorithm to seek a space/area related to 
### the feature mapping in Hilbert space for the dataset of the parameters

### Idea is just:
###     1) Generate points between the centers of Voronoi diagrams related to the feature mapping
###     2) Following strategy to puzzle out of SUDOKU: delete all points 
###         that do not match feature map. 
###     Output: The remaining points should be related to feature mapping.  

### Get pairs from Data Frame
### Output is a list of pairs
get_pairs_of_data_frame  <-  function( DF ){
    ### Input
    ### DF is data frame to get pairs
    
    ### get distances between all points
    DS  =  as.matrix( dist( DF, diag = FALSE, upper = FALSE ) )
    cl  =  max.col( DS )
    
    GEN_DF  =  vector( mode = "list", length = length( cl ) )
    for( i in 1:length( cl ) ){
        GEN_DF[[ i ]]  =  rbind( DF[ i, ], DF[ cl[i], ] )
    }
    
    return( GEN_DF )
}

### The function to generate points between pair of given points
### Output is data frame of points
generate_points_between_two_points  <-  function( pair, n = 10 ){
    ### Input
    ### pair is a data.frame of two points
    DF = pair
    DF[ 1:n,] = 0
    row.names( DF ) = 1:n
    for( i in 1:length( DF ) ){
        dlt = ( pair[ 2, i ] - pair[ 1, i ] ) / ( n + 1 )
        DF[ , i ]  =  pair[ 1, i ] + dlt * ( 1:n )
    }
    
    return( DF )
}

### The function to get 'tracer bullets' - generated points between all the farthest distance points
get_tracer_bullets  <-  function( DF , n_bullets = 20 ){
    ### Input:
    ### DF is a data frame
    
    ### get list of pairs of the farthest distance points in data frame
    pairs_list  =  get_pairs_of_data_frame( DF )
    
    tracer_bullets  =  NULL
    
    for( i in 1:length( pairs_list ) ){
        points_between_two_points  =  
            generate_points_between_two_points( pair = pairs_list[[ i ]], n = n_bullets )
        
        tracer_bullets  =  rbind( tracer_bullets,  points_between_two_points)
        
    }
    
    return( tracer_bullets )
}

### The function to get the best tracer bullets related to kernel mean embedding.
### The calculation performs ONLY for parameters dataset DT = par.sim 
sudoku  <-  function( DT , iKernelABC, n_bullets = 20, n_best = 10, halfwidth = 0.5 ){
    ### Input:
    ### DT is a data set of parameters
    ### iKernelABC is a result of calculations based on Isolation Kernel ABC
    ### n_bullets is a number of tracer bullets / additional points between the TWO farthest distance points
    ### n_best is a number of the best tracer bullets / points
    
    if ( FALSE){
        sbst_feature_Y  =  get_subset_of_feature_map( dtst  =  stat.sim, 
                                                  Matrix_Voronoi = iKernelABC$Matrix_Voronoi, 
                                                  iFeature_point = iKernelABC$iFeature_point )
    }
    
    sbst_feature_Param  =  get_subset_of_feature_map( dtst  =  DT, 
                                                      Matrix_Voronoi = iKernelABC$parameters_Matrix_Voronoi, 
                                                      iFeature_point = iKernelABC$kernel_mean_embedding )
    
    tracer_bullets   =  get_tracer_bullets( DF = sbst_feature_Param, n_bullets = n_bullets )
    
    feature_tracers  =  get_voronoi_feature_PART_dataset( data = rbind( DT, tracer_bullets ), 
                                        talkative = FALSE, start_row = nrow( DT ) + 1 ,  
                                        Matrix_Voronoi = iKernelABC$parameters_Matrix_Voronoi )
    
    ### Get similarity vector between kernel mean embedding and tracer bullets / new points
    sim_tracer_bullets  =  iKernel_point_dataset( Matrix_iKernel = feature_tracers$M_iKernel, 
                                                  t = iKernelABC$t, nr = nrow( feature_tracers$M_iKernel ), 
                                                  iFeature_point = iKernelABC$kernel_mean_embedding )
    
    ### Let take the best points:
    criterion  =  sort(sim_tracer_bullets, decreasing = TRUE )[ n_best ]
    best_tracer_bullets  =  tracer_bullets[ which(sim_tracer_bullets >= criterion), ]
    cr_out  =  criterion
    
    # criterion  =  sort(sim_tracer_bullets, decreasing = TRUE )[ n_best * 4 ]
    surroundings_best_points  =  tracer_bullets[ which(sim_tracer_bullets >= halfwidth ), ]
    
    return( list( tracer_bullets  =  tracer_bullets, 
                  criterion  =  cr_out,  
                  best_tracer_bullets  =  best_tracer_bullets,
                  surroundings_best_points  =  surroundings_best_points, 
                  feature_tracers  =  feature_tracers, 
                  similarity_to_mean  =  sim_tracer_bullets ) )
}







# SPIDERWEB algorithm -----------------------------------------------------

spiderweb  <-  function( psi = 4, t = 35, param = par.sim, 
                         stat.sim = stat.sim, stat.obs = stat.obs, 
                         talkative = FALSE, check_pos_def = FALSE ,
                         n_bullets = 10, n_best = 20, halfwidth = 0.5, 
                         epsilon = 0.001 ){
    
    input.parameters  =  list( psi = psi, t = t, param = param, 
                               stat.sim = stat.sim, stat.obs = stat.obs, 
                               talkative = talkative, check_pos_def = check_pos_def,
                               n_bullets = n_bullets, n_best = n_best, 
                               halfwidth = halfwidth, epsilon = epsilon )
    
    iKernelABC  = iKernelABC( psi = psi, t = t, param = param, 
                              stat.sim = stat.sim, stat.obs = stat.obs, 
                              talkative = talkative, check_pos_def = check_pos_def )
    
    ### The function to apply SUDOKU algorithm to get the best tracer bullets
    rslt  =  sudoku( DT = param , iKernelABC = iKernelABC, 
                      n_bullets = n_bullets, n_best = n_best, halfwidth = halfwidth )
    
    ### Get the top:
    tracers = rslt$tracer_bullets
    
    tracers_all  =  tracers
    sim.tracers_all  =  rslt$similarity_to_mean
    
    par.top = tracers[ order( rslt$similarity_to_mean , decreasing = TRUE)[1:n_best], ]
    par.best  =  par.top[ 1, ]
    par.top   = par.top[2:n_best, ]
    rm( tracers )
    sim_previous  =  0  # max( rslt$similarity_to_mean )
    
    sim.top   =  NULL
    sim.best  =  -1
    while( TRUE ){
        ### Reflect par.top through par.best 
        par.reflect  =  par.top
        for( i in 1:nrow( par.top ) )  par.reflect[ i, ]  =  2 * par.best - par.reflect[ i ,  ] 
        
        ### Generate points between par.top and par.reflect:
        tracers  =  rbind( par.best, par.top, par.reflect )
        for( i in 1:nrow( par.top  ) ){
            gen_tr  =  generate_points_between_two_points( pair = rbind( par.top[ i ,] , 
                                                    par.reflect[ i, ] ), n = n_bullets )
            tracers  =  rbind( tracers, gen_tr )
        }
        
        ### calculate the similarity for new points:
        feature_tracers  =  get_voronoi_feature_PART_dataset( data = rbind( param, tracers ), 
                                                              talkative = talkative, start_row = nrow( param ) + 1 ,  
                                                              Matrix_Voronoi = iKernelABC$parameters_Matrix_Voronoi )
        
        sim_tracers  =  iKernel_point_dataset( Matrix_iKernel = feature_tracers$M_iKernel, 
                                               t = iKernelABC$t, nr = nrow( feature_tracers$M_iKernel ), 
                                               iFeature_point = iKernelABC$kernel_mean_embedding )
        
        tracers_all  =  rbind( tracers_all, tracers )
        sim.tracers_all  =  c( sim.tracers_all, sim_tracers )
        # new best point
        par.best     =  tracers[ which.max(sim_tracers ), ]
        
        if ( abs( max( sim_tracers ) - sim_previous ) < epsilon ) break
        sim_previous   =   max( sim_tracers )
        
        ### rename new tracers:
        sim.top  =  sort( sim_tracers , decreasing = TRUE)[1:n_best]
        sim.best =  sim.top[ 1 ]
        sim.top  =  sim.top[2:n_best]
        
        par.top  =  tracers[ order( sim_tracers , decreasing = TRUE)[1:n_best], ]
        par.best  =  par.top[ 1, ]
        par.top   = par.top[2:n_best, ]
        
        rm( tracers )
    }
    
    return( list( input.parameters = input.parameters,  par.best = par.best, par.top = par.top, 
                  sim.top = sim.top, sim.best = sim.best, tracers_all = tracers_all, 
                  sim.tracers_all = sim.tracers_all, iKernelABC = iKernelABC ) )
}



# MSE ---------------------------------------------------------------------
# MSE is mean square error:
# for parameters (if true parameter is known):
MSE_parameters   <-   function( par.truth, par.top = NULL, par.best ){
    names( par.truth )  =  names( par.best )
    if ( nrow(par.best) > 1 ) stop( 'Please, use par.top for multi row data. par.best should be single row data.frame')
    # diff between par.top and par.truth
    if ( !is.null( par.top ) ){
        df  =  sapply( X = 1:ncol( par.truth ), FUN = function( x ) par.top[,x] - par.truth[,x] )
        mse_top  =  sapply( X = 1:nrow(df) , FUN = function( x ) sum( df[x,] ** 2 ) )
    } else { 
        mse_top = NULL 
    }
    
    # diff between par.best and par.truth
    df  =  sapply( X = 1:ncol( par.truth ), FUN = function( x ) ( par.best[,x] - par.truth[,x] ) )
    mse_best  =  sum( df ** 2 ) 
    
    if ( is.null( mse_top ) ) {
        return( list( mse_top = NULL, mse_best = mse_best ) ) 
    } else {
        return( list( mse_top = mean( mse_top ), mse_best = mse_best ) )
    }
}

# for output of simulation:
# MSE is mean square error:
MSE_sim   <-   function( stat.obs, stat.sim ){
    if ( nrow( stat.obs ) > 1 ) stop( 'Please, use stat.obs for each row data independently. 
                                            stat.obs should be single row data.frame')
    # diff between stat.obs and stat.sim
    df   =  sapply( X = 1:ncol( stat.obs ), FUN = function( x )  stat.sim[,x] - stat.obs[1,x]  )
    if ( !is.data.frame( df ) & !is.matrix( df ) ) df = as.data.frame( matrix( df, nrow = 1 ) )
    mse  =  sapply( X = 1:nrow(df) , FUN = function( x ) sum( df[x,] ** 2 ) )
    
    return( mse ) # mean( mse ) )
}

# PLOT --------------------------------------------------------------------

plot_sudoku_2D   <-  function( stat.sim, par.sim, stat.obs, iKernelABC, rslt, 
                              ind_X, ind_Y, names = c( 'X1', 'X2' ), draw_Y = FALSE, xlim, ylim,
                              show_tracer = TRUE, show_obs = TRUE, show_appropriate = TRUE, show_best = TRUE,
                              show_u_point = TRUE, show_legend = FALSE ){
    
    sbst_feature_Param  =  get_subset_of_feature_map( dtst  =  par.sim,
                                                      Matrix_Voronoi = iKernelABC$parameters_Matrix_Voronoi,
                                                      iFeature_point = iKernelABC$kernel_mean_embedding )
    if ( draw_Y ){
        sbst_feature_Y  =  get_subset_of_feature_map( dtst  =  stat.sim,
                                                  Matrix_Voronoi = iKernelABC$Matrix_Voronoi,
                                                  iFeature_point = iKernelABC$iFeature_point )
    }
    pr = par()
    l = par.sim[ , c( ind_X, ind_Y ) ]
    par( mgp = c(2.2, 0.5, 0), font.axis = 2, font.lab = 2 )
    plot( l[,1], l[,2], pch = 18, xlab = names[1], ylab = names[2], 
          axes = FALSE, cex.lab = 1.5, xlim = xlim, ylim = ylim )
          
          # cex.lab = 1.3 , cex.axis = 1.4 , mar = c(4,4,2,2),
          # tck = 0.04) ,
    # box()
    
    axis( 1, font = 2, tck = 0.03, cex.axis = 1.4 )
    axis( 2, font = 2, tck = 0.03, cex.axis = 1.4)    
    axis( 3, font = 2, tck = 0.03, cex.axis = 1.4, labels = FALSE )
    axis( 4, font = 2, tck = 0.03, cex.axis = 1.4, labels = FALSE ) 
    
    
    
    if ( draw_Y )    points( sbst_feature_Y[,ind_X], sbst_feature_Y[,ind_Y], col = 'red', cex = 1.6 )
    if ( show_u_point ) points( sbst_feature_Param[,ind_X], sbst_feature_Param[,ind_Y], col = 'blue', cex = 1.4, pch = 0 )
    if ( show_tracer ) points( rslt$tracer_bullets[,ind_X], rslt$tracer_bullets[,ind_Y], col = alpha('grey', 0.2 ), cex = 0.8, pch = 20 )
    if ( show_appropriate ) points( rslt$surroundings_best_points[,ind_X], rslt$surroundings_best_points[,ind_Y], 
            col = 'blue', cex = 0.5, pch = 20 )
    if ( show_best ) points( rslt$best_tracer_bullets[,ind_X], rslt$best_tracer_bullets[,ind_Y], col = 'red', cex = 0.8, pch = 5 )
    if (show_obs ) points( stat.obs[,ind_X], stat.obs[,ind_Y], col = 'red', cex = 2.4, pch = 4, lwd = 3 )
    
    if ( show_legend ) legend( x = 7.7, y = 10.3, bg = '#E0FFFF', box.col = 'grey',
                col = c('black', 'red', 'blue', 'grey', 'blue', 'red' ), 
                pch = c( 18,      4,       0 ,  1,    20,   5  ), bty = "o", 
                legend = c('Simulation','Truth observation', 
                           expression( paste('Voronoi sites ', upsilon[j]^'*' ) ),
                            'Tracers', 'iKernel > 0.5', 
                           'Top 20 points' )  
                )

}



plot_web_2D   <-  function( stat.sim, par.sim, stat.obs, iKernelABC, web, ind_X, ind_Y, 
                                names = c( 'P1', 'P2' ), draw_Y = FALSE, xlim, ylim,
                               show_tracer = TRUE, show_obs = TRUE, show_top = TRUE, show_best = TRUE,
                               show_u_point = TRUE, show_legend = FALSE ){
    
    sbst_feature_Param  =  get_subset_of_feature_map( dtst  =  par.sim,
                                                      Matrix_Voronoi = iKernelABC$parameters_Matrix_Voronoi,
                                                      iFeature_point = iKernelABC$kernel_mean_embedding )
    if ( draw_Y ){
        sbst_feature_Y  =  get_subset_of_feature_map( dtst  =  stat.sim,
                                                      Matrix_Voronoi = iKernelABC$Matrix_Voronoi,
                                                      iFeature_point = iKernelABC$iFeature_point )
    }
    pr = par()
    l = par.sim[ , c( ind_X, ind_Y ) ]
    par( mgp = c(2.2, 0.5, 0), font.axis = 2, font.lab = 2 )
    plot( l[,1], l[,2], pch = 18, xlab = names[1], ylab = names[2], 
          axes = FALSE, cex.lab = 1.5, xlim = xlim, ylim = ylim )
    
    # cex.lab = 1.3 , cex.axis = 1.4 , mar = c(4,4,2,2),
    # tck = 0.04) ,
    # box()
    
    axis( 1, font = 2, tck = 0.03, cex.axis = 1.4 )
    axis( 2, font = 2, tck = 0.03, cex.axis = 1.4)    
    axis( 3, font = 2, tck = 0.03, cex.axis = 1.4, labels = FALSE )
    axis( 4, font = 2, tck = 0.03, cex.axis = 1.4, labels = FALSE ) 
    
    
    
    if ( draw_Y )    points( sbst_feature_Y[,ind_X], sbst_feature_Y[,ind_Y], col = 'red', cex = 1.6 )
    if ( show_u_point ) points( sbst_feature_Param[,ind_X], sbst_feature_Param[,ind_Y], col = 'blue', cex = 1.4, pch = 0 )
    if ( show_tracer ) points( web$tracers_all[,ind_X], web$tracers_all[,ind_Y], col = alpha('grey', 0.2 ), cex = 0.8, pch = 20 )
    if ( show_top ) points( web$par.top[,ind_X], web$par.top[,ind_Y], 
                                    col = 'blue', cex = 0.5, pch = 20 )
    if ( show_best ) points( web$par.best[,ind_X], web$par.best[,ind_Y], col = 'red', cex = 0.8, pch = 5 )
    if (show_obs ) points( stat.obs[,ind_X], stat.obs[,ind_Y], col = 'red', cex = 2.4, pch = 4, lwd = 3 )
    
    if ( show_legend ) legend( x = 7.7, y = 10.3, bg = '#E0FFFF', box.col = 'grey',
                               col = c('black', 'red', 'blue', 'grey', 'blue', 'red' ), 
                               pch = c( 18,      4,       0 ,  1,    20,   5  ), bty = "o", 
                               legend = c('Simulation','Truth observation', 
                                          expression( paste('Voronoi sites ', upsilon[j]^'*' ) ),
                                          'Tracers', 'Top 20 points', 'The best point' )  
    )
    
}


# DRAFT_FUCTIONS -----------------------------------------------------------



### May be wrong:
Get_Mode <- function(v) {
    uniqv <- unique(v)
    return( uniqv[ which.max(tabulate(match(v, uniqv))) ] )
}
### May be wrong:
Mean_iKrnl  <- function( param, sm ){
    Mean_iKrnl = rep( 0, length( param ) )
    sum_sm  =  sum( sm ) 
    for (i in 1:length( param )) {
        Mean_iKrnl[i] = sum( sm * param[,i]) / sum_sm 
    }
    names(Mean_iKrnl) <- names( param)
    return( Mean_iKrnl )
}


