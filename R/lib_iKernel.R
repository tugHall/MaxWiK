
# Initialization of data and norm definition ------------------------------

#' The norm function for vector
#'
#' @param x numeric vector
#'
#' @return The squared root of sum of squared elements of the vector x or Euclid length of the vector x
#'
#' @keywords internal
#'
#' @examples
#' NULL
norm_vec     <-  function(x) {
    sqrt(sum(x^2))
}

#' @describeIn norm_vec The squared norm or the sum of squared elements of the vector x
#' @return The squared Euclid norm or the sum of squared elements of the vector x
#' 
#' @keywords internal
#' 
#' @examples
#' NULL
norm_vec_sq  <-  function(x) {
    (sum(x^2))
}


# Isolation kernel --------------------------------------------------------
### ISOLATION KERNEL BASED ON VORONOI DIAGRAM 

### The common coefficients:
### t is a number of trees in iKernel or dimension of RKHS
### psi is a size of each Voronoi diagram
### nr is a size of dataset (number of rows in Dataset)

### The functions

#' The function to get feature representation in RKHS based on Voronoi diagram for WHOLE dataset
#'
#' @param psi Integer number related to the size of each Voronoi diagram
#' @param t Integer number of trees in Isolation Kernel or dimension of RKHS
#' @param data dataset of points, rows - points, columns - dimensions of a point
#' @param talkative logical. If TRUE then print messages, FALSE for the silent execution
#' @param Matrix_Voronoi Matrix of Voronoi diagrams, if it is NULL then the function will calculate Matrix_Voronoi
#'
#' @return Feature representation in RKHS based on Voronoi diagram for WHOLE dataset
#' 
#' @keywords internal
#'
#' @examples
#' NULL
get_voronoi_feature  <-  function( psi = 40, t = 350, data, talkative = FALSE, 
                                   Matrix_Voronoi = NULL ){
    ### Input parameters
    ### psi is number of points in a subset / dimension of RKHS 
    ### psi is also number of areas in any Voronoi diagram
    ### t is number of trees or number of Voronoi diagrams
    ### new is the logical parameter - Is Matrix_Voronoi new ?
    ### if new = FALSE then you should define Matrix_Voronoi
    
    new = is.null(Matrix_Voronoi)
    
    ### Check the data format:
    if (talkative ) cat( 'Check the data format')
    if ( talkative & check_numeric_format( data ) ) cat( 'OK' )
    if ( new & talkative ) cat( 'This is a new Voronoi diagram and Matrix') 
    if ( !new ) {
        if (talkative ) cat( 'This is a calculation based on a GIVEN Voronoi diagram and Matrix')
        psi = ncol( Matrix_Voronoi )
        t   = nrow( Matrix_Voronoi )
    }
    
    if (talkative ) cat( 'Finding matrix of distances between all points in data')
    dissim  <-  dist( data, method = 'euclidean', diag = TRUE, upper = TRUE)
    dissim  <-  as.matrix( dissim )
    nr  =  nrow(dissim)  ## number of columns and rows in matrix of distances
    if (talkative ) cat( 'Done' )
    
    if (talkative ) cat( 'Transform data to the Hilbert space related to Isolation Kernel' )
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
    if (talkative ) cat( 'Done' )
    
    return( list( M_Voronoi = Matrix_Voronoi, 
                  M_iKernel = Matrix_iKernel,
                  M_dissim  = dissim ) )
}

#' @describeIn get_voronoi_feature The function to get RKHS mapping based on Isolation Kernel for a new point
#'
#' @param d1 Data point - usually it is an observation data point
#' @param dissim Matrix of dissimilarity or distances between all points.
#' @param nr Integer number of rows in matrix of distances (dissim) and also the size of dataset
#'
#' @return RKHS mapping for a new point based on Isolation Kernel mapping
#' 
#' @keywords internal
#' 
#' @examples
#' NULL
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



#' The function to get subset with size psi for Voronoi diagram
#'
#' @param data_set Data.frame of Voronoi diagram
#' @param pnts Integer vector of indexes of columns of the data_set
#'
#' @return Subset of data_set with columns pnts
#' 
#' @keywords internal
#' 
#' @examples
#' NULL
#' 
GET_SUBSET  <-  function( data_set, pnts ){
    voronoi_subset  <-  data_set[  , pnts]
    return( voronoi_subset )
}

### The function to calculate a similarity between two points
### in the Matrix of real numbers.

#' Function returns the value of similarity or Isolation KERNEL for TWO points
#' 
#' @description \code{iKernel()} function returns value of similarity or Isolation KERNEL 
#' for TWO points that is number in the range \code{[0,1]}
#' 
#' @param Matrix_iKernel Matrix of indexes of Voronoi cells for each point and each tree based on Isolation Kernel calculation
#' @param pnt_1 The first point of dataset
#' @param pnt_2 The second point of dataset
#' @param t is a number of columns of Matrix_iKernel or dimension of Matrix_iKernel (corresponding to the number of trees t)
#' 
#' @return The function \code{iKernel()} returns a value of similarity or Isolation KERNEL for TWO points
#' 
#' @keywords internal
#'
#' @examples
#' NULL
iKernel  <-  function( Matrix_iKernel, pnt_1, pnt_2, t ){
    ### Input data:
    ### t is a number of columns of Matrix_iKernel
    ### pnt_1 and pnt_2 are IDs of points 1 and 2 in the Matrix_iKernel
    ### 
    
    # t  <-  ncol( Matrix_iKernel )
    smlrt  <- sum( Matrix_iKernel[ pnt_1, ]  == Matrix_iKernel[ pnt_2, ] ) / t
    return( smlrt )
}

#' @describeIn iKernel The function to get Isolation Kernel between a new point and dataset
#'
#' @description  \code{iKernel_point_dataset()} function returns vector of values of similarity based on Isolation Kernel between a new point and all the points of dataset
#' 
#' @param iFeature_point Feature mapping in RKHS for a new point, that can be gotten via \code{add_new_point_iKernel()} function
#' @param nr is number of rows in Matrix_iKernel or size of dataset
#' 
#' @return The function \code{iKernel_point_dataset()} returns a value of Isolation Kernel between a new point and dataset represented via Matrix_iKernel
#' 
#' @keywords internal
#' 
#' @examples
#' NULL
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


#' @describeIn iKernel The function to get weights from Feature mapping
#'
#' @description  \code{get_weights_iKernel()} function returns list of two objects: 
#' the first object is numeric vector of weights for RKHS space, and 
#' the second object is numeric vector of weights of similarity for iFeature_point 
#' corresponding observation point
#'
#' @param GI The inverse Gram matrix
#'
#' @return The function \code{get_weights_iKernel()} returns the 
#' list of weights for RKHS space and weights of similarity for iFeature_point
#' 
#' @keywords internal
#'
#' @examples
#' NULL 
get_weights_iKernel   <-   function( GI, Matrix_iKernel, t, nr, iFeature_point ){
    
    smlrt_p_data    =  iKernel_point_dataset( Matrix_iKernel, t, nr, iFeature_point )
    
    PHI_T_k_S       =  matrix( data = smlrt_p_data, ncol = 1 ) 
    
    WGHTS  =  GI  %*%  PHI_T_k_S
    
    return( list( weights_RKHS = WGHTS, weights_similarity = PHI_T_k_S ) )
}


#' The function to get feature representation in RKHS based on Voronoi diagram for PART of dataset
#'
#' @description \code{get_voronoi_feature_PART_dataset()} function returns 
#' the feature (mapping) representation in RKHS based on Voronoi diagram for NEW PART of dataset. 
#' The \code{Matrix_Voronoi} is based on the PREVIOUS dataset. 
#' The NEW PART of dataset will appear at the end of PREVIOUS dataset
#' 
#' @param data Data.frame of new points
#' @param talkative Logical parameter to print or do not print messages
#' @param start_row Row number from which a new data should be added
#' @param Matrix_Voronoi Matrix of Voronoi diagrams based on the PREVIOUS dataset
#' 
#' @return List of three matrices: Matrix_Voronoi, Matrix_iKernel and dissim
#' 
#' @keywords internal
#' 
#' @examples
#' NULL
get_voronoi_feature_PART_dataset  <-  function( data, talkative = FALSE, start_row,  Matrix_Voronoi ){
    ### Input parameters
    ### psi is a number of points in a subset / dimension of RKHS 
    ### psi is also number of areas in any Voronoi diagram
    ### t is number of trees or number of Voronoi diagrams
    ### start_row is a row number from which a new data have added
    ### Matrix_Voronoi is a matrix of Voronoi diagrams based on the PREVIOUS dataset
    
    ### Check the data format:
    if (talkative ) cat( 'Check the data format')
    if ( talkative & check_numeric_format( data ) ) cat( 'OK' )
    if (talkative ) cat( 'This is a calculation based on a GIVEN Voronoi diagram and Matrix')
    psi = ncol( Matrix_Voronoi )
    t   = nrow( Matrix_Voronoi )
    
    if (talkative ) cat( 'Finding matrix of distances between all points in data')
    dissim  <-  as.matrix( dist( data, method = 'euclidean', diag = TRUE, upper = TRUE) )
    nnw  =  nrow(dissim) - start_row  + 1  ### number of new points
    nr   =  nrow(dissim)                   ### number of rows 
    if (talkative ) cat( 'Done' )
    
    if (talkative ) cat( 'Transform a NEW data to the Hilbert space related to Isolation Kernel' )
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
    if (talkative ) cat( 'Done' )
    
    return( list( M_Voronoi = Matrix_Voronoi, 
                  M_iKernel = Matrix_iKernel,
                  M_dissim  = dissim ) )
}


#' The function to calculate Maxima weighted kernel mean mapping for Isolation Kernel in RKHS related to parameters space
#'
#' @param parameters_Matrix_iKernel Matrix of all the points represented in RKHS related to parameters space
#' @param Hilbert_weights Maximal weights in RKHS to get related part of kernel mean embedding from parameters_Matrix_iKernel
#'
#' @return Maxima weighted kernel mean mapping in the form of integer vector with length t (number of trees). 
#' Each element of the vector is index of Voronoi cell with maximal weight in the Voronoi diagram 
#' 
#' @keywords internal
#'
#' @examples
#' NULL
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


#' The function to get subset of points based on feature mapping
#'
#' @param dtst Dataset of all the original points
#' @param Matrix_Voronoi Matrix of Voronoi diagrams based on the Isolation Kernel algorithm
#' @param iFeature_point Feature mapping in RKHS for a point, 
#' that can be gotten via \code{add_new_point_iKernel()} function
#'
#' @return The subset of dtst that has points 
#' extracted with feature mapping of an observation point (iFeature_point)
#' 
#' @keywords internal
#'
#' @examples
#' NULL
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

#' @describeIn iKernel The function to calculate Gram matrix for Isolation Kernel method
#'
#' @description \code{GRAM_iKernel()} is the function to calculate Gram matrix for Isolation Kernel method based on Voronoi diagrams
#' 
#' @param check_pos_def Logical parameter to check the Gram matrix is positive definite or do not check
#'
#' @return  The function \code{GRAM_iKernel()} returns Gram matrix of Isolation Kernel
#' 
#' @keywords internal
#' 
#' @examples
#' NULL
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

#' The function to get inverse Gram matrix
#'
#' @description Function \code{get_inverse_GRAM()} allows to get inverse Gram matrix based on given
#' positive regularization constant lambda
#'  
#' @param G Gram matrix gotten via \code{GRAM_iKernel()} function
#' @param l Lambda parameter or positive regularization constant
#' @param check_pos_def Logical parameter to check the Gram matrix is positive definite or do not check
#'
#' @return Function \code{get_inverse_GRAM()} returns the inverse Gram matrix 
#' based on the given positive regularization constant lambda l
#' 
#' @keywords internal
#' 
#' @examples
#' NULL
get_inverse_GRAM  <-  function( G, l = 1E-6, check_pos_def = FALSE ){
    ### l is a regularization constant
    n = nrow( G )
    I = diag( n ) * n * l
    G = G - I
    d = det( G )
    # if ( d == 0 & !check_pos_def ) cat( 'The Gram matrix has determinant = 0 or close to zero' )    # stop( 'The Gram matrix has determinant = 0 \n' )
    
    GI  =  solve( G ) 
    if ( check_pos_def ) {
        if ( !check_positive_definite( GI, n = 50 ) ) stop( 'The Gram matrix is NOT positive definite !')
    }
    
    
    return( GI )
}


#' @describeIn get_inverse_GRAM The function to check the positive definite property of Gram matrix
#'
#' @description Function \code{check_positive_definite()} returns logical value about n trials on 
#' 'is Gram matrix positive definite or not?' Just incorrect trial returns FALSE
#' 
#' @param n Number of iterations to check the positive definite property
#'
#' @return Function \code{check_positive_definite()} returns logical value: \cr 
#' TRUE if Gram matrix is positive definite, and FALSE if it is not
#' 
#' @keywords internal
#' 
#' @examples
#' NULL
check_positive_definite  <-  function( G , n = 10 ){
    ### G - Gram matrix, 
    ### n - number of iterations to check the positive definite property
    nr = nrow( G )
    check_pf  =  TRUE
    
    for( i in 1:n ){
        rand_vect  =  as.matrix( runif( nr, -1, 0 )  )
        x = t( rand_vect)  %*%  G  %*%  rand_vect
        if ( as.numeric( x ) < 0 )  check_pf = FALSE
    }
    
    return( check_pf )
}




