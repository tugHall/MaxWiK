
# Initialization of data and norm definition ------------------------------

#' The norm function for vector
#'
#' @param x numeric vector
#'
#' @return The squared root of sum of squared elements of the vector x or Euclid length of the vector x
#' @export
#'
#' @examples
#' norm_vec(c(3,4)) # that returns 5
#' norm_vec(c(12,5)) # that returns 13
#' 
norm_vec     <-  function(x) {
    sqrt(sum(x^2))
}

#' @describeIn norm_vec The squared norm or the sum of squared elements of the vector x
#' @return The squared Euclid norm or the sum of squared elements of the vector x
#' 
#' @export
#' @examples
#' norm_vec_sq(c(3,4)) # that returns 25
#' norm_vec_sq(c(12,5)) # that returns 169
norm_vec_sq  <-  function(x) {
    (sum(x^2))
}


#' Function to check DATA.FRAME 
#'
#' Check that DATA.FRAME has numeric format for ALL the columns and it has NO 'NA' values
#' @param l DATA.FRAME that should have data of numeric type
#'
#' @return TRUE if data.frame has ONLY numeric data and FALSE vice verse  
#' @export
#'
#' @examples
#' \dontrun{
#' check_numeric_format( data.frame( A= c(9,0), B = c(4,6)) )  # TRUE
#' check_numeric_format( data.frame( A= c(9,0), B = c(4,NA)) )  # Error due to NA value
#' check_numeric_format( data.frame( A= c(9,'0'), B = c(4,6)) ) # Error due to character in the data
#' }
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

#' Kernel function for two vectors x1 and x2
#'
#' @param x1 Numeric vector x1
#' @param x2 Numeric vector x2
#' @param sgm Sigma parameter in Gaussian kernel
#'
#' @return Value of Gaussian kernel for two vectors and given sigma
#' @export
#'
#' @examples
#' Gaussian_kernel(x1 = c(3,4), x2 = c(3.4,5), sgm = 1 )
#' Gaussian_kernel(x1 = c(3,4), x2 = c(3.4,5), sgm = 5 )
#' Gaussian_kernel(x1 = c(3,4), x2 = c(4,5), sgm = 1 )
Gaussian_kernel  <-  function( x1, x2 , sgm = 1){
    exp( -1 * norm_vec_sq(x1 - x2) / (2 * sgm**2 ) )
}

#' @describeIn Gaussian_kernel Value of Laplacian kernel for two vectors and given width
#' @param b The width of Laplacian kernel
#'
#' @return Value of Laplacian kernel for two vectors and given width
#' @export
#'
#' @examples
#' Laplacian_kernel(x1 = c(3,4), x2 = c(4,5), b = 1) 
#' Laplacian_kernel(x1 = c(3,4), x2 = c(4,5), b = 5) 
#' Laplacian_kernel(x1 = c(3,4), x2 = c(3,5), b = 1) 
Laplacian_kernel  <-  function( x1, x2 , b = 1 ){
    exp( -1 * norm_vec(x1 - x2)    / (2 * b ) )
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
#' @param new logical. Is Matrix_Voronoi new ? If TRUE then function will calculate Matrix_Voronoi, if FALSE function will use input Matrix_Voronoi.
#' @param Matrix_Voronoi Matrix of Voronoi diagrams that is used only if new = FALSE
#'
#' @return Feature representation in RKHS based on Voronoi diagram for WHOLE dataset
#' @export
#'
#' @examples
#' NULL
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

#' @describeIn get_voronoi_feature The function to get RKHS mapping based on Isolation Kernel for a new point
#'
#' @param d1 Data point - usually it is an observation data point
#' @param dissim Matrix of dissimilarity or distances between all points.
#' @param nr Integer number of rows in matrix of distances (dissim) and also the size of dataset
#'
#' @return RKHS mapping for a new point based on Isolation Kernel mapping
#' @export
#'
#' @examples
#' 
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
#' @export
#'
#' @examples
#' 
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
#' @export
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
#' @export
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
#' @export 
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
#' @export
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


#' The function to calculate Maxima weighted kernel mean mapping for Isolation Kernel in RKHS related to parameters space
#'
#' @param parameters_Matrix_iKernel Matrix of all the points represented in RKHS related to parameters space
#' @param Hilbert_weights Maximal weights in RKHS to get related part of kernel mean embedding from parameters_Matrix_iKernel
#'
#' @return Maxima weighted kernel mean mapping in the form of integer vector with length t (number of trees). 
#' Each element of the vector is index of Voronoi cell with maximal weight in the Voronoi diagram 
#' @export
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
#' @export
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
#' @export
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
#' @export
#'
#' @examples
#' NULL
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


#' @describeIn get_inverse_GRAM The function to check the positive definite property of Gram matrix
#'
#' @description Function \code{check_positive_definite()} returns logical value about n trials on 
#' 'is Gram matrix positive definite or not?' Just incorrect trial returns FALSE
#' 
#' @param n Number of iterations to check the positive definite property
#'
#' @return Function \code{check_positive_definite()} returns logical value: \cr 
#' TRUE if Gram matrix is positive definite, and FALSE if it is not
#' @export
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
        # print( paste( as.numeric( x ) , check_pf, sep = ' - ' ) )
    }
    
    return( check_pf )
}



# Kernel ABC --------------------------------------------------------------

#' Function to get Approximate Bayesian Computation based on Maxima Weighted Isolation Kernel mapping
#'
#' @description The function \code{iKernelABC()} is used to get Approximate Bayesian Computation 
#' based on Maxima Weighted Isolation Kernel mapping.
#' On given data frame of parameters, statistics of the simulations and an observation, 
#' using the internal parameters psi and t, 
#' the function \code{iKernelABC()} returns the estimation of a parameter corresponding to
#' Maxima weighted Isolation Kernel ABC method. 
#'
#' @param psi Integer number. Size of each Voronoi diagram or number of areas/points in the Voronoi diagrams
#' @param t Integer number of trees in the Isolation Forest
#' @param param Data frame of parameters of the model
#' @param stat.sim Summary statistics of the simulations (model output)
#' @param stat.obs Summary statistics of the observation point
#' @param talkative Logical parameter to print or do not print messages
#' @param check_pos_def Logical parameter to check the Gram matrix is positive definite or do not check
#'
#' @return The function \code{iKernelABC()} returns the list of :
#' - kernel_mean_embedding is a maxima weighted kernel mean embedding (mapping) related to the observation point;
#' - parameters_Matrix_Voronoi is a matrix of information about Voronoi trees (rows - trees, columns - Voronoi points/areas IDs) for parameters data set;
#' - parameters_Matrix_iKernel is a matrix of of all points of PARAMETERS in a Hilbert space (rows - points, columns - isolation trees);
#' - Hilbert_weights is a weights in Hilbert space to get maxima weighted kernel mean embedding for parameters_Matrix_iKernel;
#' - Matrix_iKernel is a matrix of all points of simulations in a Hilbert space (rows - points, columns - isolation trees);
#' - iFeature_point is a feature embedding mapping for the OBSERVATION point;
#' - similarity is a vector of similarities between the simulation points and observation point;
#' - Matrix_Voronoi is a matrix of information about Voronoi trees (rows - trees, columns - Voronoi points/areas IDs);
#' - t is a number of trees in the Isolation Forest; 
#' - psi is a number of areas/points in the Voronoi diagrams
#' 
#' @export
#'
#' @examples
#' NULL
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

#' The function to get the best tracer bullets related to kernel mean embedding
#' 
#' @description The function \code{sudoku()} allows to get the best tracer bullets related to kernel mean embedding.
#' The calculation performs ONLY for parameters dataset DT = par.sim. 
#' This function performs a heuristic algorithm to seek a space/area related to 
#' the feature mapping in Hilbert space for the dataset of the parameters. \cr
#' The main idea of the algorithm is just: \cr
#' 1) Generate points between the centers of Voronoi diagrams related to 
#' the Maxima weighted feature mapping based on Isolation Kernel \cr
#' 2) Following strategy to puzzle out of SUDOKU: delete all points that do not match feature mapping \cr
#' 3) Output: The remaining points should be corresponding to the feature mapping.  
#' 
#' @param DT Whole dataset of parameters
#' @param iKernelABC Result of calculations based on Isolation Kernel ABC 
#' that can be gotten by the function \code{iKernelABC()}
#' @param n_bullets Integer number of tracer bullets / additional points between the TWO most distant points
#' @param n_best Integer number of the best tracer bullets / points 
#' to consider them at the next algorithmic step
#' @param halfwidth Criterion to choose the best tracer points like: \cr
#' \code{if similarity_of_point >= halfwidth} then it is the point to be included to the pool of the best points 
#' 
#' @return The function \code{sudoku()} returns the list of next objects:
#' - tracer_bullets that is all the points generated during the run of the algorithm, 
#' - criterion that is a value of the similarity that is used to choose the best tracer points,  
#' - best_tracer_bullets that is the best tracer points that have similarity more or equal than \strong{criterion} value,
#' - surroundings_best_points that is the best tracer points that have similarity more or equal than \strong{halfwidth} value, 
#' - feature_tracers that is results of the function \code{get_voronoi_feature_PART_dataset()} applied to the new tracer points, 
#' - similarity_to_mean that is numeric vector of similarities of all the tracers points.
#' 
#' 
#' @export
#'
#' @examples
#' NULL
#' 
sudoku  <-  function( DT , iKernelABC, n_bullets = 20, n_best = 10, halfwidth = 0.5 ){
    ### Input:
    ### DT is a data set of parameters
    ### iKernelABC is a result of calculations based on Isolation Kernel ABC
    ### n_bullets is a number of tracer bullets / additional points between the TWO farthest distance points
    ### n_best is a number of the best tracer bullets / points
    
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
    
    ### Lets take the best points:
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

#' @describeIn sudoku The function to get pairs from Data Frame
#'
#' @description The function \code{get_pairs_of_data_frame()} is used to get pairs of points 
#' from the Data Frame that is the most distant each other.
#' In other words, the algorithm seeks the most distant coupled point to each point from the data frame
#'
#' @param DF Data Frame that is usually part of whole data frame of parameters, 
#' and this part is corresponding also to Voronoi cites/seeds
#'
#' @return The function \code{get_pairs_of_data_frame()} returns the list of the pairs of points
#' @export
#'
#' @examples
#' NULL 
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


#' @describeIn sudoku The function to generate points between the pair of given points
#'
#' @description The function \code{generate_points_between_two_points()} is used to generate 
#' points between two given points
#'
#' @param pair Data frame of two points
#' @param n Integer number of points that should be located between two input points
#'
#' @return The function \code{generate_points_between_two_points()} 
#' returns data frame of generated points between two given points, 
#' including them as the first and last rows
#' @export
#'
#' @examples
#' NULL
generate_points_between_two_points  <-  function( pair, n = 10 ){
    ### Input
    ### pair is a data.frame of two points
    DF = pair
    DF[ 1:n,] = 0
    row.names( DF ) = 1:n
    for( i in 1:length( DF ) ){
        dlt = ( pair[ 2, i ] - pair[ 1, i ] ) / ( n - 1 )  
        DF[ , i ]  =  pair[ 1, i ] + dlt * ( 0 : ( n - 1 ) )
    }
    return( DF )
}


#' @describeIn sudoku The function to get 'tracer bullets' or tracer points 
#'
#' @description The function \code{get_tracer_bullets()} is used to to get 'tracer bullets' or tracer points 
#' generated between all the pairs of the most distant points
#'
#' @param DF Data frame of oints that is used for generation of tracer points, 
#' so it is usually a subset of points corresponding to Voronoi cites/seeds
#' @param n_bullets Integer number of tracer points between each pair of points from DF
#'
#' @return The function \code{get_tracer_bullets()} returns data frame of generated tracer points
#' @export
#'
#' @examples
#' NULL 
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









# SPIDERWEB algorithm -----------------------------------------------------

#' @describeIn iKernelABC The function to get the best value of parameter corresponding to 
#' Maxima Weighted Isolation Kernel mapping which is related to an observation point
#' 
#' @description The function \code{spiderweb()} itteratively generates tracer points gotten 
#' from \code{sudoku()} algorithm, based on the simple procedure: \cr
#' - making a reflection of the top points from the best point, \cr 
#' - and then generating the point tracers between them, \cr
#' - finally, the algorithm chooses again the top points and the best point (\code{sudoku()} function is used),
#' - repeat all the steps until condition to be \code{TRUE}: \cr
#' \code{abs( max( sim_tracers ) - sim_previous ) < epsilon }
#'
#' @param n_bullets Integer number of tracer bullets / additional points between the TWO most distant points
#' @param n_best Integer number of the best tracer bullets / points 
#' to consider them at the next algorithmic step
#' @param halfwidth Criterion to choose the best tracer points like: \cr
#' \code{if similarity_of_point >= halfwidth} then it is the point to be included to the poool of the best points
#' @param epsilon Criterion to stop algorithm \code{spiderweb()} that isused to check: \cr
#' \code{if ( abs( max( sim_tracers ) - sim_previous ) < epsilon ) break}
#'
#' @return The function \code{spiderweb()} returns the list of the next objects:
#' - input.parameters the list of all the input parameters for Isolation Kernel ABC method;
#' - par.best that is data frame of one point that is the best from all the generated tracer points;
#' - par.top that is data frame of n_best points that are the top from all the generated tracer points; 
#' - sim.top that is numeric vecor of similarities of the top points;
#' - sim.best that is numeric value of the similarity of the best tracer point;
#' tracers_all that is data frame of all the generated tracer points; 
#' - sim.tracers_all that is numeric vector of similarities of all the generated tacer points;
#' - iKernelABC that is result of the function \code{iKernelABC()} given on \code{input parameters}.
#' 
#' @export
#' 
#' @examples
#' NULL
spiderweb_old  <-  function( psi = 4, t = 35, param = param, 
                         stat.sim = stat.sim, stat.obs = stat.obs, 
                         talkative = FALSE, check_pos_def = FALSE ,
                         n_bullets = 5, n_best = 10, halfwidth = 0.5, 
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
    
    par.top   =  tracers[ order( rslt$similarity_to_mean , decreasing = TRUE)[1:n_best], ]
    par.best  =  par.top[ 1, ]
    par.top   =  par.top[2:n_best, ]
    rm( tracers )
    sim_previous  =  0  # max( rslt$similarity_to_mean )
    
    sim.top   =  NULL
    sim.best  =  -1
    while( TRUE ){
        ### Reflect par.top through par.best 
        par.reflect  =  par.top
        for( i in 1:nrow( par.top ) )  par.reflect[ i, ]  =  2 * par.best - par.top[ i ,  ] 
        
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

#' @describeIn iKernelABC The function to get the best value of parameter corresponding to 
#' Maxima Weighted Isolation Kernel mapping which is related to an observation point
#' 
#' @description The function \code{spiderweb()} itteratively generates tracer points gotten 
#' from \code{sudoku()} algorithm, based on the simple procedure: \cr
#' - making a reflection of the top points from the best point, \cr 
#' - and then generating the point tracers between them, \cr
#' - finally, the algorithm chooses again the top points and the best point (\code{sudoku()} function is used),
#' - repeat all the steps until condition to be \code{TRUE}: \cr
#' \code{abs( max( sim_tracers ) - sim_previous ) < epsilon }
#'
#' @param n_bullets Integer number of tracer bullets / additional points between the TWO most distant points
#' @param n_best Integer number of the best tracer bullets / points 
#' to consider them at the next algorithmic step
#' @param halfwidth Criterion to choose the best tracer points like: \cr
#' \code{if similarity_of_point >= halfwidth} then it is the point to be included to the poool of the best points
#' @param epsilon Criterion to stop algorithm \code{spiderweb()} that isused to check: \cr
#' \code{if ( abs( max( sim_tracers ) - sim_previous ) < epsilon ) break}
#'
#' @return The function \code{spiderweb()} returns the list of the next objects:
#' - input.parameters the list of all the input parameters for Isolation Kernel ABC method;
#' - par.best that is data frame of one point that is the best from all the generated tracer points;
#' - par.top that is data frame of n_best points that are the top from all the generated tracer points; 
#' - sim.top that is numeric vecor of similarities of the top points;
#' - sim.best that is numeric value of the similarity of the best tracer point;
#' tracers_all that is data frame of all the generated tracer points; 
#' - sim.tracers_all that is numeric vector of similarities of all the generated tacer points;
#' - iKernelABC that is result of the function \code{iKernelABC()} given on \code{input parameters}.
#' 
#' @export
#' 
#' @examples
#' NULL
spiderweb  <-  function( psi = 4, t = 35, param = param, 
                             stat.sim = stat.sim, stat.obs = stat.obs, 
                             talkative = FALSE, check_pos_def = FALSE ,
                             n_bullets = 5, n_best = 10, halfwidth = 0.5, 
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
    
    par.top   =  tracers[ order( rslt$similarity_to_mean , decreasing = TRUE)[1:n_best], ]
    par.best  =  par.top[ 1, ]
    par.top   =  par.top[2:n_best, ]
    rm( tracers )
    sim_previous  =  0  # max( rslt$similarity_to_mean )
    
    sim.top   =  NULL
    sim.best  =  -1
    while( TRUE ){
        ### Reflect par.top through par.best 
        par.reflect  =  par.top
        for( i in 1:nrow( par.top ) )  par.reflect[ i, ]  =  2 * par.best - par.top[ i ,  ] 
        
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


#' The function to get the mean square error values for statistics of simulations
#'
#' @description The function \code{MSE_sim()} allows to get 
#' the mean square error values for statistics of simulations
#'
#' @param stat.sim Summary statistics of the simulations (model output)
#' @param stat.obs Summary statistics of the observation point
#'
#' @return The function \code{MSE_sim()} returns numeric vector of
#' the mean square error values for statistics of simulations
#' 
#' @export
#'
#' @examples
#' NULL
MSE_sim   <-   function( stat.obs, stat.sim ){
    if ( nrow( stat.obs ) > 1 ) stop( 'Please, use stat.obs for each row data independently. 
                                            stat.obs should be single row data.frame')
    # diff between stat.obs and stat.sim
    df   =  sapply( X = 1:ncol( stat.obs ), FUN = function( x )  stat.sim[,x] - stat.obs[1,x]  )
    if ( !is.data.frame( df ) & !is.matrix( df ) ) df = as.data.frame( matrix( df, nrow = 1 ) )
    mse  =  sapply( X = 1:nrow(df) , FUN = function( x ) sum( df[x,] ** 2 ) )
    
    return( mse ) # mean( mse ) )
}


#' @describeIn MSE_sim The function calculates mean square error (MSE) value 
#' for parameters as differences between them and already the known truth parameter 
#'
#' @description The function \code{MSE_parameters()} allows to get MSE for parameters if the truth parameter is known
#'
#' @param par.truth The truth parameter
#' @param par.top Parameters from the top of similarities of \code{iKernelABC()} algorithm
#' @param par.best The best parameter from \code{iKernelABC()} algorithm
#'
#' @return The function \code{MSE_parameters()} returns list of two numbers: \cr
#' - mean of MSE values for all the points from par.top; \cr
#' - MSE value for the point of par.best 
#' 
#' @export
#' 
#' @examples
#' NULL 
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



# DRAFT_FUCTIONS -----------------------------------------------------------



### To check correctness:

#' The function to get mode of the input vector
#' 
#' @description The function \code{Get_Mode} returns the mode of vector. \cr 
#' The mode is the value that appears most often in a dataset. 
#' 
#' @param v Input vector of the integer, numeric numbers or characters 
#'
#' @return The function \code{Get_Mode} returns the mode of vector
#' 
#' @export
#'
#' @examples
#' NULL 
Get_Mode <- function(v) {
    uniqv <- unique(v)
    return( uniqv[ which.max(tabulate(match(v, uniqv))) ] )
}



#' The function returns the weighted mean of the parameter based on Isolation Kernel ABC method
#'
#' @description The function \code{Mean_iKernel_parameters()} returns the weighted mean of the parameter 
#' that was calculated with \code{KernelABC()} function that represents Isolation Kernel ABC method 
#' 
#' @param param Data frame of parameters
#' @param sm Numeric vector of weights gotten from \code{iKernelABC()} function
#'
#' @return The function \code{Mean_iKernel_parameters()} returns the weighted mean of the parameter 
#' @export
#'
#' @examples 
#' NULL 
Mean_iKernel_parameters  <- function( param, sm ){
    Mean_iKrnl = rep( 0, length( param ) )
    sum_sm  =  sum( sm ) 
    for (i in 1:length( param )) {
        Mean_iKrnl[i] = sum( sm * param[,i]) / sum_sm 
    }
    names(Mean_iKrnl) <- names( param )
    return( Mean_iKrnl )
}


