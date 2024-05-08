# Kernel ABC --------------------------------------------------------------

#' Function to get Approximate Bayesian Computation based on Maxima Weighted Isolation Kernel mapping
#'
#' @description The function \code{get.MaxWiK()} is used to get Approximate Bayesian Computation 
#' based on Maxima Weighted Isolation Kernel mapping.
#' On given data frame of parameters, statistics of the simulations and an observation, 
#' using the internal parameters psi and t, 
#' the function \code{get.MaxWiK()} returns the estimation of a parameter corresponding to
#' Maxima weighted Isolation Kernel ABC method. 
#'
#' @param psi Integer number. Size of each Voronoi diagram or number of areas/points in the Voronoi diagrams
#' @param t Integer number of trees in the Isolation Forest
#' @param param or \code{par.sim} - data frame of parameters of the model
#' @param stat.sim Summary statistics of the simulations (model output)
#' @param stat.obs Summary statistics of the observation point
#' @param talkative Logical parameter to print or do not print messages
#' @param check_pos_def Logical parameter to check the Gram matrix is positive definite or do not check
#' @param Matrix_Voronoi is a predefined matrix of information about Voronoi trees 
#' (rows - trees, columns - Voronoi points/areas IDs). By default it is NULL and is generated randomly.
#' 
#' @return The function \code{get.MaxWiK()} returns the list of :
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
get.MaxWiK  <-  function( psi = 40, t = 350, param, 
                          stat.sim, stat.obs, talkative = FALSE, 
                          check_pos_def = TRUE, Matrix_Voronoi = NULL ){
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
                                   talkative = talkative, Matrix_Voronoi = Matrix_Voronoi )
    
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
                                         Matrix_Voronoi = Matrix_Voronoi )
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

