
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
#' that can be gotten by the function \code{get.MaxWiK()}
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
#' 
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
    
    tracer_bullets   =  unique.data.frame( tracer_bullets )
    
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
#' and this part is corresponding also to Voronoi sites/seeds
#'
#' @return The function \code{get_pairs_of_data_frame()} returns the list of the pairs of points
#' 
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
#' including given points as the first and the last rows
#' 
#' 
#'
#' @examples
#' NULL
generate_points_between_two_points  <-  function( pair, n = 10 ){
    ### Input
    ### pair is a data.frame of two points
    DF = pair
    DF[ 1:n,] = 0
    row.names( DF ) = 1:n
    for( i in 1:ncol( DF ) ){
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
#' so it is usually a subset of points corresponding to Voronoi sites/seeds
#' @param n_bullets Integer number of tracer points between each pair of points from DF
#'
#' @return The function \code{get_tracer_bullets()} returns data frame of generated tracer points
#' 
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


