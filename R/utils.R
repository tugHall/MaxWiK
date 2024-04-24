

#' Function to read file
#'
#' @param file_name Name of file to read
#' @param stringsAsFactors Parameter for read.table function, by default stringsAsFactors = FALSE
#' @param header Logical type to read or do not read head of a file
#'
#' @return data.frame of data from a file
#' 
#'
#' @examples
#' # fl = system.file('extdata/Input', 'gene_map.txt',package = 'tugHall.3', mustWork = TRUE )
#' # read_file(file_name = fl, stringsAsFactors = FALSE )
#' # fl = system.file('extdata/Input', 'CF.txt',package = 'tugHall.3', mustWork = TRUE )
#' # read_file(file_name = fl, stringsAsFactors = FALSE, header = FALSE )
#' NULL
read_file  <-  function( file_name = '', stringsAsFactors = FALSE, header = TRUE ){
    if ( file.size( file_name )  < 10 ) return( NULL )
    return( read.table( file = file_name, stringsAsFactors  =  stringsAsFactors ,
                        sep="\t", header = header ))
}



#' Check the installation of a package for some functions
#'
#' @param pkg Package name
#'
#' @return if the package is installed then it returns NULL else it returns error message
#'
#' 
#' @examples
#' check_pkg( pkg = 'grDevices' )
check_pkg  <-  function( pkg ){
    msg  =  paste0( 'Package ', pkg, ' must be installed to use this function. \n ' )
    if ( !requireNamespace( pkg , quietly = TRUE ) )    stop( msg, call. = FALSE )
}


#' Check the installation of packages and attach them with corresponding functions
#'
#' @param pkgs List of package names with related function names, 
#' by default (or when pkgs = NULL) the list of packages are described in Namespace file of the package or
#' 'R/MaxWiK-package.R' file
#' 
#' @return if the packages are installed then it returns NULL else it returns error message
#'
#' @export
#' 
#' @examples
#' check_packages(  )
check_packages  <-  function( pkgs = NULL ){
    
    if ( is.null( pkgs ) ) {
        pkgs  =  list(  abc = 'abc',
                        bayestestR    =  'point_estimate',
                        ggplot2       =  c( 'ggplot', 'geom_line', 'scale_color_manual', 'ggtitle', 'ylab', 'xlab', 'aes' ),
                        graphics      =  c('axis', 'legend', 'lines', 'par', 'plot', 'text', 'title', 'points' ),
                        grDevices     =  c('dev.off', 'pdf', 'rgb'),
                        kernlab       =  c( 'kernelMatrix', 'laplacedot', 'rbfdot' ),
                        magrittr      =  '%>%',
                        methods       =  'new',
                        parallel      =  'mclapply',
                        plotly        =  c( 'plot_ly', 'add_trace', 'add_surface' ),
                        randomcoloR   =  'randomColor',
                        scales        =  'alpha',
                        stats         =  c('aggregate', 'rbinom', 'rexp', 'rnorm', 'runif', 'dist', 'complete.cases' ),
                        utils         =  c('read.delim', 'read.table', 'write.table', 'globalVariables' )
                        )
    }
    
    ### Attach the packages
    for( pck in names( pkgs ) ){
        check_pkg( pkg = pck )
        require( package = pck, character.only = TRUE, include.only = pkgs[[ pck ]])
    }
}



#' Function to check DATA.FRAME 
#'
#' Check that DATA.FRAME has numeric format for ALL the columns and it has NO 'NA' values
#' @param l DATA.FRAME that should have data of numeric type
#'
#' @return TRUE if data.frame has ONLY numeric data and FALSE vice verse  
#' 
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




#' Function to restrict data in the size to accelerate the calculations 
#' 
#' @description \code{restrict_data()} is based on rejection ABC method to restrict original dataset
#'
#' @param size Integer number of points to leave from original dataset
#' @param par.sim Data frame of parameters
#' @param stat.sim Data frame of outputs of simulations
#' @param stat.obs Data frame of observation point
#'
#' @return \code{restrict_data()} returns the list of: \cr
#' par.sim - restricted parameters which are close to observation point \cr
#' stat.sim - restricted stat.sim which are close to observation point
#' 
#' @export
#'
#' @examples
#' NULL
restrict_data  <-  function( par.sim, stat.sim, stat.obs, size = 300 ){
    l  =  nrow( par.sim )
    if ( l != nrow( stat.sim ) ) stop( "The parameters and statistics of simulations have different number of rows" )
    rej_abc  =  abc( target = stat.obs, param = par.sim, sumstat = stat.sim, 
                     tol = size / l, method = 'rejection' )
    # rej_abc$region
    new_par  =   par.sim[ rej_abc$region, ]
    new_sim  =  stat.sim[ rej_abc$region, ]
    return( list( par.sim   =  new_par, 
                  stat.sim  =  new_sim ) )
}



#' Function to copy the templates from extdata folder in the library to /Templates/ folder in the working directory
#'
#' @param dir Folder to where files should be save, by default dir = './'
#'
#' @return List of logic numbers for each copied file, TRUE - success, FALSE - not success
#' @export
#'
#' @examples
#' copy_templates( dir = 'Input' )
copy_templates  <-  function( dir = './' ){
    
    dir_pck =  system.file('extdata', 'templates', package = 'MaxWiK', mustWork = TRUE )
    files  =  list.files( dir_pck )
    fls  =  lapply( X = files,
                    FUN = function( x ) system.file('extdata/templates', x, package = 'MaxWiK', mustWork = TRUE ) )
    
    if ( !file.exists( dir ) ) dir.create( dir )
    lapply( X = 1:length( fls ) , FUN = function( x ){
        file.copy( fls[[ x ]],  dir, overwrite = TRUE, recursive = TRUE, copy.mode = TRUE )
    } )
    
}
