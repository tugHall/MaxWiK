

#' Function to read file
#'
#' @param file_name Name of file to read
#' @param stringsAsFactors Parameter for read.table function, by default stringsAsFactors = FALSE
#' @param header Logical type to read or do not read head of a file
#'
#' @returns data.frame of data from a file
#' 
#' @export
#'
#' @examples
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
#' @returns if the package is installed then it returns NULL else it returns error message
#' 
#' @keywords internal
#' 
#' @examples
#' NULL
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
#' @returns if the packages are installed then it returns NULL else it returns error message
#'
#' @keywords internal
#' 
#' @examples
#' NULL
check_packages  <-  function( pkgs = NULL ){
    
    if ( is.null( pkgs ) ) {
        pkgs  =  list(  abc = 'abc',
                        ggplot2       =  c( 'aes', 'aes_string', 'geom_density', 'geom_line', 'geom_vline', 'ggplot', 'ggtitle', 'scale_color_manual', 'ylab', 'xlab' ),
                        methods       =  'new',
                        parallel      =  'mclapply',
                        scales        =  'alpha',
                        stats         =  c( 'complete.cases', 'dist', 'runif' ),
                        utils         =  c( 'read.table', 'setTxtProgressBar', 'txtProgressBar', 'write.table' )
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
#' @returns TRUE if data.frame has ONLY numeric data and FALSE vice verse  
#' 
#' @keywords internal
#' 
#' @examples
#' NULL
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
#' @returns \code{restrict_data()} returns the list of: \cr
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
#' @returns List of logic numbers for each copied file, TRUE - success, FALSE - not success
#' @export
#'
#' @examples
#' MaxWiK_templates( dir = './' )
MaxWiK_templates  <-  function( dir = './' ){
    
    dir_pck =  system.file('extdata', 'templates', package = 'MaxWiK', mustWork = TRUE )
    files  =  list.files( dir_pck )
    fls  =  lapply( X = files,
                    FUN = function( x ) system.file('extdata/templates', x, package = 'MaxWiK', mustWork = TRUE ) )
    
    if ( !file.exists( dir ) ) dir.create( dir )
    lapply( X = 1:length( fls ) , FUN = function( x ){
        file.copy( fls[[ x ]],  dir, overwrite = TRUE, recursive = TRUE, copy.mode = TRUE )
    } )
    
}



#' The function to get the mean square error values for statistics of simulations
#'
#' @description The function \code{MSE_sim()} allows to get 
#' the mean square error values for statistics of simulations
#'
#' @param stat.sim Summary statistics of the simulations (model output)
#' @param stat.obs Summary statistics of the observation point
#'
#' @returns The function \code{MSE_sim()} returns numeric vector of
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
#' @param par.top Parameters from the top of similarities of \code{get.MaxWiK()} algorithm
#' @param par.best The best parameter from \code{get.MaxWiK()} algorithm
#'
#' @returns The function \code{MSE_parameters()} returns list of two numbers: \cr
#' - mean of MSE values for all the points from par.top; \cr
#' - MSE value for the point of par.best 
#' 
#' @keywords internal
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



#' Function to read hyperparameters and their values from the file
#'
#' @param input File name to input 
#'
#' @returns Parameters and their values
#' 
#' @export
#'
#' @examples
#' NULL
read_hyperparameters  <- function( input ){
    
    check_packages()
    input.table <- read.table( input, header=F, sep="\t" )
    
    out = list()
    suppressWarnings( expr = {
        for( ii in 1:nrow( input.table ) ) {
            V1 <- input.table[ ii, 1 ]
            V2 <- input.table[ ii, 2 ]
            
            # add quotes to string 
            if ( is.na( as.numeric( V2 ) ) & 
                 is.na( as.logical( V2 ) ) &
                 V2 != 'list()'
            ) V2 <- paste("\'", V2, "\'", sep='')
            text <- paste( 'out$', V1, ' <- ', V2, sep='' )
            eval( parse( text = text ) )
    }
    }, classes = "warning" )
    
    # out  =  sapply(ls(),function(x)get(x),simplify=F,USE.NAMES=T)
    
    return( out )
}



#' Function to restrict values of the data according with the range for each dimension
#'
#' @param diapason Vector of min and max values or data frame with two rows (min and max) for each dimension of input data
#' @param input.data Data frame of input where values will be corrected
#'
#' @returns The same data frame with corrected values according to the diapason
#' 
#' @export
#'
#' @examples
#' NULL
apply_range  <- function( diapason, input.data ){
    
    if ( is.vector( diapason) ) if ( length( diapason ) != 2 ) 
        stop( 'diapason should be vector with length = 2' )
    if ( is.data.frame( diapason ) ) if ( nrow( diapason) !=2 ) 
        stop( 'diapason should be data frame with number of rows = 2 for min and max values for each dimension' )
    if ( is.null( input.data ) ) return( NULL)
    
    out.data  =  input.data
    
    if ( !is.data.frame( input.data ) ){
        if ( !is.vector( input.data ) ) stop( 'input.data should be data frame or vector.')
        if ( !is.vector( diapason ) )   stop( 'If input data is vector, diapason also shopuld be vector.')
        w = which( input.data > diapason[ 2 ] )
        if ( length( w ) > 0 ) out.data[ w ] = diapason[ 2 ] 
        w = which( input.data < diapason[ 1 ] )
        if ( length( w ) > 0 ) out.data[ w ] = diapason[ 1 ] 
        return( out.data )
    }
    
    f_minmax = function( icol, input.data, diapason ) {
                out.data = input.data
                if ( is.vector( diapason ) ) {
                    min_max = diapason
                } else {
                    min_max = diapason[ , icol ]
                }
                # Max check
                w = which( input.data[ , icol ] > min_max[ 2 ] )
                if ( length( w ) > 0 ) out.data[ w, icol ]  =  min_max[ 2 ]
                # Min check
                w = which( input.data[ , icol ] < min_max[ 1 ] )
                if ( length( w ) > 0 ) out.data[ w, icol ]  =  min_max[ 1 ]
                return( out.data[ , icol] )
    }
    
    for( icol in 1:ncol( input.data ) ){
        out.data[ , icol ]  =  f_minmax( icol = icol, input.data = input.data, diapason = diapason )
    }
    return( out.data )
}


