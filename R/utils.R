
# Define global variables:
# utils::globalVariables( c( 'var1', 'Var2' ) )



#' Function to read file
#'
#' @param file_name Name of file to read
#' @param stringsAsFactors Parameter for read.table function, by default stringsAsFactors = FALSE
#' @param header Logical type to read or do not read head of a file
#'
#' @return data.frame of data from a file
#' @export
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


#' Function to make a large number of colors
#'
#' @param nm Number of colors
#'
#' @return Vector of colors with length more than nm
#'
#' @export
#' @examples
#' clrs = gen_colors( nm = 120 )
gen_colors  <-  function(nm = 12){
    # nm is a number of colors
    w <- (nm^(1/3)) %/% 1 +1

    st <-  w^3 %/% nm

    sq <- seq(0,1-1/w,1/w)

    cr <- 1:nm

    l <- 0
    R <- 1:(w^3)
    G <- R
    B <- R

    for (i in 1:w) {
        for (j in 1:w) {
            for (k in 1:w) {
                l <- l+1
                R[l] <- sq[i]
                G[l] <- sq[j]
                B[l] <- sq[k]
            }
        }
    }

    # seq(1,w^3,st) # is consequence of each color to make a high diversity of colors
    jColor <- data.frame( number = 1:length( seq( 1,w^3, st ) ),
                          color  = rgb( R[seq( 1, w^3, st ) ], G[seq( 1, w^3, st)],
                                        B[seq( 1, w^3, st ) ] ), stringsAsFactors = FALSE )

    return(jColor)

}



#' Check the installation of a package for some functions
#'
#' @param pkg Package name
#'
#' @return if the package is installed then it returns NULL else it returns error message
#'
#' @export
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
                        ggplot2       =  c( 'ggplot', 'geom_line', 'scale_color_manual', 'ggtitle', 'ylab', 'xlab' ),
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
