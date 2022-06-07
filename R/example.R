#' Example of simulation for lazy start
#'
#' @param verbose Logical type to show or do not show messages during execution
#' @param to_plot Logical type to plot or do not plot graphical results of a simulation
#' @param seed Numeric type to set seed for a simulation, if seed = NA (by default) then it will be skipped
#'
#' @return List of results of simulation with default values for all the parameters
#' @export
#'
#' @examples
#' NULL
#' # it takes a time for a simulation and then it will demonstrates results, \cr
#' # so, please, wait for a while
#' simulation_example( verbose = FALSE , to_plot = FALSE )
simulation_example  <-  function( verbose = TRUE , to_plot = TRUE, seed = NA ){

    if ( !is.na( seed ) ) set.seed( seed = seed )

    if ( verbose ) print('This code will be executed: ')
    if ( verbose ) print( simulation_example )

    # Attach packages from import list
    packages  =  list(  actuar = 'rztpois',
                        randomcoloR = 'randomColor',
                        methods = 'new',
                        stats = c('aggregate', 'rbinom', 'rexp', 'rnorm', 'runif' ),
                        stringr = c('str_length', 'str_split', 'str_sub', 'str_trim'),
                        utils = c('read.delim', 'read.table', 'write.table', 'globalVariables' ),
                        grDevices = c('dev.off', 'pdf', 'rgb'),
                        graphics = c('axis', 'legend', 'lines', 'par', 'plot', 'text', 'title' ) )

    for( pck in names( packages ) ){
        library( package = pck, character.only = TRUE, include.only = packages[[ pck ]])
    }


    return( NULL )
}
