
# Define global variables in tugHall.3:
# utils::globalVariables( c( 'CF', 'Compaction_factor' ) )



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
#' fl = system.file('extdata/Input', 'gene_map.txt',package = 'tugHall.3', mustWork = TRUE )
#' read_file(file_name = fl, stringsAsFactors = FALSE )
#' fl = system.file('extdata/Input', 'CF.txt',package = 'tugHall.3', mustWork = TRUE )
#' read_file(file_name = fl, stringsAsFactors = FALSE, header = FALSE )
read_file  <-  function( file_name = '', stringsAsFactors = FALSE, header = TRUE ){
    if ( file.size( file_name )  < 10 ) return( NULL )
    return( read.table( file = file_name, stringsAsFactors  =  stringsAsFactors ,
                        sep="\t", header = header ))
}


#' Function to copy the files by default from extdata folder in the library to Input folder in the working directory
#'
#' @param files Files to copy, vector of names of files by default: \cr
#' files = c( 'CCDS.current.txt', 'CF.txt', 'cloneinit.txt', 'gene_hallmarks.txt','gene_map.txt','parameters.txt' )
#' @param dir Folder to where files should be save, by default dir = 'Input'
#'
#' @return List of logic numbers for each copied file, TRUE - success, FALSE - not success
#' @export
#'
#' @examples
#' files = c('CF.txt', 'cloneinit.txt','gene_hallmarks.txt','gene_map.txt','parameters.txt' )
#' copy_files_to_Input( files, dir = 'Input' )
copy_files_to_Input  <-  function( files = c('CCDS.current.txt', 'CF.txt', 'cloneinit.txt',
                                             'gene_hallmarks.txt','gene_map.txt','parameters.txt' ) ,
                                   dir = 'Input' ){

    fls  =  lapply( X = files,
            FUN = function( x ) system.file('extdata/Input', x, package = 'tugHall.3', mustWork = TRUE ) )
    # fls  =  unlist( fls )

    if ( !file.exists( dir ) ) dir.create( dir )
    lapply( X = 1:length( fls ) , FUN = function( x ){
        file.copy( fls[[ x ]],  dir, overwrite = TRUE, recursive = TRUE, copy.mode = TRUE )
        } )

}



#' Function to copy the files of an example of simulation or from '/extdata/Output/' folder in the library to '/Output/' folder in the working directory
#'
#' @param files Files to copy, vector of names of files by default: \cr
#' files = c( 'cloneout.txt', 'CNA_mutations.txt', 'point_mutations.txt', 'gene_MAP.txt','geneout.txt','log.txt', 'order_genes_dysfunction.txt', 'VAF_data.txt', 'VAF.txt', 'weights.txt' )
#' @param dir Folder to where files should be save, by default dir = 'Output'
#'
#' @return List of logic numbers for each copied file, TRUE - success, FALSE - not success
#' @export
#'
#' @examples
#' files = c( 'cloneout.txt', 'CNA_mutations.txt', 'point_mutations.txt', 'gene_MAP.txt')
#' copy_files_to_Output( files )
#' files = c('geneout.txt','log.txt', 'VAF_data.txt', 'VAF.txt', 'weights.txt' )
#' copy_files_to_Output( files )
copy_files_to_Output  <-  function( files = c( 'cloneout.txt', 'CNA_mutations.txt', 'point_mutations.txt',
                                               'gene_MAP.txt','geneout.txt','log.txt', 'order_genes_dysfunction.txt',
                                               'VAF_data.txt', 'VAF.txt', 'weights.txt' ),
                                   dir = 'Output' ){

    fls  =  lapply( X = files,
                    FUN = function( x ) system.file('extdata/Output', x, package = 'tugHall.3', mustWork = TRUE ) )
    # fls  =  unlist( fls )

    if ( !file.exists( dir ) ) dir.create( dir )
    lapply( X = 1:length( fls ) , FUN = function( x ){
        file.copy( fls[[ x ]],  dir, overwrite = TRUE, recursive = TRUE, copy.mode = TRUE )
    } )

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



#' Function to change par for plots and after plotting it returns par values
#'
#' @param change_par_back logical. If TRUE it changes par options back after finishing a function call
#'
#' @return Change par() options and returns it (after finishing the function) to values before a function has started
#'
#' @export
#' @examples
#' define_par_for_plot( change_par_back = TRUE )
define_par_for_plot  <-  function( change_par_back ){

    if ( change_par_back ){
        op = par( no.readonly=TRUE )
        par(xpd=TRUE, cex.lab=2, lwd = 2, mar = c(5, 5, 5, 5),
            tcl = 0.5, cex.axis = 1.75,  mgp = c(3, 0.6, 0))
        on.exit( par( op ) )
    } else {
        par(xpd=TRUE, cex.lab=2, lwd = 2, mar = c(5, 5, 5, 5),
            tcl = 0.5, cex.axis = 1.75,  mgp = c(3, 0.6, 0))
    }
}

