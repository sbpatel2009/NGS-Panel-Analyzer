
#' Run the NGS Panel Shiny app
#'
#' This function runs the main shiny app. PRIOR to running this
#' function, the user will have needed to download the COSMIC tsv
#' file. That file must currently be called "cosmic.tsv" and live
#' in the current working directory.
#'
#' #@param cosmic_tsv The path to the manually-downloaded cosmic tsv
#' #    file.
#' 
#' @import shiny
#' @import rmarkdown
#' @import stringr
#' @import GenomicRanges
#' @import AnnotationHub
#' @import ggplot2
#'
#' 
#' @export
runNGSPanelApp = function() {
    if(!file.exists('cosmic.tsv')) {
        stop('COSMIC file must have been downloaded and present as "cosmic.tsv" in the R working directory')
    }
    shiny::runApp(system.file(package='NGSPanelAnalyzer','inst/app/app.R'))
}
        
