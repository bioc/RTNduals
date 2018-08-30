##------------------------------------------------------------------------------
setClassUnion("TNInull", members = c("TNI", "NULL"))

#' MBR objects
#' 
#' @description 
#' MBR: an S4 class for co-regulation analysis and inference of 'dual 
#' regulons'. 
#'
#' @details 
#' The MBR class is a container for results from the MBR methods. 
#' The class slots are used to store information of different transcriptional networks, 
#' regulator annotation, infered 'dual regulons' and parameters used in the analysis. 
#' All the information is stored in nine slots.
#' 
#' @slot TNI a 'TNI' object created by the RTN package.
#' @slot regulatoryElements regulatory elements listed in the TNI.
#' @slot dualRegulons all possible 'duals regulons' computed by 
#' \code{\link{mbrAssociation}}
#' @slot results a list, results from the MBR methods.
#' @slot para a list, parameters used in the MBR methods.
#' @slot summary a list, summary for 'para' and 'results'.
#' @slot status a character vector specifying the status of the MBR object 
#' based on the available methods.
#' 
#' @section Constructor:
#' 
#' There are two constructors to create an MBR object, users can opt for one of the following: 
#' (1) \code{\link{mbrPreprocess}}; 
#' (2) \code{\link{tni2mbrPreprocess}}.
#' \itemize{
#' \item {(1)}: {It is used to create an MBR object without any pre-computed transcriptional network.}
#' \item {(2)}: {It is used to create an MBR object using available transcriptional networks.} 
#' }
#' 
#' @aliases MBR
#' @exportClass MBR

##Class MBR (Motifs Between Regulons)
setClass(
  "MBR",
  slots=c(
    TNI="TNInull",
    regulatoryElements="character",
    dualRegulons="character",
    results="list",
    para='list',
    summary='list',
    status="character"
  ), prototype=list(
    TNI=NULL,
    regulatoryElements=character(),
    dualRegulons=character(),
    results=list(),
    para=list(),
    summary=list(),
    status=character()
  )
)
