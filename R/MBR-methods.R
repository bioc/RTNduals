################################################################################
##########################         MBR-methods      ############################
################################################################################

#########################################
###class builder
#########################################
.mbr.constructor <- function(tni){
  
  #---create an MBR-object
  object <- new(Class="MBR", TNI=tni)
  
  #---status
  status <- rep('[ ]', 1, 5)
  names(status) <- c('Preprocess', 'Permutation','Bootstrap', 
                     'DPI.filter', 'Association')
  #---parameters
  sum.info.para <- list()
  sum.info.para$TNI$perm <- NA
  sum.info.para$TNI$boot <- NA
  sum.info.para$TNI$dpi <- NA
  sum.info.para$MBR$association <- as.data.frame(matrix(NA, 1, 4))
  colnames(sum.info.para$MBR$association) <- c('minRegulonSize', 'estimator',
                                               'pAdjustMethod', 'nPermutations')
  rownames(sum.info.para$MBR$association) <- 'Parameter'
  #---summary
  sum.info.summary <- list()
  sum.info.summary$MBR$Duals <- as.data.frame(matrix(NA, 1, 2))
  colnames(sum.info.summary$MBR$Duals) <- c('Tested','Predicted')
  rownames(sum.info.summary$MBR$Duals) <- 'Duals'
  
  #---set
  object <- .mbr.set(name="status", para=status, object=object)
  object <- .mbr.set(name="para", para=sum.info.para, object=object)
  object <- .mbr.set(name="summary", para=sum.info.summary, object=object)
  return(object)
}

#' A preprocessing function for objects of class MBR.
#'
#' This function converts a TNI class objects and into one MBR class object.
#'
#' @param tni A 'TNI' class object.
#' @param regulatoryElements An optional character vector specifying which 
#' 'TNI' regulatory elements should be evaluated. If 'NULL' all regulatory 
#' elements will be evaluated.
#' @return An \linkS4class{MBR} object.
#' @examples
#' ##--- load a dataset for demonstration
#' data("tniData", package = "RTN")
#' gexp <- tniData$expData
#' annot <- tniData$rowAnnotation
#' tfs <- c("IRF8","IRF1","PRDM1","E2F3","STAT4","LMO4","ZNF552")
#' 
#' ##--- construct a tni object
#' rtni <- tni.constructor(gexp, regulatoryElements = tfs, rowAnnotation=annot)
#'
#' ##--- compute regulons 
#' ## set nPermutations>=1000
#' rtni <- tni.permutation(rtni, nPermutations=30)
#' ## set nBootstrap>=100
#' rtni <- tni.bootstrap(rtni, nBootstrap=30)
#' ## 'eps=NA' estimates threshold from empirical null
#' rtni <- tni.dpi.filter(rtni, eps=NA)
#' 
#' ##--- construct a mbr object
#' rmbr <- tni2mbrPreprocess(rtni)
#'
#' @importClassesFrom RTN TNI
#' @import methods
#' @docType methods
#' @rdname tni2mbrPreprocess-methods
#' @aliases tni2mbrPreprocess
#' @export

##------------------------------------------------------------------------------
setMethod("tni2mbrPreprocess",
          "TNI",
          function (tni, regulatoryElements=NULL){
            if(tni@status["Permutation"]!="[x]")
              stop("NOTE: input 'tni' object needs permutation/bootstrep analysis!")  
            mbr.checks(name="regulatoryElements", para=regulatoryElements)
            
            #---check compatibility
            tni <- upgradeTNI(tni)
            if(tni@status["Preprocess"]!="[x]")
              stop("NOTE: TNI object is not compleate: requires preprocessing!")
            if(tni@status["Permutation"]!="[x]")
              stop("NOTE: TNI object is not compleate: requires permutation/bootstrap and DPI filter!")  
            if(tni@status["DPI.filter"]!="[x]")
              stop("NOTE: TNI object is not compleate: requires DPI filter!")
            
            #--- creates an MBR object
            object <- .mbr.constructor(tni) 
            
            #--- match regulatoryElements
            if(is.null(regulatoryElements)) {
              regulatoryElements <- tni.get(tni, "regulatoryElements")
            } else {
              regulatoryElements <- .checkRegel(tni, regulatoryElements)
            }
              
            #--- get/set status
            tni_status <- tni.get(tni, what="status")
            status <- names(tni_status[tni_status=="[x]"])
            object <- .mbr.set(name="statusUpdate", para=status, object=object)
            
            #--- get updates
            tni_summary <- tni.get(tni, what="summary")
            mbr_summary <- mbrGet(object, what="summary")
            mbr_para <- mbrGet(object, what="para")
            
            ##--- mbr_summary
            mbr_summary$TNI$tnet <- tni_summary$results$tnet
            mbr_summary$TNI$regulonSize <- tni_summary$results$regulonSize
            colnames(mbr_summary$TNI$tnet) <- c('Regulators', 'Targets', 'Edges')
            
            ##--- mbr_para
            mbr_para$perm <- tni_summary$para$perm
            ##--- bootstrap
            mbr_para$boot <- tni_summary$para$boot
            ##---summary dpi.filter
            mbr_para$dpi <- tni_summary$para$dpi
            
            #---set
            object <- .mbr.set(name="regulatoryElements", 
                               para=regulatoryElements, object=object)
            object <- .mbr.set(name="para", para=mbr_para, object=object)
            object <- .mbr.set(name="summary", para=mbr_summary, object=object)
            return (object)
            
          }
)

#' Motifs analysis and inference of 'dual regulons'.
#'
#' This function takes an MBR object and compares the shared regulon 
#' targets in order to test whether regulon pairs agree on the predicted 
#' downstream effects.
#'
#' @param object A processed object of class \linkS4class{MBR}
#' @param regulatoryElements An optional character vector specifying which 
#' 'TNI' regulatory elements should be evaluated. If 'NULL' all regulatory 
#' elements will be evaluated.
#' @param minRegulonSize A single integer or numeric value specifying the 
#' minimum number of elements in a regulon. Gene sets with fewer than this 
#' number are removed from the analysis.
#' @param doSizeFilter a logical value. If TRUE, negative and positive targets are 
#' independently verified by the 'minRegulonSize' argument.
#' @param pValueCutoff a single numeric value specifying the cutoff for p-values 
#' considered significant.
#' @param pAdjustMethod A single character value specifying the p-value 
#' adjustment method to be used (see 'p.adjust' function for details).
#' @param estimator A character value specifying the estimator used in the 
#' association analysis. One of "spearman" (default), "kendall", or "pearson".
#' @param nPermutations A single integer value specifying the number of 
#' permutations for deriving p-values associating regulon pairs.
#' @param miFilter A single logical value specifying to apply the 'miFilter' 
#' between two regulators.
#' @param verbose A single logical value specifying to display detailed 
#' messages (when verbose=TRUE) or not (when verbose=FALSE).
#' @return An \linkS4class{MBR} object with two data.frames in the slot 
#' 'results' listing the inferred 'dual regulons' and correspoding statistics.
#' @examples
#' ##--- load a dataset for demonstration
#' data("tniData", package = "RTN")
#' gexp <- tniData$expData
#' annot <- tniData$rowAnnotation
#' tfs <- c("IRF8","IRF1","PRDM1","E2F3","STAT4","LMO4","ZNF552")
#' 
#' ##--- construct a tni object
#' rtni <- tni.constructor(gexp, regulatoryElements = tfs, rowAnnotation=annot)
#'
#' ##--- compute regulons 
#' ## set nPermutations>=1000
#' rtni <- tni.permutation(rtni, nPermutations=30)
#' ## set nBootstrap>=100
#' rtni <- tni.bootstrap(rtni, nBootstrap=30)
#' ## 'eps=NA' estimates threshold from empirical null
#' rtni <- tni.dpi.filter(rtni, eps=NA)
#' 
#' ##--- construct a mbr object
#' rmbr <- tni2mbrPreprocess(rtni)
#' 
#' ##--- run mbrAssociation 
#' ## set nPermutations>=1000
#' rmbr <- mbrAssociation(rmbr, pValueCutoff = 0.05, nPermutations=30)
#'
#' @import RTN 
#' @importFrom stats p.adjust phyper median pnorm sd
#' @importFrom stats cor quantile pt complete.cases
#' @import methods
#' @docType methods
#' @rdname mbrAssociation-methods
#' @aliases mbrAssociation
#' @export

##------------------------------------------------------------------------------
##Inference of duals
setMethod("mbrAssociation",
          "MBR",
          function(object, regulatoryElements=NULL, minRegulonSize=15, 
                   doSizeFilter=FALSE,  pValueCutoff=0.001, 
                   pAdjustMethod="bonferroni", estimator="spearman", 
                   nPermutations=1000, miFilter=TRUE, verbose=TRUE){
            ##--- input check
            if(object@status["Preprocess"]!="[x]")
              stop("NOTE: MBR object is not complete: requires preprocessing!")
            if(object@status["Permutation"]!="[x]")
              stop("NOTE: MBR object is not complete: requires permutation/bootstrap and DPI filter!")  
            if(object@status["DPI.filter"]!="[x]")
              stop("NOTE: MBR object is not complete: requires DPI filter!")
            ##--- gets
            tni <- mbrGet(object, what="TNI")
            tni_gexp <- tni.get(tni, what="gexp")
            tni_para <- tni.get(tni, what="para")
            ##--- checks
            mbr.checks(name="regulatoryElements", para=regulatoryElements)
            mbr.checks(name="minRegulonSize", para=minRegulonSize)
            mbr.checks(name="doSizeFilter",para=doSizeFilter)
            mbr.checks(name="estimator", para=estimator)
            mbr.checks(name="pValueCutoff", para=pValueCutoff)
            mbr.checks(name="pAdjustMethod", para=pAdjustMethod)
            mbr.checks(name="nPermutations", para=nPermutations)
            mbr.checks(name="miFilter", para=miFilter)
            mbr.checks(name="verbose", para=verbose)
            
            if(verbose) cat("-Checking regulons and regulatory elements...\n")
            
            if(is.null(regulatoryElements)) {
              regulatoryElements <- mbrGet(object, "regulatoryElements")
            } else {
              temp <- mbrGet(object, "regulatoryElements")
              regulatoryElements <- temp[temp%in%regulatoryElements|
                                           temp%in%names(regulatoryElements)]
            }
            mbr.checks(name="numberRegElements", para=regulatoryElements)
            
            ##--- get all regulons (from tni)
            regulonsAndMode <- tni.get(tni, what="regulons.and.mode")
            regulons <- tni.get(tni, what="regulons")
            
            ##--- get regulons for 'regulatoryElements'
            regulonsAndMode <- regulonsAndMode[regulatoryElements]
            regulons <- regulons[regulatoryElements]
            
            ##-----check regulon size (both clouds)
            gs.size.max <- unlist(lapply(regulonsAndMode, function(reg){
              max(sum(reg>0),sum(reg<0))
            }))
            gs.size.min <- unlist(lapply(regulonsAndMode, function(reg){
              min(sum(reg>0),sum(reg<0))
            }))
            ##-----stop when no subset passes the size requirement
            if(all(gs.size.max<minRegulonSize)){
              stop(paste("NOTE: no partial regulon has minimum >= ", minRegulonSize, sep=""))
            }
            ##-----get filtered list
            if(doSizeFilter){
              regulatoryElements<-regulatoryElements[which(gs.size.min>=minRegulonSize)]
              if(length(regulatoryElements)<2){
                stop("NOTE: at least two regulons have to pass the 'doSizeFilter' requirement!")
              }
            } else {
              regulatoryElements<-regulatoryElements[which(gs.size.max>=minRegulonSize)]
              if(length(regulatoryElements)<2){
                stop("NOTE: at least two regulons have to pass the 'minRegulonSize' requirement!")
              }
            }
            ##---
            regulonsAndMode <- regulonsAndMode[regulatoryElements]
            regulons <- regulons[regulatoryElements]
            
            ##--- group of regulons and regulatory elements
            targets <- unique(c(unlist(regulons)))
            targets <- unique(c(regulatoryElements, targets))
            
            ##--- get mi from refnet
            if(verbose) {
                cat("-Extrating inferred regulatory associations...\n")
            }
            mitnet <- abs(tni.get(tni, what="refnet"))
            mitnet <- mitnet[targets,regulatoryElements, drop=FALSE]
            ntargets <- tni.get(tni,"summary")$results$tnet["tnet.ref","Targets"]
            
            if(verbose) {
              cat("-Mapping network triplets between regulons...\n")
            }
            
            ##--- compute correlation between regulators and all targets
            cortnet <- .tni.cor(tni_gexp, mitnet, estimator=estimator, dg=0, 
                                asInteger=FALSE, mapAssignedAssociation=FALSE)
            
            ##--- map interactions for all potential dual regulons
            uniregulons <- .getUnionRegulons(regulons)
            
            ##--- compute correlation between regulons
            ## the transformation with double 'cor' calls combines the two distributions
            ## into a single one, making it more symmetrical, while preserveing the 
            ## directionality of the approach.
            regcor <- sapply(regulatoryElements, function(r1){
              sapply(regulatoryElements, function(r2){
                tar12 <- uniregulons[[r1]][[r2]]
                r1 <- cortnet[tar12,r1]
                r2 <- cortnet[tar12,r2]
                c1 <- cor(r1, r2, method=estimator)
                c2 <- cor(abs(r1), abs(r2), method=estimator)
                abs(c2)*sign(c1)
              })
            })
            regcor[is.na(regcor)] <- 0
            regcor <- t(regcor)
            
            #--- 'cormat' and 'sigmat' will be returned at the end 
            cormat <- regcor
            sigmat <- cormat
            sigmat[,] <- ""
            
            #--- set NAs for regs eventually included in both lists
            diag(regcor) <- NA
            regcor[upper.tri(regcor)] <- NA
            testedDuals <- sum(!is.na(regcor))

            ##--- get statlist
            statlist <- .getStatlist(regcor=regcor, regulatoryElements)
            statlist <- .getMI(statlist, mitnet)
            if(miFilter) statlist <- statlist[which(statlist$MI != 0), ]

            ##--- further assessing inferred duals
            if(nrow(statlist) > 0){
              
              ##--- n. tests before any filter
              n.tests <- sum(!is.na(regcor))
              
              ##--- assessing overlap between regulons
              if(verbose) {
                tp <- paste("-Assessing overlap between",
                            length(regulons),"regulons...\n")
                cat(tp)
              }
              overlapStats <- .mbr.overlap(statlist, regulons, ntargets, verbose)
              ##--- adjust Pvalue for n.tests
              overlapStats$Adjusted.Pvalue <- p.adjust(overlapStats$Pvalue, 
                                                       method=pAdjustMethod, 
                                                       n=n.tests)
              overlapStats <- overlapStats[overlapStats$Adjusted.Pvalue<pValueCutoff,]
              statlist <- statlist[rownames(overlapStats),]
              
              ##--- running permutation on R.Regulons
              if(verbose) {
                tp <- paste("-Assessing correlation between",
                            length(regulons),"regulons...\n")
                cat(tp)
              }
              
              if(nrow(statlist) > 0){
                
                ##--- Permutation on correlation
                corrStats <- .permCorPval(cortnet, statlist, uniregulons, 
                                          estimator, nper=nPermutations, 
                                          verbose=verbose)
                ##--- adjust Pvalue for n.tests
                corrStats$Adjusted.Pvalue <- p.adjust(corrStats$Pvalue, 
                                                      method=pAdjustMethod, n=n.tests)
                corrStats <- corrStats[corrStats$Adjusted.Pvalue<pValueCutoff,]
                statlist <- statlist[rownames(corrStats),]
                
                #--- align results
                overlapStats <- overlapStats[rownames(statlist),]
                corrStats <- corrStats[rownames(statlist),]
                
                ##--- map significant duals in 'sigmat' 
                if(nrow(statlist)>0){
                  for(i in 1:nrow(statlist)){
                    r1 <- statlist$Regulon1[i]
                    r2 <- statlist$Regulon2[i]
                    sigmat[r1,r2] <- "*"
                    if( r1%in%colnames(sigmat) && r2%in%rownames(sigmat)){
                      sigmat[r2,r1] <- "*"
                    }
                  }
                }
              }
            }
            
            ##--- organize results
            predictedDuals <- nrow(statlist)
            if(predictedDuals>0){
              labs <- c("Regulon1","Regulon2", "Universe.Size", "Regulon1.Size", "Regulon2.Size",
                        "Expected.Overlap","Observed.Overlap","Pvalue", "Adjusted.Pvalue")
              overlapStats <- overlapStats[,labs]
              labs <- c("Regulon1","Regulon2", "MI.Regulators", "R.Regulons", 
                        "Pvalue", "Adjusted.Pvalue")
              corrStats <- corrStats[,labs]
            } else {
              warning("No 'dual regulon' has been observed for the input parameters.",
                      call.=FALSE)
              nm <- as.numeric()
              ch <- as.numeric()
              overlapStats <- data.frame(Regulon1=ch, egulon2=ch, Universe.Size=nm, Regulon1.Size=nm, 
                                         Regulon2.Size=nm, Expected.Overlap=nm, Observed.Overlap=nm,
                                         Pvalue=nm, Adjusted.Pvalue=nm)
              corrStats <- data.frame(Regulon1=ch, Regulon2=ch, MI.Regulators=nm, 
                                      R.Regulons=nm, Pvalue=nm, Adjusted.Pvalue=nm)
            }
            
            ##--- para
            mbr_para <- mbrGet(object,what="para")
            sum.info.par <- data.frame(minRegulonSize, estimator, pAdjustMethod, nPermutations, 
                                       stringsAsFactors = FALSE)
            mbr_para$MBR$association['Parameter', ] <- sum.info.par
            
            ##---
            mbr_summary <- mbrGet(object, what="summary")
            mbr_summary$MBR$Duals[,'Tested'] <- testedDuals
            mbr_summary$MBR$Duals[,'Predicted'] <- predictedDuals
            ##--- set
            object <- .mbr.set(name="statusUpdate", 
                               para="Association", object=object)
            object <- .mbr.set(name="para", 
                               para=mbr_para, object=object)
            object <- .mbr.set(name="summary", 
                               para=mbr_summary, object=object)
            object <- .mbr.set(name="regulatoryElements", 
                               para=regulatoryElements, object=object)
            object <- .mbr.set(name="dualRegulons", 
                               para=rownames(statlist), object=object)
            object <- .mbr.set(name="dualsCorrelation", 
                               para=corrStats, object=object)
            object <- .mbr.set(name="dualsOverlap", 
                               para=overlapStats, object=object)
            object@results$dualsmat <- list(cormat=cormat, sigmat=sigmat)
            return(object)
          }
)


#' Entry point for external evidences.
#'
#' If available, this function adds external evidences to an 'MBR' object.
#'
#' @param object A processed object of class \linkS4class{MBR} evaluated by the 
#' method \code{\link{mbrAssociation}}.
#' @param priorEvidenceTable An 'data.frame' with three columns 
#' representing (1) regulatory elements 1, (2) regulatory elements 2,
#' and (3) external evidences between the regulatory elements.
#' @param evidenceColname A single character value specifying a column in 
#' the 'priorEvidenceTable'.
#' @param verbose A single logical value specifying to display detailed 
#' messages (when verbose=TRUE) or not (when verbose=FALSE).
#' @return An \linkS4class{MBR} object with an updated 'data.frame' in the slot 
#' 'results' listing the input additional evidences.
#' @examples
#' ##--- load a dataset for demonstration
#' data("tniData", package = "RTN")
#' gexp <- tniData$expData
#' annot <- tniData$rowAnnotation
#' tfs <- c("IRF8","IRF1","PRDM1","E2F3","STAT4","LMO4","ZNF552")
#' 
#' ##--- construct a tni object
#' rtni <- tni.constructor(gexp, regulatoryElements = tfs, rowAnnotation=annot)
#'
#' ##--- compute regulons 
#' ## set nPermutations>=1000
#' rtni <- tni.permutation(rtni, nPermutations=30)
#' ## set nBootstrap>=100
#' rtni <- tni.bootstrap(rtni, nBootstrap=30)
#' ## 'eps=NA' estimates threshold from empirical null
#' rtni <- tni.dpi.filter(rtni, eps=NA)
#' 
#' ##--- construct a mbr object
#' rmbr <- tni2mbrPreprocess(rtni)
#' 
#' ##--- run mbrAssociation 
#' ## set nPermutations>=1000
#' rmbr <- mbrAssociation(rmbr, pValueCutoff = 0.05, nPermutations=30)
#'
#' ##--- check results
#' results <- mbrGet(rmbr, what="dualsCorrelation")
#' 
#' ##--- add supplementary evidence table
#' ## here we build a 'toy' example using the 'rnorm' function 
#' ## for demonstration purposes only!
#' priorEvidenceTable <- results[,c("Regulon1","Regulon2")]
#' priorEvidenceTable$ToyEvidence <- rnorm(nrow(results))
#' priorEvidenceTable
#' 
#' ##--- add supplementary evidences
#' # rmbr <- mbrPriorEvidenceTable(rmbr, priorEvidenceTable=priorEvidenceTable, evidenceColname = "ToyEvidence")
#' 
#' ##--- check updated results
#' # mbrGet(rmbr, what="dualsCorrelation")
#'
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @import methods
#' @docType methods
#' @rdname mbrPriorEvidenceTable-methods
#' @aliases mbrPriorEvidenceTable
#' @export

##------------------------------------------------------------------------------
##organize duals
setMethod( "mbrPriorEvidenceTable",
           "MBR",
           function(object, priorEvidenceTable, evidenceColname, 
                    verbose=TRUE){
             if(object@status["Association"]!="[x]")
               stop("NOTE: MBR object is not compleate: requires 'mbrAssociation' analysis!")
             ##--- initial checks
             mbr.checks(name="object", para=object)
             dualsCorrelation <- mbrGet(object, what="dualsCorrelation")
             if(is.null(dim(dualsCorrelation)))
               stop("NOTE: empty results in the input 'object'!",call.=FALSE)
             
             #--- update object
             object <- .mbr.set(name="dualsCorrelation", 
                                para=dualsCorrelation, object=object)
             object <- .mbr.set(name="dualRegulons", 
                                para=rownames(dualsCorrelation), object=object)
             
             #--- add priorEvidenceTable if avaible
             ##--- general checks
             if(missing(evidenceColname)){
               stop("'evidenceColname' should be listed in colnames of 'priorEvidenceTable'!",call.=FALSE)
             }
             priorEvidenceTable <- mbr.checks(name="priorEvidenceTable", 
                                              para=priorEvidenceTable, 
                                              paraSuppl=evidenceColname)
             ##--- check consistency
             if(verbose) cat("-Checking annotation consistency in 'priorEvidenceTable'...\n")
             .checkConsistencySuppTable(object, priorEvidenceTable, verbose=verbose)
             ##--- update duals
             object <- .updateEvidenceTable(object, priorEvidenceTable, verbose=verbose)
             return(object)
           }
)

##------------------------------------------------------------------------------
##show summary information on screen
setMethod( "show",
           "MBR",
           function(object){
             cat("An MBR (Motifs Between Regulons) object:\n")
             message("--status:")
             print(object@status, quote=FALSE)
           }
)

#' Get information from individual slots in MBR object.
#' 
#' Get information from individual slots in an MBR object and any available 
#' results from previous analysis.
#' 
#' @param object A preprocessed object of class \linkS4class{MBR}
#' @param what a single character value specifying which information should be 
#' retrieved from the slots. Options: "TNI", "regulatoryElements", 
#' "dualRegulons", "results", "para", "summary", 
#' "status", "dualsCorrelation", "dualsOverlap", and "dualsCorMatrix"
#' @return Content from slots in the \linkS4class{MBR} object
#' @examples
#' ##--- load a dataset for demonstration
#' data("tniData", package = "RTN")
#' gexp <- tniData$expData
#' annot <- tniData$rowAnnotation
#' tfs <- c("IRF8","IRF1","PRDM1","E2F3","STAT4","LMO4","ZNF552")
#' 
#' ##--- construct a tni object
#' rtni <- tni.constructor(gexp, regulatoryElements = tfs, rowAnnotation=annot)
#'
#' ##--- compute regulons 
#' ## set nPermutations>=1000
#' rtni <- tni.permutation(rtni, nPermutations=30)
#' ## set nBootstrap>=100
#' rtni <- tni.bootstrap(rtni, nBootstrap=30)
#' ## 'eps=NA' estimates threshold from empirical null
#' rtni <- tni.dpi.filter(rtni, eps=NA)
#' 
#' ##--- construct a mbr object
#' rmbr <- tni2mbrPreprocess(rtni)
#'
#' ##--- get the 'TNI' slot using 'mbrGet'
#' tni <- mbrGet(rmbr, what="TNI")
#' 
#' @import methods
#' @docType methods
#' @rdname mbrGet-methods
#' @aliases mbrGet
#' @export
##------------------------------------------------------------------------------
##get slots from MBR object
setMethod( "mbrGet",
           "MBR", 
           function(object, what="status"){
             ##---check input arguments
             mbr.checks(name="mbrGet", para=what)
             ##---options that need 'mbrAssociation' evaluation!
             optsAssoci <- c("regulatoryElements", 
                             "results", "dualRegulons", "dualsCorrelation",
                             "dualsOverlap","dualsCorMatrix")
             ##---get query
             if(what=="TNI"){
               query <- object@TNI
             } else if(what=="para"){
               query <- object@para
             } else if(what=="summary"){
               query <- object@summary
             } else if(what=="status"){
               query <- object@status
             } else if(what=="regulatoryElements"){
               query <- object@regulatoryElements               
             } else if(what%in%optsAssoci){
               if(object@status["Association"] != "[x]"){
                 warning("Input 'object' needs 'mbrAssociation' evaluation!", call.=FALSE)
                 query <- NULL
               } else {
                 if(what=="dualRegulons"){
                   query <- object@dualRegulons
                 } else if(what=="results"){
                   query <- object@results
                 } else if(what=="dualsCorrelation"){
                   query <- object@results$dualsCorrelation
                   query <- query[sort.list(query$Pvalue),]
                 } else if(what=="dualsOverlap"){
                   query <- object@results$dualsOverlap
                   query <- query[sort.list(query$Pvalue),]
                 } else if(what=="dualsCorMatrix"){
                   query <- object@results$dualsmat
                 }
               }
             }
             return(query)
           }
)
