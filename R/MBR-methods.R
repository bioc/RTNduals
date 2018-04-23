################################################################################
##########################         MBR-methods      ############################
################################################################################

#########################################
###class builder
#########################################
.MBRmaker <- function(gexp, regulatoryElements1, regulatoryElements2){
  #---creates the 'MBR' object
  object <- newMBR(gexp=gexp,
                   regulatoryElements1=regulatoryElements1,
                   regulatoryElements2=regulatoryElements2
  )
  #---status
  status <- rep('[ ]', 1, 5)
  names(status) <- c('Preprocess', 'Permutation','Bootstrap', 
                     'DPI.filter', 'Association')
  #---parameters
  sum.info.para <- list()
  sum.info.para$TNIs$perm <- NA
  sum.info.para$TNIs$boot <- NA
  sum.info.para$TNIs$dpi <- NA
  sum.info.para$MBR$association <- as.data.frame(matrix(NA, 1, 4))
  colnames(sum.info.para$MBR$association) <- c('minRegulonSize', 'estimator',
                                               'pAdjustMethod', 'nPermutations')
  rownames(sum.info.para$MBR$association) <- 'Parameter'
  #---summary
  sum.info.summary <- list()
  sum.info.summary$MBR$Duals <- as.data.frame(matrix(NA, 1))
  colnames(sum.info.summary$MBR$Duals) <- c('testedDuals')
  rownames(sum.info.summary$MBR$Duals) <- 'duals'
  
  #---set
  object <- .mbr.set(name="status", para=status, object=object)
  object <- .mbr.set(name="para", para=sum.info.para, object=object)
  object <- .mbr.set(name="summary", para=sum.info.summary, object=object)
  return(object)
}

#----------------------------------------------------------
#it creates the 'MBR' class object
newMBR <- function(gexp, regulatoryElements1, regulatoryElements2){
  #---checks
  if(missing(gexp)) 
    stop("NOTE: 'gexp' is missing ", call.=FALSE)
  if(missing(regulatoryElements1)) 
    stop("NOTE: 'regulatoryElements1' is missing", call.=FALSE)
  if(missing(regulatoryElements2)) 
    stop("NOTE: 'regulatoryElements2' is missing", call.=FALSE)
  mbr.checks(name='gexp', gexp)
  mbr.checks(name='regulatoryElements1', regulatoryElements1)
  mbr.checks(name='regulatoryElements2', regulatoryElements2)
  
  #---creating TNIs
  regulonsTNI1 <- new("TNI", gexp=gexp, 
                      regulatoryElements=regulatoryElements1)
  regulonsTNI2 <- new("TNI", gexp=gexp, 
                      regulatoryElements=regulatoryElements2)
  
  #---creating MBR-object
  new(Class="MBR", TNI1=regulonsTNI1, TNI2=regulonsTNI2)
}

##------------------------------------------------------------------------------
#' A preprocessing function for objects of class MBR.
#'
#' @param gexp A numerical matrix, typically with mRNA and/or miRNA expression 
#' values.
#' @param regulatoryElements1 A named vector with regulatory elements listed in 
#' 'gexp' rownames.
#' @param regulatoryElements2 A named vector with regulatory elements listed in 
#' 'gexp' rownames.
#' @param verbose A single logical value specifying to display detailed 
#' messages (when verbose=TRUE) or not (when verbose=FALSE).
#' @param ... Additional arguments passed on to 
#' \code{\link[RTN:tni.preprocess]{tni.preprocess}} function.
#' @return A preprocessed 'MBR-class' object.
#' @examples
#' ##--- load a dataset for demonstration
#' data("dt4rtn", package = "RTN")
#' gexp <- dt4rtn$gexp
#' annot <- dt4rtn$gexpIDs
#' tfs1 <- dt4rtn$tfs[c("IRF8","IRF1","PRDM1","AFF3","E2F3")]
#' tfs2 <- dt4rtn$tfs[c("HCLS1","STAT4","STAT1","LMO4","ZNF552")]
#' 
#' ##--- run mbrPreprocess
#' rmbr <- mbrPreprocess(gexp=gexp, regulatoryElements1 = tfs1, 
#' regulatoryElements2=tfs2, rowAnnotation=annot)
#'
#' @import RTN
#' @import methods
#' @docType methods
#' @rdname mbrPreprocess-methods
#' @aliases mbrPreprocess
#' @export

##Regulons pre-processing method
setMethod("mbrPreprocess",
          "matrix",
          function(gexp, regulatoryElements1, regulatoryElements2, 
                   verbose=TRUE,...){
            ##---
            mbr.checks(name="verbose", para=verbose)  
            object <- .MBRmaker(gexp=gexp, 
                                regulatoryElements1=regulatoryElements1,
                                regulatoryElements2=regulatoryElements2)  
            ##---get TNIs
            TNI1 <- mbrGet(object, what="TNI1")
            TNI2 <- mbrGet(object, what="TNI2")
            ##---pre-processing TNIs
            if(verbose) cat("-Preprocessing TNI objects...\n\n")
            TNI1 <- tni.preprocess(TNI1, verbose=verbose,...=...)
            TNI2 <- tni.preprocess(TNI2, verbose=verbose,...=...)
            regulatoryElements1 <- tni.get(TNI1, what="regulatoryElements")
            regulatoryElements2 <- tni.get(TNI2, what="regulatoryElements")
            mbr.checks(name="regulatoryElements", para=regulatoryElements1)
            mbr.checks(name="regulatoryElements", para=regulatoryElements2)
            
            ##---set
            object <- .mbr.set(name="TNI1", para=TNI1, object=object)
            object <- .mbr.set(name="TNI2", para=TNI2, object=object)
            object <- .mbr.set(name="statusUpdate", para="Preprocess", 
                               object=object)
            
            return(object)
          }
)
##------------------------------------------------------------------------------
#' Inference of transcriptional networks.
#'
#' This function takes an MBR object and computes two transcriptional networks 
#' inferred 
#' by mutual information (with multiple hypothesis testing corrections).
#' 
#' @param object A preprocessed object of class \linkS4class{MBR}.
#' @param verbose A single logical value specifying to display detailed 
#' messages (when verbose=TRUE) or not (when verbose=FALSE).
#' @param ... Additional arguments passed on to the 
#' \code{\link[RTN:tni.permutation]{tni.permutation}} function.
#' @return An \linkS4class{MBR} object with two mutual information matrices, 
#' one in each "TNI" slot.
#' @examples
#' ##--- load a dataset for demonstration
#' data("dt4rtn", package = "RTN")
#' gexp <- dt4rtn$gexp
#' annot <- dt4rtn$gexpIDs
#' tfs1 <- dt4rtn$tfs[c("IRF8","IRF1","PRDM1","AFF3","E2F3")]
#' tfs2 <- dt4rtn$tfs[c("HCLS1","STAT4","STAT1","LMO4","ZNF552")]
#' 
#' ##--- run mbrPreprocess
#' rmbr <- mbrPreprocess(gexp=gexp, regulatoryElements1 = tfs1, 
#' regulatoryElements2=tfs2, rowAnnotation=annot)
#' 
#' ##--- run mbrPermutation (set nPermutations>=1000)
#' rmbr <- mbrPermutation(rmbr, nPermutations=10)
#'
#' @import RTN 
#'
#' @import methods
#' @docType methods
#' @rdname mbrPermutation-methods
#' @aliases mbrPermutation
#' @export

## permutation
setMethod("mbrPermutation",
          "MBR",
          function(object, verbose=TRUE, ...){
            ##--- input check
            if(object@status["Preprocess"]!="[x]")
              stop("NOTE: MBR object is not compleate: requires preprocessing!")
            ##---checks
            mbr.checks(name="object", para=object)
            mbr.checks(name="verbose", para=verbose)
            ##---get TNIs
            TNI1 <- mbrGet(object, what="TNI1")
            TNI2 <- mbrGet(object, what="TNI2")
            
            ##---permutation TNIs
            if(verbose)
              cat("-Performing permutation analysis for two TNI objects...\n\n")
            TNI1 <- tni.permutation(TNI1, verbose=verbose,...=...)
            TNI2 <- tni.permutation(TNI2, verbose=verbose,...=...)
            #---get
            tni1_summary <- tni.get(TNI1, what="summary")
            tni2_summary <- tni.get(TNI2, what="summary")
            mbr_summary <- mbrGet(object, what="summary")
            mbr_para <- mbrGet(object, what="para")
            
            #---changes
            mbr_para$TNIs$perm <- tni1_summary$para$perm
            mbr_summary$TNIs$TNI1 <- tni1_summary$results$tnet
            colnames(mbr_summary$TNIs$TNI1) <- c('RE', 'Targets', 'Edges')
            mbr_summary$TNIs$TNI2 <- tni2_summary$results$tnet
            colnames(mbr_summary$TNIs$TNI2) <- c('RE', 'Targets', 'Edges')
            
            #---set
            object <- .mbr.set(name="TNI1", para=TNI1, object=object)
            object <- .mbr.set(name="TNI2", para=TNI2, object=object)
            object <- .mbr.set(name="statusUpdate", para="Permutation", 
                               object=object)
            object <- .mbr.set(name="para", para=mbr_para, object=object)
            object <- .mbr.set(name="summary", para=mbr_summary, object=object)
            return(object)
          }
)

#' Inference of consensus transcriptional networks.
#'
#' This function takes an MBR object and computes two consensus transcriptional 
#' networks.
#'
#' @param object A processed objec of class \linkS4class{MBR} evaluated by the 
#' method \code{\link[RTNduals:mbrPermutation]{mbrPermutation}}.
#' @param verbose A single logical value specifying to display detailed 
#' messages (when verbose=TRUE) or not (when verbose=FALSE).
#' @param ... Additional arguments passed to the 
#' \code{\link[RTN:tni.bootstrap]{tni.bootstrap}} function.
#' @return An \linkS4class{MBR} object with two consensus mutual information 
#' matrices, one in each "TNI" slot.
#' @examples
#' ##--- load a dataset for demonstration
#' data("dt4rtn", package = "RTN")
#' gexp <- dt4rtn$gexp
#' annot <- dt4rtn$gexpIDs
#' tfs1 <- dt4rtn$tfs[c("IRF8","IRF1","PRDM1","AFF3","E2F3")]
#' tfs2 <- dt4rtn$tfs[c("HCLS1","STAT4","STAT1","LMO4","ZNF552")]
#' 
#' ##--- run mbrPreprocess
#' rmbr <- mbrPreprocess(gexp=gexp, regulatoryElements1 = tfs1, 
#' regulatoryElements2=tfs2, rowAnnotation=annot)
#' 
#' ##--- run mbrPermutation (nPermutations>=1000)
#' rmbr <- mbrPermutation(rmbr, nPermutations=10)
#' 
#' ##--- run mbrBootstrap (nPermutations>=100)
#' rmbr <- mbrBootstrap(rmbr, nBootstrap=10)
#'
#' @import RTN 
#'
#' @import methods
#' @docType methods
#' @rdname mbrBootstrap-methods
#' @aliases mbrBootstrap
#' @export

##------------------------------------------------------------------------------
## bootstrap method
setMethod("mbrBootstrap",
          "MBR",
          function(object, verbose=TRUE, ...){
            ##--- input check
            if(object@status["Preprocess"]!="[x]")
              stop("NOTE: MBR object is not compleate: requires preprocessing!")
            if(object@status["Permutation"]!="[x]")
              stop("NOTE: MBR object is not compleate: requires permutation analysis!") 
            ##---checks
            mbr.checks(name="object", para=object)
            mbr.checks(name="verbose", para=verbose)
            ##---get TNIs
            TNI1 <- mbrGet(object, what="TNI1")
            TNI2 <- mbrGet(object, what="TNI2")
            
            ##---bootstrap TNIs
            if(verbose) cat("-Performing bootstrap analysis for two TNI objects...\n\n")
            TNI1 <- tni.bootstrap(TNI1, verbose=verbose,...=...)
            TNI2 <- tni.bootstrap(TNI2, verbose=verbose,...=...)
            
            #---get
            tni1_summary <- tni.get(TNI1, what="summary")
            tni2_summary <- tni.get(TNI2, what="summary")
            mbr_summary <- mbrGet(object, what="summary")
            mbr_para <- mbrGet(object, what="para")
            
            #---changes
            mbr_para$TNIs$boot <- tni1_summary$para$boot
            
            mbr_summary$TNIs$TNI1 <- tni1_summary$results$tnet
            colnames(mbr_summary$TNIs$TNI1) <- c('RE', 'Targets', 'Edges')
            mbr_summary$TNIs$TNI2 <- tni2_summary$results$tnet
            colnames(mbr_summary$TNIs$TNI2) <- c('RE', 'Targets', 'Edges')
            
            #---set
            object <- .mbr.set(name="TNI1", para=TNI1, object=object)
            object <- .mbr.set(name="TNI2", para=TNI2, object=object)
            object <- .mbr.set(name="statusUpdate", para="Bootstrap", 
                               object=object)
            object <- .mbr.set(name="para", para=mbr_para, object=object)
            object <- .mbr.set(name="summary", para=mbr_summary, object=object)
            return(object)
          }
)

#' A filter based on the Data Processing Inequality (DPI) algorithm.
#'
#' This function takes an MBR object and computes two transcriptional networks 
#' filtered by the data processing inequality algorithm.
#'
#' @param object A processed object of class \linkS4class{MBR} evaluated by 
#' the methods
#'  \code{\link[RTNduals:mbrPermutation]{mbrPermutation}} and 
#'  \code{\link[RTNduals:mbrBootstrap]{mbrBootstrap}}.
#' @param verbose A single logical value specifying to display detailed 
#' messages (when verbose=TRUE) or not (when verbose=FALSE).
#' @param ... Additional arguments passed to the 
#' \code{\link[RTN:tni.dpi.filter]{tni.dpi.filter}} function.
#' @return An \linkS4class{MBR} object with two DPI-filtered mutual information 
#' matrices, one in each "TNI" slot.
#' @examples
#' ##--- load a dataset for demonstration
#' data("dt4rtn", package = "RTN")
#' gexp <- dt4rtn$gexp
#' annot <- dt4rtn$gexpIDs
#' tfs1 <- dt4rtn$tfs[c("IRF8","IRF1","PRDM1","AFF3","E2F3")]
#' tfs2 <- dt4rtn$tfs[c("HCLS1","STAT4","STAT1","LMO4","ZNF552")]
#' 
#' ##--- run mbrPreprocess
#' rmbr <- mbrPreprocess(gexp=gexp, regulatoryElements1 = tfs1, 
#' regulatoryElements2=tfs2, rowAnnotation=annot)
#' 
#' ##--- run mbrPermutation (set nPermutations>=1000)
#' rmbr <- mbrPermutation(rmbr, nPermutations=10)
#' 
#' ##--- run mbrBootstrap (set nBootstrap>=100)
#' rmbr <- mbrBootstrap(rmbr, nBootstrap=10)
#' 
#' ##--- run mbrDpiFilter
#' rmbr <- mbrDpiFilter(rmbr)
#'
#' @import RTN 
#'
#' @import methods
#' @docType methods
#' @rdname mbrDpiFilter-methods
#' @aliases mbrDpiFilter
#' @export

##------------------------------------------------------------------------------
## dpi filter method
setMethod("mbrDpiFilter",
          "MBR",
          function(object, verbose=TRUE, ...){
            ##--- input check
            if(object@status["Preprocess"]!="[x]")
              stop("NOTE: MBR object is not compleate: requires preprocessing!")
            if(object@status["Permutation"]!="[x]")
              stop("NOTE: MBR object is not compleate: requires permutation/bootstrap analysis!")
            ##---checks
            mbr.checks(name="object", para=object)
            mbr.checks(name="verbose", para=verbose)
            ##---get TNIs
            TNI1 <- mbrGet(object, what="TNI1")
            TNI2 <- mbrGet(object, what="TNI2")
            
            ##---Dpi filter TNIs
            if(verbose) cat("-Applying dpi filter for two TNI objects...\n\n")
            TNI1 <-tni.dpi.filter(TNI1, verbose=verbose, ...=...)
            TNI2 <-tni.dpi.filter(TNI2, verbose=verbose, ...=...)
            #---get
            tni1_summary <- tni.get(TNI1, what="summary")
            tni2_summary <- tni.get(TNI2, what="summary")
            mbr_summary <- mbrGet(object, what="summary")
            mbr_para <- mbrGet(object, what="para")
            
            #---changes
            mbr_para$TNIs$dpi <- tni1_summary$para$dpi
            
            mbr_summary$TNIs$TNI1 <- tni1_summary$results$tnet
            colnames(mbr_summary$TNIs$TNI1) <- c('RE', 'Targets', 'Edges')
            mbr_summary$TNIs$TNI2 <- tni2_summary$results$tnet
            colnames(mbr_summary$TNIs$TNI2) <- c('RE', 'Targets', 'Edges')
            
            #---set
            object <- .mbr.set(name="TNI1", para=TNI1, object=object)
            object <- .mbr.set(name="TNI2", para=TNI2, object=object)
            object <- .mbr.set(name="statusUpdate", para="DPI.filter", 
                               object=object)
            object <- .mbr.set(name="para", para=mbr_para, object=object)
            object <- .mbr.set(name="summary", para=mbr_summary, object=object)
            return(object)
          }
)

#' Motifs analysis and inference of 'dual regulons'.
#'
#' This function takes an MBR object and compares the shared regulon 
#' targets in order to test whether regulon pairs agree on the predicted 
#' downstream effects.
#'
#' @param object A processed object of class \linkS4class{MBR} evaluated by the 
#' methods \code{\link[RTNduals:mbrPermutation]{mbrPermutation}}, 
#' \code{\link[RTNduals:mbrBootstrap]{mbrBootstrap}} 
#' and \code{\link[RTNduals:mbrDpiFilter]{mbrDpiFilter}}.
#' @param regulatoryElements1 An optional character vector specifying which 
#' 'TNI1' regulatory elements should be evaluated. If 'NULL' all regulatory 
#' elements will be evaluated.
#' @param regulatoryElements2 An optional character vector specifying which 
#' 'TNI2' regulatory elements should be evaluated. If 'NULL' all regulatory 
#' elements will be evaluated.
#' @param minRegulonSize A single integer or numeric value specifying the 
#' minimum number of elements in a regulon. Gene sets with fewer than this 
#' number are removed from the analysis.
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
#' data("dt4rtn", package = "RTN")
#' gexp <- dt4rtn$gexp
#' annot <- dt4rtn$gexpIDs
#' tfs1 <- dt4rtn$tfs[c("IRF8","IRF1","PRDM1","AFF3","E2F3")]
#' tfs2 <- dt4rtn$tfs[c("HCLS1","STAT4","STAT1","LMO4","ZNF552")]
#' 
#' ##--- run mbrPreprocess
#' rmbr <- mbrPreprocess(gexp=gexp, regulatoryElements1 = tfs1, 
#' regulatoryElements2=tfs2, rowAnnotation=annot)
#' 
#' ##--- run mbrPermutation (set nPermutations>=1000)
#' rmbr <- mbrPermutation(rmbr, nPermutations=30)
#' 
#' ##--- run mbrBootstrap (nBootstrap>=100)
#' rmbr <- mbrBootstrap(rmbr, nBootstrap=30)
#' 
#' ##--- run mbrDpiFilter
#' rmbr <- mbrDpiFilter(rmbr)
#' 
#' ##--- run mbrAssociation (set nPermutations>=1000)
#' rmbr <- mbrAssociation(rmbr, pValueCutoff = 0.05, nPermutations=30)
#'
#' @import RTN 
#' @importFrom stats p.adjust phyper median pnorm sd
#' @importFrom stats cor quantile pt complete.cases
#'
#' @import methods
#' @docType methods
#' @rdname mbrAssociation-methods
#' @aliases mbrAssociation
#' @export

##------------------------------------------------------------------------------
##Inference of duals
setMethod("mbrAssociation",
          "MBR",
          function(object, regulatoryElements1=NULL, regulatoryElements2=NULL, 
                   minRegulonSize=15, pValueCutoff=0.001, pAdjustMethod="BH", 
                   estimator="spearman", nPermutations=1000, miFilter=TRUE, 
                   verbose=TRUE){
            ##--- input check
            if(object@status["Preprocess"]!="[x]")
              stop("NOTE: MBR object is not complete: requires preprocessing!")
            if(object@status["Permutation"]!="[x]")
              stop("NOTE: MBR object is not complete: requires permutation/bootstrap and DPI filter!")  
            if(object@status["DPI.filter"]!="[x]")
              stop("NOTE: MBR object is not complete: requires DPI filter!")
            ##--- gets
            TNI1 <- mbrGet(object, what="TNI1")
            TNI2 <- mbrGet(object, what="TNI2")
            tni_gexp <- tni.get(TNI1, what="gexp")
            tni_para <- tni.get(TNI1, what="para")
            ##--- checks
            mbr.checks(name="minRegulonSize", para=minRegulonSize)
            mbr.checks(name="estimator", para=estimator)
            mbr.checks(name="pValueCutoff", para=pValueCutoff)
            mbr.checks(name="pAdjustMethod", para=pAdjustMethod)
            mbr.checks(name="nPermutations", para=nPermutations)
            mbr.checks(name="miFilter", para=miFilter)
            mbr.checks(name="verbose", para=verbose)
            
            if(verbose) cat("-Checking regulons and regulatory elements...\n")
            
            if(is.null(regulatoryElements1)) {
              regulatoryElements1 <- tni.get(TNI1, "regulatoryElements")
            } else {
              regulatoryElements1 <- .checkRegel(TNI1, regulatoryElements1)
            } 
            if(is.null(regulatoryElements2)) {
              regulatoryElements2 <- tni.get(TNI2, "regulatoryElements")
            } else {
              regulatoryElements2 <- .checkRegel(TNI2, regulatoryElements2)
            }
            mbr.checks(name="numberRegElements", para=regulatoryElements1)
            mbr.checks(name="numberRegElements", para=regulatoryElements2)
            
            ##--- get all regulons (from tnet)
            regulons1 <- tni.get(TNI1, what="regulons")
            regulons2 <- tni.get(TNI2, what="regulons")
            refregulons1 <- tni.get(TNI1, what="refregulons")
            refregulons2 <- tni.get(TNI2, what="refregulons")
            
            ##--- get selected regulons
            regulons1 <- regulons1[regulatoryElements1]
            regulons2 <- regulons2[regulatoryElements2]
            size1 <- sapply(regulons1, length)
            size2 <- sapply(regulons2, length)
        
            ##--- check size
            if( sum(size1>=minRegulonSize)==0 | sum(size2>=minRegulonSize)==0) {
              stop("NOTE: at least one input regulon should be above the 
                   'minRegulonSize' in both TNIs!",call.=FALSE)
            }
            ##--- filter those smaller than minRegulonSize
            idx <- size1 >= (minRegulonSize)
            regulons1 <- regulons1[idx]
            regulatoryElements1 <- regulatoryElements1[idx]
            size1 <- size1[idx]
            ##---
            idx <- size2 >= (minRegulonSize)
            regulons2 <- regulons2[idx]
            regulatoryElements2 <- regulatoryElements2[idx]
            size2 <- size2[idx]
            ##---
            refregulons1 <- refregulons1[regulatoryElements1]
            refregulons2 <- refregulons2[regulatoryElements2]  
            
            ##--- group of regulons and regulatory elements
            regel <- unique(c(regulatoryElements1, regulatoryElements2))
            targets <- unique(c(unlist(regulons1), unlist(regulons2)))
            targets <- unique(c(regel, targets))
            
            ##--- get mi from refnet
            if(verbose) {
                cat("-Extrating inferred regulatory associations...\n")
            }
            mitnet1 <- abs(tni.get(TNI1, what="refnet"))
            mitnet2 <- abs(tni.get(TNI2, what="refnet"))
            mitnet1 <- mitnet1[targets,regulatoryElements1]
            mitnet2 <- mitnet2[targets,regulatoryElements2]
            
            if(verbose) {
              cat("-Mapping network motifs between regulons...\n")
            }
            
            ##--- compute correlation between regulators and all targets
            cortnet1 <- .tni.cor(tni_gexp, mitnet1, estimator=estimator, dg=0, 
                                 asInteger=FALSE, mapAssignedAssociation=FALSE)
            cortnet2 <- .tni.cor(tni_gexp, mitnet2, estimator=estimator, dg=0, 
                                 asInteger=FALSE, mapAssignedAssociation=FALSE)
            
            ##--- map interactions for all potential dual regulons
            interegulons <- .getInterRegulons(refregulons1, refregulons2, regulons1, regulons2)
            
            ##--- compute correlation between regulons
            ## the transformation with double 'cor' calls combines the two distributions
            ## into a single one, making it more symmetrical, while preserveing the 
            ## directionality of the approach.
            regcor <- sapply(regulatoryElements1, function(r1){
              sapply(regulatoryElements2, function(r2){
                tar12 <- interegulons[[r1]][[r2]]
                r1 <- cortnet1[tar12,r1]
                r2 <- cortnet2[tar12,r2]
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
            intregs <- intersect(names(regulatoryElements1),names(regulatoryElements2))
            if(length(intregs)>1){
              diag(regcor[intregs,intregs]) <- NA
              regcor[intregs,intregs][upper.tri(regcor[intregs,intregs])] <- NA
            } else {
              regcor[intregs,intregs] <- NA
            }
            testedDuals <- sum(!is.na(regcor))

            ##--- get statlist
            statlist <- .getStatlist(regcor=regcor, regulatoryElements1, regulatoryElements2)
            statlist <- .getMI(statlist, mitnet1)
            if(miFilter) statlist <- statlist[which(statlist$MI != 0), ]

            ##--- further assessing inferred duals
            if(nrow(statlist) > 0){
              
              ##--- n. tests before any filter
              n.tests <- sum(!is.na(regcor))
              
              ##--- assessing overlap between regulons using 'interegulons'
              if(verbose) {
                cat("-Assessing overlap between regulons...\n")
              }
              overlapStats <- .mbr.overlap(statlist, regulons1, regulons2, refregulons1, 
                                           refregulons2, verbose)
              ##--- adjust Pvalue for n.tests
              overlapStats$Adjusted.Pvalue <- p.adjust(overlapStats$Pvalue, 
                                                       method=pAdjustMethod, n=n.tests )
              overlapStats <- overlapStats[overlapStats$Adjusted.Pvalue<pValueCutoff,]
              statlist <- statlist[rownames(overlapStats),]
              
              ##--- running permutation on R.Regulons
              if(verbose) {
                cat("-Assessing correlation between regulons...\n")
              }
              
              ##--- Permutation on correlation
              corrStats <- .permCorPval(cortnet1, cortnet2, statlist, interegulons, 
                                        estimator, nper=nPermutations, verbose=verbose)
              ##--- adjust Pvalue for n.tests
              corrStats$Adjusted.Pvalue <- p.adjust(corrStats$Pvalue, method=pAdjustMethod, n=n.tests )
              corrStats <- corrStats[corrStats$Adjusted.Pvalue<pValueCutoff,]
              statlist <- statlist[rownames(corrStats),]
              
              #--- align results
              overlapStats <- overlapStats[rownames(statlist),]
              corrStats <- corrStats[rownames(statlist),]
              
              ##--- organize corrStats
              labs <- c("Regulon1","Regulon2", "Universe.Size", "Intersect.Size", "Effect.Size",
                        "Expected.Overlap","Observed.Overlap","Pvalue", "Adjusted.Pvalue")
              overlapStats <- overlapStats[,labs]
              labs <- c("Regulon1","Regulon2", "MI.Regulators", "R.Regulons", 
                        "Pvalue", "Adjusted.Pvalue")
              corrStats <- corrStats[,labs]
              
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
              
            } else {
              warning("No 'dual regulon' has been observed for the input parameters.",
                      call.=FALSE)
              nm <- as.numeric()
              ch <- as.numeric()
              corrStats <- data.frame(Regulon1=ch, Regulon2=ch, MI.Regulators=nm, R.Regulons=nm, 
                                      Pvalue=nm, Adjusted.Pvalue=nm)
              overlapStats <- data.frame(Regulon1=ch, egulon2=ch, Universe.Size=nm, R1.Size=nm, 
                                         R2.Size=nm, Expected.Overlap=nm, Observed.Overlap=nm,
                                         Pvalue=nm, Adjusted.Pvalue=nm)
            }
            
            ##--- para
            mbr_para <- mbrGet(object,what="para")
            sum.info.par <- data.frame(minRegulonSize, estimator, pAdjustMethod, nPermutations, 
                                       stringsAsFactors = FALSE)
            mbr_para$MBR$association['Parameter', ] <- sum.info.par
            
            ##---
            mbr_summary <- mbrGet(object, what="summary")
            mbr_summary$MBR$Duals[,'testedDuals'] <- testedDuals
            ##--- set
            object <- .mbr.set(name="statusUpdate", 
                               para="Association", object=object)
            object <- .mbr.set(name="para", 
                               para=mbr_para, object=object)
            object <- .mbr.set(name="summary", 
                               para=mbr_summary, object=object)
            object <- .mbr.set(name="regulatoryElements1", 
                               para=regulatoryElements1, object=object)
            object <- .mbr.set(name="regulatoryElements2", 
                               para=regulatoryElements2, object=object)
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
#' method \code{\link[RTNduals:mbrAssociation]{mbrAssociation}}.
#' @param priorEvidenceTable An optional 'data.frame' with three columns 
#' representing (1) regulatory elements of 'TNI1', (2) regulatory elements 
#' of 'TNI2', and (3) external evidences between the regulatory elements.
#' @param evidenceColname A single character value specifying a column in 
#' the 'priorEvidenceTable'.
#' @param verbose A single logical value specifying to display detailed 
#' messages (when verbose=TRUE) or not (when verbose=FALSE).
#' @return An \linkS4class{MBR} object with an updated 'data.frame' in the slot 
#' 'results' listing the input additional evidences.
#' @examples
#' ##--- load a dataset for demonstration
#' data("dt4rtn", package = "RTN")
#' gexp <- dt4rtn$gexp
#' annot <- dt4rtn$gexpIDs
#' tfs1 <- dt4rtn$tfs[c("IRF8","IRF1","PRDM1","AFF3","E2F3")]
#' tfs2 <- dt4rtn$tfs[c("HCLS1","STAT4","STAT1","LMO4","ZNF552")]
#' 
#' ##--- run mbrAssociation (set nPermutations>=1000 and nBootstrap>=100)
#' rmbr <- mbrPreprocess(gexp=gexp, regulatoryElements1 = tfs1, 
#' regulatoryElements2=tfs2, rowAnnotation=annot)
#' rmbr <- mbrPermutation(rmbr, nPermutations=10)
#' rmbr <- mbrBootstrap(rmbr, nBootstrap=10)
#' rmbr <- mbrDpiFilter(rmbr)
#' rmbr <- mbrAssociation(rmbr, pValueCutoff = 0.05, nPermutations=10)
#' 
#' ##--- check results
#' results <- mbrGet(rmbr, what="dualsCorrelation")
#' 
#' ##--- add supplementary evidence table
#' ## here we build a 'toy' example using the 'rnorm' function for 
#' ## demonstration purposes only!
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
           function(object, priorEvidenceTable=NULL, evidenceColname, 
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
             if(!is.null(priorEvidenceTable)){
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
             }
             return(object)
           }
)

#' A preprocessing function for objects of class MBR.
#'
#' This function converts one or two TNI class objects and into one MBR class object.
#'
#' @param TNI1 A 'TNI' class object.
#' @param TNI2 Optional 'TNI' class object
#' @param regulatoryElements1 A character vector specifying which 
#' 'TNI1' regulatory elements should be evaluated.
#' @param regulatoryElements2 Optional character vector specifying which 
#' 'TNI2' regulatory elements should be evaluated.
#' @param verbose A single logical value specifying to display detailed messages 
#' (when verbose=TRUE) or not (when verbose=FALSE).
#' @return An \linkS4class{MBR} object.
#' @examples
#' #--- load a dataset for demonstration
#' data("dt4rtn", package = "RTN")
#' gexp <- dt4rtn$gexp
#' annot <- dt4rtn$gexpIDs
#' tfs1 <- dt4rtn$tfs[c("IRF8","IRF1","PRDM1","AFF3","E2F3")]
#' tfs2 <- dt4rtn$tfs[c("HCLS1","STAT4","STAT1","LMO4","ZNF552")]
#' 
#' \dontrun{
#' 
#' ##--- compute a TNI for tfs1
#' tni1 <- tni.constructor(expData=gexp, regulatoryElements=tfs1, rowAnnotation=annot)
#' tni1 <-tni.permutation(tni1)
#' tni1 <-tni.bootstrap(tni1)
#' tni1 <-tni.dpi.filter(tni1)
#' 
#' ##--- compute a TNI for tfs2
#' tni2 <- tni.constructor(expData=gexp, regulatoryElements=tfs2, rowAnnotation=annot)
#' tni2 <-tni.permutation(tni2)
#' tni2 <-tni.bootstrap(tni2)
#' tni2 <-tni.dpi.filter(tni2)
#' 
#' ##--- run tni2mbrPreprocess
#' rmbr <- tni2mbrPreprocess(tni1, tni2)
#' }
#'
#' @import RTN
#' @import methods
#' @docType methods
#' @rdname tni2mbrPreprocess-methods
#' @aliases tni2mbrPreprocess
#' @export

##------------------------------------------------------------------------------
##Combine two TNIs produced separately
setMethod("tni2mbrPreprocess",
          "TNI",
          function (TNI1,  TNI2, 
                    regulatoryElements1=NULL, regulatoryElements2=NULL, 
                    verbose=TRUE){
            if(missing(TNI1)) stop("NOTE: 'tni1' is missing ", call.=FALSE)
            if(missing(TNI2)) TNI2 <- TNI1
            TNI1 <- upgradeTNI(TNI1)
            TNI2 <- upgradeTNI(TNI2)
            mbr.checks(name='tni', para=TNI1)
            mbr.checks(name='tni', para=TNI2)
            if(!is.null(regulatoryElements1)){
              mbr.checks(name="regulatoryElements1", para=regulatoryElements1)
              regulatoryElements1 <- .checkRegel(TNI1, regulatoryElements1)
              TNI1@regulatoryElements <- regulatoryElements1
              if(is.null(regulatoryElements2))
                regulatoryElements2 <- regulatoryElements1
            }
            if(!is.null(regulatoryElements2)){
              mbr.checks(name="regulatoryElements2", para=regulatoryElements2)
              regulatoryElements2 <- .checkRegel(TNI2, regulatoryElements2)
              TNI2@regulatoryElements <- regulatoryElements2
            }
            
            ##---- check processing
            .checkTNIsProcessing(tni1=TNI1, tni2=TNI2, verbose=verbose)
            #---get
            gexp <- tni.get(TNI1, what="gexp")
            regulatoryElements1 <- tni.get(TNI1, what="regulatoryElements")
            regulatoryElements2 <- tni.get(TNI2, what="regulatoryElements")
            ##---- creates MBR object
            object <- .MBRmaker(gexp=gexp,
                                regulatoryElements1=regulatoryElements1,
                                regulatoryElements2=regulatoryElements2)
            #--- Update TNIs
            object <- .mbr.set(name="TNI1", para=TNI1, object=object)
            object <- .mbr.set(name="TNI2", para=TNI2, object=object)
            #---statu update
            TNI_status <- tni.get(TNI1, what="status")
            status <- names(TNI_status[TNI_status=="[x]"])
            object <- .mbr.set(name="statusUpdate", para=status, object=object)
            #---get Updates
            TNI1_summary <- tni.get(TNI1, what="summary")
            TNI2_summary <- tni.get(TNI2, what="summary")
            mbr_summary <- mbrGet(object, what="summary")
            mbr_para <- mbrGet(object, what="para")
            ##--- TNI1 and TNI2
            mbr_summary$TNIs$TNI1 <- TNI1_summary$results$tnet
            colnames(mbr_summary$TNIs$TNI1) <- c('RE', 'Targets', 'Edges')
            mbr_summary$TNIs$TNI2 <- TNI2_summary$results$tnet
            colnames(mbr_summary$TNIs$TNI2) <- c('RE', 'Targets', 'Edges')
            ##--- permutation
            mbr_para$TNIs$perm <- TNI1_summary$para$perm
            ##--- bootstrap
            mbr_para$TNIs$boot <- TNI1_summary$para$boot
            ##---summary dpi.filter
            mbr_para$TNIs$dpi <- TNI1_summary$para$dpi
            #---set
            object <- .mbr.set(name="para", para=mbr_para, object=object)
            object <- .mbr.set(name="summary", para=mbr_summary, object=object)
            return (object)
          }
)

##------------------------------------------------------------------------------
##show summary information on screen
setMethod( "show",
           "MBR",
           function(object){
             cat("an MBR (Motifs Between Regulons) object:\n")
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
#' retrieved from the slots. Options: "TNI1", "TNI2", "regulatoryElements1", 
#' "regulatoryElements2", "dualRegulons", "results", "para", "summary", 
#' "status", "dualsCorrelation", "dualsOverlap", and "dualsCorMatrix"
#' @return A slot content from a object of class 'MBR' \linkS4class{MBR} object
#' @examples
#' ##--- load a dataset for demonstration
#' data("dt4rtn", package = "RTN")
#' gexp <- dt4rtn$gexp
#' annot <- dt4rtn$gexpIDs
#' tfs1 <- dt4rtn$tfs[c("IRF8","IRF1","PRDM1","AFF3","E2F3")]
#' tfs2 <- dt4rtn$tfs[c("HCLS1","STAT4","STAT1","LMO4","ZNF552")]
#' 
#' ##--- run mbrPreprocess
#' rmbr <- mbrPreprocess(gexp=gexp, regulatoryElements1 = tfs1, 
#' regulatoryElements2=tfs2, rowAnnotation=annot)
#' 
#' ##--- get the 'TNI1' slot using 'mbrGet'
#' tni1 <- mbrGet(rmbr, what="TNI1")
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
             optsAssoci <- c("regulatoryElements1", "regulatoryElements2", 
                             "results", "dualRegulons", "dualsCorrelation",
                             "dualsOverlap","dualsCorMatrix")
             ##---get query
             if(what=="TNI1"){
               query <- object@TNI1
             } else if(what=="TNI2"){
               query <- object@TNI2
             } else if(what=="para"){
               query <- object@para
             } else if(what=="summary"){
               query <- object@summary
             } else if(what=="status"){
               query <- object@status
             } else if(what%in%optsAssoci){
               if(object@status["Association"] != "[x]"){
                 warning("Input 'object' needs 'mbrAssociation' evaluation!", call.=FALSE)
                 query <- NULL
               } else {
                 if(what=="regulatoryElements1"){
                   query <- object@regulatoryElements1
                 } else if(what=="regulatoryElements2"){
                   query <- object@regulatoryElements2
                 } else if(what=="dualRegulons"){
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
