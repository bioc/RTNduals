##generic functions
##------------------------------------------------------------------------------
setGeneric("mbrPreprocess",
           function(gexp, regulatoryElements1, regulatoryElements2, 
                    verbose=TRUE,...)
             standardGeneric("mbrPreprocess"), package="RTNduals")
setGeneric("mbrPermutation",
           function(object, verbose = TRUE, ...)
             standardGeneric("mbrPermutation"), package="RTNduals")
setGeneric("mbrBootstrap",
           function(object, verbose = TRUE, ...)
             standardGeneric("mbrBootstrap"), package="RTNduals")
setGeneric("mbrDpiFilter",
           function(object, verbose = TRUE, ...)
             standardGeneric("mbrDpiFilter"), package="RTNduals")
setGeneric("mbrAssociation",
           function(object, regulatoryElements1=NULL, 
                    regulatoryElements2=NULL, minRegulonSize=15, 
                    pValueCutoff=0.001, pAdjustMethod="BH", 
                    estimator='spearman', nPermutations=1000, 
                    miFilter=TRUE, verbose=TRUE)
             standardGeneric("mbrAssociation"), package="RTNduals")
setGeneric("mbrPriorEvidenceTable",
           function(object, priorEvidenceTable = NULL,
                    evidenceColname=NULL, verbose = TRUE)
             standardGeneric("mbrPriorEvidenceTable"), package="RTNduals")
setGeneric("tni2mbrPreprocess",
           function(TNI1, TNI2, regulatoryElements1=NULL, 
                    regulatoryElements2=NULL, verbose = TRUE)
             standardGeneric("tni2mbrPreprocess"), package="RTNduals")
setGeneric("mbrGet",
           function(object, what="status")
             standardGeneric("mbrGet"), package="RTNduals")
