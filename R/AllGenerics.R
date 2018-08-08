##generic functions
##------------------------------------------------------------------------------
setGeneric("mbrPreprocess",
           function(gexp, regulatoryElements, verbose=TRUE,...)
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
           function(object, regulatoryElements=NULL, minRegulonSize=15, 
                    doSizeFilter=FALSE, pValueCutoff=0.001, 
                    pAdjustMethod="bonferroni", 
                    estimator='spearman', nPermutations=1000, 
                    miFilter=TRUE, verbose=TRUE)
             standardGeneric("mbrAssociation"), package="RTNduals")
setGeneric("mbrPriorEvidenceTable",
           function(object, priorEvidenceTable, evidenceColname, verbose = TRUE)
             standardGeneric("mbrPriorEvidenceTable"), package="RTNduals")
setGeneric("tni2mbrPreprocess",
           function(tni, regulatoryElements=NULL, verbose = TRUE)
             standardGeneric("tni2mbrPreprocess"), package="RTNduals")
setGeneric("mbrGet",
           function(object, what="status")
             standardGeneric("mbrGet"), package="RTNduals")
