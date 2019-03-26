##generic functions
##------------------------------------------------------------------------------
setGeneric("tni2mbrPreprocess",
           function(tni, regulatoryElements=NULL)
             standardGeneric("tni2mbrPreprocess"), package="RTNduals")
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
setGeneric("mbrGet",
           function(object, what="status")
             standardGeneric("mbrGet"), package="RTNduals")
