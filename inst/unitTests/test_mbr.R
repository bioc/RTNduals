# Unit tests for MBR-class methods
test_mbr <- function(){
  data("tniData", package = "RTN")
  tfs <- c("FOXM1", "E2F2")
  rtni <- tni.constructor(tniData$expData, regulatoryElements = tfs, 
                          rowAnnotation=tniData$rowAnnotation)
  rtni <- tni.permutation(rtni, nPermutations=10)
  rtni <- tni.bootstrap(rtni, nBootstrap=10)
  rtni <- tni.dpi.filter(rtni, eps=NA)
  rmbr <- tni2mbrPreprocess(rtni)
  rmbr <- mbrAssociation(rmbr, estimator="pearson", nPermutations=10)
  status <- mbrGet(rmbr, what="status")
  checkTrue( all(status=="[x]") )
  dualsCorrelation <- mbrGet(rmbr, what="dualsCorrelation")
  checkTrue( is.data.frame(dualsCorrelation) )
}
