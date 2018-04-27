# Unit tests for MBR-class methods
test_mbr <- function()
{
    data("dt4rtn", package = "RTN")
    tfs <- dt4rtn$tfs[c("FOXM1", "E2F2")]
    ##mbrPreprocess
    rmbr <- mbrPreprocess(gexp=dt4rtn$gexp, regulatoryElements=tfs, 
                          rowAnnotation=dt4rtn$gexpIDs)
    status <- mbrGet(rmbr, what="status")
    checkTrue(status[1]=="[x]")
    ##mbrPermutation
    rmbr <- mbrPermutation(rmbr, nPermutations=10, estimator="pearson")
    status <- mbrGet(rmbr, what="status")
    checkTrue(status[2]=="[x]")
    ##mbrBootstrap
    rmbr <- mbrBootstrap(rmbr, estimator="pearson", nBootstrap=10)
    status <- mbrGet(rmbr, what="status")
    checkTrue(status[3]=="[x]")
    ##mbrDpiFilter
    rmbr <- mbrDpiFilter(rmbr)
    status <- mbrGet(rmbr, what="status")
    checkTrue(status[4]=="[x]")
    ##mbr.combine.TNIs
    tni <- mbrGet(rmbr, what="TNI")
    rmbr <- tni2mbrPreprocess(tni)
    status <- mbrGet(rmbr, what="status")
    checkTrue(all(status[1:4]=="[x]"))
    ##mbrAssociation
    rmbr <- mbrAssociation(rmbr, estimator="pearson", nPermutations=10)
    status <- mbrGet(rmbr, what="status")
    dualsCorrelation <- mbrGet(rmbr, what="dualsCorrelation")
    checkTrue(status[5]=="[x]" && is.data.frame(dualsCorrelation) )
    ##mbr.motifs
    dualsCorrelation <- mbrGet(rmbr, what="dualsCorrelation")
    checkTrue(is.data.frame(dualsCorrelation))
}
