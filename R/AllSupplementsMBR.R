################################################################################
#################  Internal functions for RTNduals-methods  ####################
################################################################################

##------------------------------------------------------------------------------
##check annotation of the input regulatory elements in the 'TNI'
.checkRegel <- function(tni, regel){
  #---check regel
  tp <- sapply(colnames(tni@rowAnnotation), function(i) {
    sum(regel%in%tni@rowAnnotation[, i])
  })
  colid <- names(tp[which.max(tp)])
  idx <- which(tni@rowAnnotation[, colid]%in%regel)
  if(length(idx) < length(regel)) {
    warning("Not all names in 'regulatoryElements' are available in the 'TNI' rowAnnotation!",
            call.=FALSE)
  }
  tp <- tni@rowAnnotation[idx,]
  idx <- tni@regulatoryElements %in% rownames(tp)
  if(sum(idx)==0){
    tp <- paste("NOTE: no names in 'regulatoryElements' has been used to call ",
                "regulons in the provided 'TNI'!", sep="")
    stop(tp, call.=FALSE)
  } else if(sum(idx) < length(regel)){
    tp <- paste("Not all names in 'regulatoryElements' have been used to call ",
                "regulons in the provided 'TNI'!", sep="")
    warning(tp, call.=FALSE)
  }
  regel <- tni@regulatoryElements[idx]
  if(length(regel)<2){
    tp <- paste("NOTE: at least two valid names in 'regulatoryElements'", 
                " are required to call dual regulons!", sep="")
    stop(tp, call.=FALSE)
  }
  return (regel)
}

##------------------------------------------------------------------------------
#return 'regcor' in a table format
.getStatlist <- function(regcor,regulatoryElements){
  coord <- which(!is.na(regcor), arr.ind = TRUE)
  rnames <- rownames(regcor)[coord[, 1]]
  cnames <- colnames(regcor)[coord[, 2]]
  rvals <- regcor[coord]
  statlist <- data.frame(Regulon1=cnames,Regulon2=rnames, 
                         ID_regulon1=regulatoryElements[cnames], ID_regulon2=regulatoryElements[rnames],
                         R.Regulons=rvals, stringsAsFactors=FALSE)
  rownames(statlist) <- paste(statlist$Regulon1, statlist$Regulon2, sep="~")
  return(statlist)
}

##------------------------------------------------------------------------------
##get mutual information inferred between regulators
.getMI <- function(statlist, mitnet){
  tb <- as.matrix(statlist[, c("ID_regulon1", "ID_regulon2")])
  mutinf <- apply(tb, 1, function(x){
    reg1 <- x[1]
    reg2 <- x[2]
    abs(mitnet[reg2, reg1])
  })
  statlist <- cbind(statlist, MI.Regulators=mutinf)
  return (statlist)
}

##------------------------------------------------------------------------------
##compute correlation between regulator and targets
.tni.cor<-function(x, tnet, estimator="pearson",dg=0, asInteger=TRUE, 
                   mapAssignedAssociation=TRUE){
  tfs<-colnames(tnet)
  tar<-rownames(tnet)
  ids<-unique(c(tfs,setdiff(tar,tfs)))
  x=x[ids,]
  x=t(x)
  #--
  pcorm=cor(x[,tfs],x[,tar], method=estimator,use="complete.obs")
  if(asInteger){
    pcorm[pcorm<0]=-1
    pcorm[pcorm>0]=1
  }
  if(length(tfs)>1)diag(pcorm[,tfs])=dg
  #--
  pcorm<-t(pcorm)
  colnames(pcorm)<-tfs
  if(mapAssignedAssociation)pcorm[tnet==0]=NA
  pcorm
}

##------------------------------------------------------------------------------
##compure correlation p-value from permutation analysis
.permCorPval <- function(cortnet, statlist, interegulons,
                         estimator, nper=1000, verbose=TRUE){
  r1names <- statlist$ID_regulon1
  r2names <- statlist$ID_regulon2
  cortnet <- abs(cortnet)
  #--- get min regulon size
  ssize <- sapply(1:nrow(statlist), function(i){
    r1 <- r1names[i]
    r2 <- r2names[i]
    length(interegulons[[r1]][[r2]])
  })
  #--- run permutation
  ntests <- nrow(statlist)
  if(verbose) pb <- txtProgressBar(style=3)
  nulls <- sapply(1:ntests, function(i){
    if(verbose) setTxtProgressBar(pb, i/ntests)
    r1 <- cortnet[,r1names[i]]
    r2 <- cortnet[,r2names[i]]
    n <- ssize[i]
    cor(replicate(nper,sample(r1,n)),sample(r2,n), method=estimator)
  })
  nulls[is.na(nulls)] <- 0
  if(verbose)close(pb)
  # #--- compute p for nominal counts
  # nulls <- t(abs(nulls))
  # obs <- abs(statlist$R.Regulons)
  # ct <- rowSums(obs < nulls)
  # pvals <- (1 + ct)/(1 + nper)
  #--- compute p using pnorm
  nulls <- rbind(nulls,statlist$R.Regulons)
  nlmd <- abs(apply(nulls,2,median))
  nlsd <- abs(apply(nulls,2,sd))
  scores <- abs(statlist$R.Regulons)
  scores <- (scores - nlmd)/nlsd
  pvals <- pnorm(scores, lower.tail=FALSE)
  statlist$Pvalue <- pvals
  return(statlist)
}


##------------------------------------------------------------------------------
##This function takes two related gene sets (i.e. dual regulons) and compute the 
## number of genes expected to occur in both regulons, the actual observed overlap,
## and the pvalues. The test assesses the whether the intersect between refregulons
## is enriched with the high stringent dpi targets.
.mbr.overlap <- function(statlist, regulons, refregulons, verbose=TRUE){
  regpairs <- statlist[, c("ID_regulon1", "ID_regulon2")]
  ##get gene sets: used to test whether the intersect between refregulons
  ##is enriched with the high stringent dpi targets
  interegulons <- .getInterRegulons(refregulons, regulons)
  unirefegulons <- .getUnionRefRegulons(refregulons)
  if(verbose) pb<-txtProgressBar(style=3)
  res <- NULL
  for(i in 1:nrow(regpairs)) {
    if(verbose) setTxtProgressBar(pb, i/nrow(regpairs))
    vecpairs <- as.character(regpairs[i, ])
    ##---
    reg1 <- vecpairs[1]
    reg2 <- vecpairs[2]
    sz_set1 <- length( regulons[[reg1]])
    sz_set2 <- length( regulons[[reg2]])
    ##---
    if(sz_set1 <= sz_set2){
      sz_overlap <- length( intersect(regulons[[reg1]], interegulons[[reg1]][[reg2]]) )
    } else {
      sz_overlap <- length( intersect(regulons[[reg2]], interegulons[[reg1]][[reg2]]) )
    }
    #note: regulons are already conditioned to pairs of regulators,
    #so the universe is the ref. regulons
    sz_universe <- length( unirefegulons[[reg1]][[reg2]] )
    ##---
    tmp <- .regulon.hyper(sz_set1, sz_set2, sz_overlap, sz_universe)
    res <- rbind(res, tmp)
  }
  if(verbose) close(pb)
  res <- cbind(statlist[,1:4], res)
  res
}
.regulon.hyper <- function(sz_set1, sz_set2, sz_overlap, sz_universe) {
  ##number of genes in the universe
  N <- sz_universe
  ##size of gene sets
  m <- sz_set1
  n <- sz_set2
  Nm <- N-m
  ##number of hits gene sets
  k <- sz_overlap
  pvals <- phyper(k-1, m, Nm, n, lower.tail = FALSE)
  ex <- (n/N)*m
  if(m == 0 | n == 0) pvals <- 1
  results <- c(N, m, n, ex, k, pvals)
  names(results) <- c("Universe.Size", "Regulon1.Size", "Regulon2.Size",
                      "Expected.Overlap", "Observed.Overlap", "Pvalue")
  return(results)
}

##------------------------------------------------------------------------------
## Map interactions for all potential dual regulons
## (1) get a broad intersect between regulons (i.e. refregulons)
## (2) then put the focus on targets of higher stringency
.getInterRegulons <- function(refregulons, regulons){
  interegulons <- lapply(names(regulons), function(r1){
    inter <- lapply(names(regulons), function(r2){
      reftar12 <- intersect(refregulons[[r1]],refregulons[[r2]])
      tar12 <- union(regulons[[r1]],regulons[[r2]])
      intersect(tar12,reftar12)
    })
    names(inter) <- names(regulons)
    inter
  })
  names(interegulons) <- names(regulons)
  return(interegulons)
}
.getUnionRefRegulons <- function(refregulons){
  interefegulons <- lapply(names(refregulons), function(r1){
    inter <- lapply(names(refregulons), function(r2){
      union(refregulons[[r1]],refregulons[[r2]])
    })
    names(inter) <- names(refregulons)
    inter
  })
  names(interefegulons) <- names(refregulons)
  return(interefegulons)
}

##------------------------------------------------------------------------------
##check the consistency of priorEvidenceTable
.checkConsistencySuppTable <- function(object, priorEvidenceTable, verbose)
  {
  ##---
  evidenceColname <- colnames(priorEvidenceTable)[3]
  dualsCorrelation <- mbrGet(object, what="dualsCorrelation")
  colnms <- colnames(dualsCorrelation)
  if(evidenceColname%in%colnms)
  {
    cat("-NOTE: evidence table has been already provided, overwriting information...\n")
  }
  ##-----check consistency
  regs <- c(dualsCorrelation$Regulon1,dualsCorrelation$Regulon2)
  regs <- unique(regs)
  regsStab <- c(priorEvidenceTable$Regulon1,priorEvidenceTable$Regulon2)
  regsStab <- unique(regsStab)
  ##---consistency between priorEvidenceTable and annotation
  consc <- (sum(regs%in%regsStab)/length(regs))*100
  if(consc<70){
    tp <- paste("Only ",round(consc,1),"% of the regulatory elements ",
                "are represented in the 'priorEvidenceTable'!\n", sep="")
    warning(tp)
  } else if(consc>=70 && verbose){
    tp <- paste("-",round(consc,1),"% of the regulatory elements ",
                "are represented in the 'priorEvidenceTable'!\n", sep="")
    cat(tp)
  }
}

##------------------------------------------------------------------------------
##update 'priorEvidenceTable'
.updateEvidenceTable <- function (object, priorEvidenceTable, verbose=TRUE){
  evidenceColname <- colnames(priorEvidenceTable)[3]
  dualsCorrelation <- mbrGet(object, what="dualsCorrelation")
  dualsCorrelation[, evidenceColname] <- NA
  ##---
  if(verbose){
    tp <- paste("-Transferring evidences from 'priorEvidenceTable'",
                "to inferred duals...\n", sep=" ")
    cat(tp)
  }
  duals <- rownames(dualsCorrelation)
  ##--- Test A-B ordering, and update if available
  rownames(priorEvidenceTable) <- paste(priorEvidenceTable$Regulon1,priorEvidenceTable$Regulon2, sep = "~")
  suppl <- priorEvidenceTable[rownames(priorEvidenceTable)%in%duals,]
  if(nrow(suppl)>0){
    dualsCorrelation[rownames(suppl),evidenceColname] <- suppl[[evidenceColname]]
  }
  ##--- Test B-A ordering, and update if available
  rownames(priorEvidenceTable) <- paste(priorEvidenceTable$Regulon2,priorEvidenceTable$Regulon1, sep = "~")
  suppl <- priorEvidenceTable[rownames(priorEvidenceTable)%in%duals,]
  if(nrow(suppl)>0){
    dualsCorrelation[rownames(suppl),evidenceColname] <- suppl[[evidenceColname]]
  }
  ##--- update object
  object <- .mbr.set(name="dualsCorrelation", para=dualsCorrelation, object=object)
  return (object)
}

##------------------------------------------------------------------------------
##set the slots of an MBR object
.mbr.set <- function(name, para, object){
  if(name=="para"){
    object@para <- para
  } else if(name=="summary"){
    object@summary <- para
  } else if(name=="status") {
    object@status <- para
  } else if(name=="TNI"){
    object@TNI <- para
  } else if(name=="statusUpdate"){
    object@status[para] <- "[x]"
  } else if(name=="regulatoryElements"){
    object@regulatoryElements <- para
  } else if(name=="dualRegulons"){
    object@dualRegulons <- para
  } else if(name=="dualsCorrelation"){
    object@results$dualsCorrelation <- para
  } else if(name=="dualsOverlap"){
    object@results$dualsOverlap <- para
  }
  return(object)
}


# .merge.tnis <- function (TNI1, TNI2){
#   elreg1 <- tni.get(TNI1, "regulatoryElements")
#   elreg2 <- tni.get(TNI2, "regulatoryElements")
#   elregs <- c(elreg1, elreg2)
#   elregs <- elregs[!duplicated(elregs)]
#   rtni_merge <- new("TNI",gexp = tni.get(TNI1, "gexp"), regulatoryElements = elregs)
#   rtni_merge@rowAnnotation <- object@TNI1@rowAnnotation
#   rtni_merge@para <- tni.get(TNI1, "para")
#   tnet1 <- tni.get(TNI1, "refnet")[, elreg1]
#   tnet2 <- tni.get(TNI2, "refnet")[, setdiff(elreg2,elreg1)]
#   rtni_merge@results$tn.ref <- cbind(tnet1, tnet2)
#   tnet1 <- tni.get(TNI1, "tnet")[, elreg1]
#   tnet2 <- tni.get(TNI2, "tnet")[, setdiff(elreg2,elreg1)]
#   rtni_merge@results$tn.dpi <- cbind(tnet1, tnet2)
#   rtni_merge@status [1:4] <- "[x]"
#   return (rtni_merge)
# }

