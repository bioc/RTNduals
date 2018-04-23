################################################################################
#################  Internal functions for RTNduals-methods  ####################
################################################################################

##------------------------------------------------------------------------------
##check annotation of the input regulatory elements in the 'TNI'
.checkRegel <- function(tni, regulatoryElements){
  #---check regulatoryElements
  tp <- sapply(colnames(tni@rowAnnotation), function(i) {
    sum(regulatoryElements%in%tni@rowAnnotation[, i])
  })
  colid <- names(tp[which.max(tp)])
  idx <- which(tni@rowAnnotation[, colid]%in%regulatoryElements)
  if(length(idx) < length(regulatoryElements)) {
    warning("Not all 'regulatory elements' are available in the 'TNI' rowAnnotation!",
            call.=FALSE)
  }
  regulatoryElements <- tni@rowAnnotation[idx,]
  idx <- match(rownames(regulatoryElements), tni@regulatoryElements)
  if(length(idx)==0){
    tp <- paste("NOTE: no 'regulatory element' has been used to call",
                "regulons in the provided 'TNI'!")
    stop(tp, call.=FALSE)
  } else if(length(idx) < length(regulatoryElements)){
    tp <- paste("Not all input 'regulatory elements' have been used to call",
                "regulons in the provided 'TNI'!")
    warning(tp, call.=FALSE)
  }
  regulatoryElements <- tni@regulatoryElements[idx]
  return (regulatoryElements)
  }

##------------------------------------------------------------------------------
#return 'regcor' in a table format
.getStatlist <- function(regcor,regulatoryElements1, regulatoryElements2){
  coord <- which(!is.na(regcor), arr.ind = TRUE)
  rnames <- rownames(regcor)[coord[, 1]]
  cnames <- colnames(regcor)[coord[, 2]]
  rvals <- regcor[coord]
  statlist <- data.frame(Regulon1=rnames,Regulon2=cnames, 
                         ID_regulon1=regulatoryElements1[rnames], ID_regulon2=regulatoryElements2[cnames],
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
.permCorPval <- function(cortnet1, cortnet2, statlist, interegulons,
                         estimator, nper=1000, verbose=TRUE){
  r1names <- statlist$ID_regulon1
  r2names <- statlist$ID_regulon2
  cortnet1 <- abs(cortnet1)
  cortnet2 <- abs(cortnet2)
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
    r1 <- cortnet1[,r1names[i]]
    r2 <- cortnet2[,r2names[i]]
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
.mbr.overlap <- function(statlist, regulons1, regulons2, refregulons1, 
                         refregulons2, verbose=TRUE){
  regpairs <- statlist[, c("ID_regulon1", "ID_regulon2")]
  ##get gene sets: used to test whether the intersect between refregulons 
  ##is enriched with the high stringent dpi targets
  interegulons <- .getInterRegulons(refregulons1, refregulons2,regulons1, regulons2)
  interefegulons <- .getInterRefRegulons(refregulons1, refregulons2)
  uniregulons <- .getUnionRegulons(regulons1, regulons2)
  unirefegulons <- .getUnionRefRegulons(refregulons1, refregulons2)
  if(verbose) pb<-txtProgressBar(style=3)
  res <- NULL
  for(i in 1:nrow(regpairs)) {
    if(verbose) setTxtProgressBar(pb, i/nrow(regpairs))
    vecpairs <- as.character(regpairs[i, ])
    ##---
    reg1 <- vecpairs[1]
    reg2 <- vecpairs[2]
    ##---
    sz_overlap <- length( interegulons[[reg1]][[reg2]] )
    sz_set1 <- length( interefegulons[[reg1]][[reg2]] )
    sz_set2 <- length( uniregulons[[reg1]][[reg2]] )
    #note: regulons are already conditioned to duals,
    #so the universe is reduced to ref. regulons
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
  names(results) <- c("Universe.Size", "Intersect.Size", "Effect.Size",
                      "Expected.Overlap", "Observed.Overlap", "Pvalue")
  return(results)
}

##------------------------------------------------------------------------------
## Map interactions for all potential dual regulons
## (1) get a broad intersect between regulons (i.e. with refregulons)
## (2) then put the focus on regulons of higher stringency
.getInterRegulons <- function(refregulons1, refregulons2, regulons1, regulons2){
  interegulons <- lapply(names(regulons1), function(r1){
    inter <- lapply(names(regulons2), function(r2){
      reftar12 <- intersect(refregulons1[[r1]],refregulons2[[r2]])
      tar12 <- union(regulons1[[r1]],regulons2[[r2]])
      intersect(tar12,reftar12)
    })
    names(inter) <- names(regulons2)
    inter
  })
  names(interegulons) <- names(regulons1)
  return(interegulons)
}
## .. and here only gets the broad intersect with refregulons
.getInterRefRegulons <- function(refregulons1, refregulons2){
  interefegulons <- lapply(names(refregulons1), function(r1){
    inter <- lapply(names(refregulons2), function(r2){
      intersect(refregulons1[[r1]],refregulons2[[r2]])
    })
    names(inter) <- names(refregulons2)
    inter
  })
  names(interefegulons) <- names(refregulons1)
  return(interefegulons)
}
.getUnionRegulons <- function(regulons1, regulons2){
  interefegulons <- lapply(names(regulons1), function(r1){
    inter <- lapply(names(regulons2), function(r2){
      union(regulons1[[r1]],regulons2[[r2]])
    })
    names(inter) <- names(regulons2)
    inter
  })
  names(interefegulons) <- names(regulons1)
  return(interefegulons)
}
.getUnionRefRegulons <- function(refregulons1, refregulons2){
  interefegulons <- lapply(names(refregulons1), function(r1){
    inter <- lapply(names(refregulons2), function(r2){
      union(refregulons1[[r1]],refregulons2[[r2]])
    })
    names(inter) <- names(refregulons2)
    inter
  })
  names(interefegulons) <- names(refregulons1)
  return(interefegulons)
}

##------------------------------------------------------------------------------
##checks to merge two 'TNIs' objects
.checkTNIsProcessing <- function (tni1, tni2, verbose = TRUE){
  ## checks gexp consistency
  if (verbose)
    cat("-Checking expression matrix consistency...\n")
  if (!identical(tni1@gexp,tni2@gexp))
    stop("NOTE: TNIs should use the same expression matrix.")
  
  ## checks parameter consistency
  if (verbose)
    cat("-Checking parameter consistency...\n")
  tp1 <- unlist(tni1@para)
  tp2 <- unlist(tni2@para)
  idx <- tp1%in%tp2
  if(!all(idx)){
    tp <- paste("TNIs were not computed using the same parameters!",
                "The following patameters seem to differ between the two:\n",
                paste(names(tp1)[idx], collapse = "\n "))
    warning(tp)
  }
  ## checks whether both TNIs have undergone all methods in RTN from
  ## Permutation to DPI filter
  if (verbose)
    cat("-Checking if all TNI methods are completed...\n")
  if(any(tni1@status[1:4] != "[x]") || any(tni2@status[1:4] != "[x]")){
    ## gives feedback on which methods were not run
    if (verbose){
      cat("TNI1: ")
      print(tni1@status)
      cat("TNI2: ")
      print(tni2@status)
    }
    tp <- paste("NOTE: both TNIs must be evaluated by the RTN pipeline,",
                "up to the DPI filter step!")
    stop(tp, call. = FALSE)
  }
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
  } else if(name=="TNI1"){
    object@TNI1 <- para
  } else if(name=="TNI2"){
    object@TNI2 <- para
  } else if(name=="statusUpdate"){
    object@status[para] <- "[x]"
  } else if(name=="regulatoryElements1"){
    object@regulatoryElements1 <- para
  } else if(name=="regulatoryElements2"){
    object@regulatoryElements2 <- para
  } else if(name=="dualRegulons"){
    object@dualRegulons <- para
  } else if(name=="dualsCorrelation"){
    object@results$dualsCorrelation <- para
  } else if(name=="dualsOverlap"){
    object@results$dualsOverlap <- para
  }
  return(object)
}

# ##------------------------------------------------------------------------------
# #compute cor p-values
# .corPval <- function(statlist, regulatoryElements1, regulatoryElements2,
#                      regulons1, regulons2, pAdjustMethod){
#   pvals <- sapply(1:nrow(statlist), function(i) {
#     reg1 <- regulatoryElements1[statlist$Regulon1[i]]
#     reg2 <- regulatoryElements2[statlist$Regulon2[i]]
#     N <- min(length(regulons1[[reg1]]), length(regulons2[[reg2]]))
#     R <- statlist$R.Regulons[i]
#     2 * pt( -abs( R*sqrt((N-2)/(1-R^2)) ), N-2)
#   })
#   adjpvals <- p.adjust(pvals, method = pAdjustMethod)
#   statlist$Pvalue2 <- pvals
#   statlist$Adjusted.Pvalue2 <- adjpvals
#   return(statlist)
#   
# }

