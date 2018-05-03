
################################################################################
##########################    plot dual regulons    ############################
################################################################################

#' Plot shared target clouds between dual regulons.
#'
#' This function plots the shared target clouds between a regulon pair.
#'
#' @param object A processed object of class \linkS4class{MBR} evaluated by 
#' the method \code{\link[RTNduals:mbrAssociation]{mbrAssociation}}.
#' @param dualreg.name A string indicating the name of a dual regulon.
#' @param filepath A string indicating the file path where the plot 
#' should be saved.
#' @param ptcols A vector of length 2 indicating colors of negative 
#' and positive correlations, respectively.
#' @return A plot showing targets of dual regulons.
#' @examples
#' ##--- load a dataset for demonstration
#' data("dt4rtn", package = "RTN")
#' gexp <- dt4rtn$gexp
#' annot <- dt4rtn$gexpIDs
#' tfs <- dt4rtn$tfs[c("IRF8","IRF1","PRDM1","AFF3","E2F3")]
#' 
#' ##--- run mbrPreprocess
#' rmbr <- mbrPreprocess(gexp=gexp, regulatoryElements = tfs, 
#' rowAnnotation=annot)
#' 
#' \dontrun{
#' 
#' ##--- run mbrPermutation
#' rmbr <- mbrPermutation(rmbr, nPermutations=1000)
#' 
#' ##--- run mbrBootstrap
#' rmbr <- mbrBootstrap(rmbr, nBootstrap=100)
#' 
#' ##-- run mbrDpiFilter
#' rmbr <- mbrDpiFilter(rmbr)
#' 
#' ##--- run mbrAssociation
#' rmbr <- mbrAssociation(rmbr, nPermutations=1000)
#'
#' ##--- get inferred duals and plot the shared cloud of targets
#' duals <- mbrGet(rmbr, what="dualRegulons")
#' mbrPlotDuals(rmbr, dualreg.name=duals[1])
#'
#'}
#' @importFrom grDevices adjustcolor dev.off pdf colorRampPalette col2rgb
#' @importFrom graphics abline axis par plot.new plot.window points title segments
#' legend
#' @export

##------------------------------------------------------------------------------
mbrPlotDuals <- function(object, dualreg.name, filepath=NULL, 
                         ptcols=c("#006400FF","#CD6600FF")){
  ##----check object class
  mbr.checks(name="object", para=object)
  mbr.checks(name="dualreg.name", para=dualreg.name)
  mbr.checks(name="filepath", para=filepath)
  mbr.checks(name="ptcols", para=ptcols)
  ##---
  mbr_para <- mbrGet(object,"para")
  mbr_estimator <- mbr_para$MBR$association$estimator
  mbr_nper <- mbr_para$MBR$association$nPermutations
  rtni <- mbrGet(object, "TNI")
  rtni_para <- tni.get(rtni, what="para")
  rtni_nper <- rtni_para$perm$nPermutations
  dualstb <- mbrGet(object, what="dualsCorrelation")
  if(!dualreg.name%in%rownames(dualstb)){
    tp <- unlist(strsplit(dualreg.name, split = "~", fixed=TRUE))
    dualreg.name <- paste(tp[2:1], collapse = "~")
    if(!dualreg.name%in%rownames(dualstb)){
      stop("NOTE: all 'dualreg.name' should be listed in 'dualsCorrelation' table!\nsee 'mbrGet' function. \n", call.=FALSE)
    }
    dualstb <- dualstb[dualreg.name,]
    reg1 <- dualstb$Regulon2
    reg2 <- dualstb$Regulon1
  } else {
    dualstb <- dualstb[dualreg.name,]
    reg1 <- dualstb$Regulon1
    reg2 <- dualstb$Regulon2
  }
  r_val <- dualstb$R.Regulons
  r_adjpval <- dualstb$Adjusted.Pvalue
  #---
  if(!is.null(filepath)){
    filename <- paste("dual_", reg1,"_vs_" ,reg2, sep="")
    filename <- paste(path.expand(filepath), filename, sep="/")
  } else {
    filename <- NULL
  }
  .tni.plot.greement(rtni=rtni, duals=c(reg1, reg2), r_val=r_val, r_adjpval=r_adjpval, 
                     filename=filename, ptcols=ptcols, mbr_estimator=mbr_estimator,
                     mbr_nper=mbr_nper, rtni_nper=rtni_nper)
}

##------------------------------------------------------------------------------
##subfunction for 'mbrPlotDuals'
.tni.plot.greement<-function(rtni, duals, r_val, r_adjpval, filename, ptcols, 
                             mbr_estimator, mbr_nper, rtni_nper){
  
  ##--- get labs
  tfs <- tni.get(rtni, "regulatoryElements")
  idx1 <- match(duals, names(tfs))
  idx2 <- match(duals, tfs)
  idxcheck<-which(is.na(idx1))
  idx1[idxcheck]<-idx2[idxcheck]
  duals<-tfs[idx1]
  
  ##--- get correlation mtx
  gexp <- tni.get(rtni, "gexp")
  tnet <- tni.get(rtni, "tnet")
  reftnet <- tni.get(rtni, "refnet")
  tnet<-tnet[,duals]
  reftnet<-reftnet[,duals]
  
  #--- get shared targets (regs intersect)
  #--- ..same idea as described at '.getInterRegulons' function!
  idx1 <- rowSums(reftnet!=0)==2
  idx2 <- rowSums(tnet!=0)>0
  sharedtar <- rownames(tnet)[idx1 & idx2]
  reftnet <- reftnet[c(duals,sharedtar),]
  tnet <- tnet[c(duals,sharedtar),]
  gexp <- gexp[c(duals,sharedtar),]
  
  #--- compute correlation for sharedtar
  xy <- .tni.cor(gexp,tnet,asInteger=FALSE,estimator=mbr_estimator, 
                 mapAssignedAssociation=FALSE)
  xy <- xy[sharedtar,]
  
  ##--- set cols and others
  substr(mbr_estimator, 1, 1) <- toupper(substr(mbr_estimator, 1, 1))
  xlab <- paste(mbr_estimator,"'s correlation between\n",names(duals)[1]," and shared targets", sep="")
  ylab <- paste(mbr_estimator,"'s correlation between\n",names(duals)[2]," and shared targets", sep="")
  xlim=c(-1.0,1.0)
  ylim=c(-1.0,1.0)
  bgcols <- ptcols
  bgcols[1]<-colorRampPalette(c(ptcols[1],"white"))(30)[15]
  bgcols[2]<-colorRampPalette(c(ptcols[2],"white"))(30)[15]
  bgcols<-adjustcolor(bgcols,alpha.f=0.8)
  
  ##---get legend
  if(!is.null(r_val) && !is.null(r_adjpval)){
    p_resolution <- 1/mbr_nper
    if(r_adjpval < p_resolution){
      r_adjpval <- paste("< ",format(p_resolution,digits=2, scientific = T),sep="")
    } else if(r_adjpval < 1e-16){
      r_adjpval <- "< 1e-16"
    } else {
      r_adjpval <- paste("= ",format(r_adjpval,digits=2, scientific = T),sep="")
    }
    tp <- paste("= ",format(r_val,digits=3,nsmall=3),sep="")
    leg1 <- c("Regulons' association",paste("R ", tp, sep=""),
              paste("Adj. Pvalue ", r_adjpval, sep=""))
  } else {
    leg1 <- NULL
  }
  
  ##-------------------------------------------------------------------
  ##---start plot
  if(!is.null(filename)){
    pdf(file=paste(filename,".pdf",sep=""), height=3, width=3)
    cexleg1 <- 0.5
    cexleg2 <- 0.6
    cexleg3 <- 0.8
  } else {
    cexleg1 <- 0.7
    cexleg2 <- 0.7
    cexleg3 <- 0.8
  }
  par(mgp=c(2.2, 0.5, 0),mar=c(3.3, 3.3, 0.2, 0.2) + 0.1)
  plot.new()
  plot.window(ylim=xlim,xlim=ylim)
  par(mgp=c(1.8, 0.5, 0))
  axis(2,cex.axis=cexleg3,las=1,tcl=-0.15,lwd=2)
  title(ylab=ylab,cex.lab=cexleg3)
  par(mgp=c(2, 0.3, 0))
  axis(1,cex.axis=cexleg3,las=1,tcl=-0.15,lwd=2)
  title(xlab=xlab,cex.lab=cexleg3)
  
  ##---add lines
  segments(0,-1,0,1,lwd=1, lty=212, col="grey70")
  segments(-1,0,1,0,lwd=1, lty=212, col="grey70")
  
  ##---add main plot
  if(r_val<0){
    ##---negative Dual
    #compete
    tpp<-xy[(sign(xy[, 1])==1 & sign(xy[, 2])==-1),]
    points(tpp, col=ptcols[1], pch=21, cex=cexleg1, bg=bgcols[1], lwd=0.70)
    #compete
    tpp<-xy[sign(xy[, 1])==-1 & sign(xy[, 2])==1,]
    points(tpp,col=ptcols[1],pch=21,cex=cexleg1,bg="white", lwd=0.70)
    #agree
    tpp<-xy[rowSums(sign(xy))!=0, ]
    points(tpp,col="grey",pch=21,cex=cexleg1,bg="grey", lwd=0.70)
    legend("topright", c(leg1), bty="n", cex = cexleg2, text.font = c(2,1,1,1))
  } else {
    ##---positive Dual
    #co-activation
    tpp<-xy[rowSums(sign(xy))==2, ]
    points(tpp,col=ptcols[2],pch=21,cex=cexleg1,bg="white", lwd=0.70)
    #co-repression
    tpp<-xy[rowSums(sign(xy))==-2, ]
    points(tpp,col=ptcols[2],pch=21,cex=cexleg1,bg=bgcols[2], lwd=0.70)
    #compete
    tpp<-xy[abs(rowSums(sign(xy)))!=2, ]
    points(tpp,col="grey",pch=21,cex=cexleg1,bg="grey", lwd=0.70)
    legend("bottomright", c(leg1), bty="n", cex = cexleg2, text.font = c(2,1,1,1))
  }
  ##--- 
  if(!is.null(filename)){
    dev.off()
    tp <- paste("- file '", filename,".pdf' has been generated!\n", sep="")
    cat(tp)
  }
  
  ##---report
  colnames(xy)<-paste(names(duals),"(R)",sep="")
  nms<-rownames(xy)
  annot<-rtni@rowAnnotation[nms,]
  report<-cbind(annot,format(round(xy,3)))
  invisible(report)
}
