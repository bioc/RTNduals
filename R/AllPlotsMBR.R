
################################################################################
##########################     plot interactions     ###########################
################################################################################

#' Plot interaction effects between two continuous variables.
#'
#' This function plots the interaction effects between two continuous variables 
#' for linear, logistic, or Cox models.
#'
#' @param model An object of class 'lm', 'glm', or 'coxph'.
#' @param vars A character vector of length 2 with the names of two varibles in the 'model'.
#' @param xlim A numeric vector of length 2, i.e. xlim = c(x1, x2), 
#' indicating the x limits of the plot. If xlim = NULL, it will be derevided  
#' from the observed data ranges.
#' @param ylim A numeric vector of length 2, i.e. ylim = c(y1, y2), 
#' indicating the y limits of the plot. If ylim = NULL, it will be derevided 
#' from the observed data ranges.
#' @param zlim A numeric vector of length 2, i.e. zlim = c(z1, z2), 
#' indicating the z limits of the plot. If zlim = NULL, it will be derevided 
#' from the observed data ranges.
#' @param xlab A string with the label for the x-axis.
#' @param ylab A string with the label for the y-axis.
#' @param zlab A string with the label for the z-axis.
#' @param zcenter A numeric value indicating a z value to center the color scale.
#' @param zlog A logical value indicating whether the z axis is to be logarithmic.
#' @param zcols A vector of length 2 indicating a diverging color scheme 
#' for the z-axis varible.
#' @param ycols A vector of length 2 indicating a diverging color scheme 
#' for the y-axis varible (only used when type='2D').
#' @param showdata A logical value indicating whether to show the original data 
#' used to fit linear model.
#' @param datacols  When 'showdata = TRUE', this can be a named vector of colors for 
#' data points (names should match samples in the 'model' object). Alternativelly, 
#' it can be a single color value.
#' @param fname A string. The name of the PDF file which will contain the plot.
#' @param fpath A string. The directory where the file will be saved.
#' @param width A numeric value. The width of the plot.
#' @param height A numeric value. The height of the plot.
#' @param plotype A string indicating '2D' of '3D' plot type. If plotype = '2D', 
#' the z-axis (and all related parameters) is transposed to the y-axis.
#' @param plotpdf A logical value.
#' @return A interaction plot.
#' @examples
#' #-- Example of simulated data, with response variable modelled by:
#' #-- (1) Main effects of 'reg1' and 'reg2'
#' #-- (2) Interaction effects between 'reg1' and 'reg2' 
#' #-- (3) Additional random uniform noise
#' reg1 <- rnorm(1000,0,2)
#' reg2 <- rnorm(1000,0,2)
#' response <- 3*reg1 + 2*reg1*reg2 + runif(1000,0,2)
#' dataset <- data.frame(reg1, reg2, response)
#' model <- lm(response ~ reg1*reg2, data=dataset)
#' mbrPlotInteraction(model, vars=c("reg1","reg2"))
#'
#' @importFrom stats model.frame predict coef terms
#' @importFrom graphics image layout mtext plot
#' @importFrom grDevices adjustcolor rgb2hsv hsv
#' @export

##------------------------------------------------------------------------------
mbrPlotInteraction <- function(model, vars, xlim = NULL, ylim = NULL, zlim = NULL, 
                               xlab = NULL, ylab = NULL, zlab = NULL, zcenter = NULL, 
                               zlog = NULL, zcols = c("#008080ff","#d45500ff"), 
                               ycols = c("#4A97C9","#D92522"),
                               showdata = FALSE, datacols = "grey50", 
                               fname = "interplot", fpath = ".", 
                               width = 4.5, height = 4, plotype = c("3D","2D"),
                               plotpdf = FALSE){
  #--- checks
  mbr.checks(name="model", para=model)
  mbr.checks(name="vars", para=vars)
  if(!is.null(xlim)) mbr.checks(name="xlim", para=xlim)
  if(!is.null(ylim)) mbr.checks(name="ylim", para=ylim)
  if(!is.null(zlim)) mbr.checks(name="zlim", para=zlim)
  if(!is.null(xlab)) mbr.checks(name="xlab", para=xlab)
  if(!is.null(ylab)) mbr.checks(name="ylab", para=ylab)
  if(!is.null(zlab)) mbr.checks(name="zlab", para=zlab)
  if(!is.null(zlog)) mbr.checks(name="zlog", para=zlog)
  if(!is.null(zcenter)) mbr.checks(name="zcenter", para=zcenter)
  mbr.checks(name="zcols", para=zcols)
  mbr.checks(name="ycols", para=ycols)
  mbr.checks(name="showdata", para=showdata)
  mbr.checks(name="datacols", para=datacols)
  mbr.checks(name="fname", para=fname)
  mbr.checks(name="fpath", para=fpath)
  mbr.checks(name="width", para=width)
  mbr.checks(name="height", para=height)
  mbr.checks(name="plotpdf", para=plotpdf)
  plotype <- match.arg(plotype)
  .mbrPlotInteraction(model, vars, xlim, ylim, zlim, xlab, ylab, zlab, 
                      zcenter, zlog, zcols, ycols, showdata, datacols, 
                      fname, fpath, width, height, plotype, plotpdf)
  
}

##------------------------------------------------------------------------------
.mbrPlotInteraction <- function(model, vars, xlim, ylim, zlim, xlab, ylab, 
                                zlab, zcenter, zlog, zcols, ycols, showdata, 
                                datacols, fname, fpath, width, height, 
                                plotype, plotpdf){
  vartypes <- attr(terms(model),"dataClasses")
  if(!all(vars%in%names(vartypes)))
    stop("All 'vars' should be listed in the 'model'!", call. = FALSE) 
  if(!all(vartypes[vars]=="numeric"))
    stop("All 'vars' in the model should be 'numeric'!", call. = FALSE)
  
  #--- set xy labels
  labs <- names(vars)
  if(is.null(labs) || !all.characterValues(labs))
    labs <- vars
  if(is.null(xlab))
    xlab = labs[1]
  if(is.null(ylab))
    ylab = labs[2]
  
  #-- get observed data and set xylim
  obdata <- stats::model.frame(model, drop.unused.levels=TRUE)
  obdata <- obdata[,vars]
  if(showdata){
    if(is.null(xlim)) 
      xlim <- range(obdata[,labs[1]])*1.04
    if(is.null(ylim))
      ylim <- range(obdata[,labs[2]])*1.04
  } else {
    if(is.null(xlim)) 
      xlim <- range(obdata)*1.04
    if(is.null(ylim))
      ylim <- xlim
  }
  
  #--- check datacols and obdata
  if(length(datacols)>1){
    if(length(datacols)!=nrow(obdata))
      stop("'datacols' should provide one color for each sample in the 'model' object!", 
           call. = FALSE)
    nms1 <- rownames(obdata)
    nms2 <- names(datacols)
    if(!is.null(nms1)){
      if(is.null(nms2))
        if(is.null(nms2) || !all(nms2%in%nms1) ){
          stop("names in 'datacols' should match sample names in the 'model' object!", 
               call. = FALSE)
        }
    }
    datacols <- datacols[nms1]
  }
  
  #--- set xy lims
  labxlim <- .prettylabs(xlim)
  labylim <- .prettylabs(ylim)
  
  #--- get vars
  dtvars <- .get.vars(vars,labxlim$lim,labylim$lim)
  
  #-- get predictions
  tp <- class(model)
  if("lm"%in%tp || "glm"%in%tp){
    prmat <- .get.predictions.glm(model,dtvars)
    if(is.null(zlab)) zlab = "Response"
    if(is.null(zlog)) zlog = FALSE
  } else if("coxph"%in%class(model)){
    prmat <- .get.predictions.coxph(model,dtvars)
    if(is.null(zlab)) zlab = "HR"
    if(is.null(zlog)) zlog = TRUE
  } else {
    stop("'model' object should be of class 'lm', 'glm', or 'coxph'!", 
         call. = FALSE) 
  }
  if(!is.null(model$pAdjustInteraction)){
    lab_pvalue <- model$pAdjustInteraction
  } else {
    summ <- coef(summary(model))
    tp1 <- paste(vars,collapse = ":")
    tp2 <- paste(rev(vars),collapse = ":")
    if(tp1%in%rownames(summ)){
      lab_pvalue <- summ[tp1,4]
    } else if(tp2%in%rownames(summ)){
      lab_pvalue <- summ[tp2,4]
    } else {
      lab_pvalue <- NA
    }
  }
  
  #-- set zlim
  if(is.null(zlim))
    zlim <- range(prmat)
  if(plotype=="3D"){
    nn <- 100
  } else {
    nn <- 10
  }
  if(zlog){
    if(is.null(zcenter))zcenter <- 1
    labzlim <- .prettycenter(zlim, n=nn, center=zcenter)
  } else {
    if(is.null(zcenter)){
      labzlim <- .prettynum(zlim, n = nn)
    } else {
      labzlim <- .prettycenter(zlim, n=nn, center=zcenter)
    }
  }
  
  #-- call .interaction.plot
  if(plotpdf){
    fname <- gsub(".pdf", '',fname, ignore.case = TRUE)
    fname <- paste(fname,".pdf", sep = "")
    pdf(file = paste(fpath, "/", fname, sep = ""), width = width, height = height)
  }
  if(plotype=="3D"){
    .plot3D(prmat=prmat, obdata=obdata, labxlim=labxlim, labylim=labylim, 
            labzlim=labzlim, xlab=xlab, ylab=ylab, zlab=zlab, zcols=zcols, 
            zlog=zlog, showdata=showdata, 
            datacols=datacols, lab_pvalue=lab_pvalue)
  } else {
    .plot2D(prmat=prmat, dtvars=dtvars, labxlim=labxlim, labylim=labylim, 
            labzlim=labzlim, xlab=xlab, ylab=ylab, zlab=zlab, 
            zlog=zlog, ycols=ycols)
  }
  if(plotpdf){
    dev.off()
    tp1 <- c("NOTE: file '",fname,"' should be available either in the\n")
    tp2 <- c("working directory or in a user's custom directory!\n")
    message(tp1,tp2)
  }
}

##------------------------------------------------------------------------------
.plot2D <- function(prmat, dtvars, labxlim, labylim, labzlim, xlab, 
                    ylab, zlab, zlog, ycols){
  
  #-- set prmat and zlim
  prmat[prmat<labzlim$lim[1]] <- NA
  prmat[prmat>labzlim$lim[2]] <- NA
  
  #-- set xy lims (note that z is at y-axis here)
  if(zlog){
    log="y"
    zlb <- labzlim$labs
    zat <- labzlim$at
  } else {
    log=""
    zlb <- zat <- labzlim$xlabs
  }
  ract1 <- cbind(dtvars$var1,prmat[,1])
  ract2 <- cbind(dtvars$var1,prmat[,ncol(prmat)])
  #-- plot
  op <- par(no.readonly = TRUE)
  par(mar=c(4,4.3,2.8,4.5), mgp=c(2, 0.4, 0.2), cex=1)
  plot(NA, ylim=labzlim$lim, xlim=labxlim$lim, axes = F, xlab="", ylab="", log = log)
  points(ract1[,1],ract1[,2], type = "l", lwd=2.5, col=ycols[1])
  points(ract2[,1],ract2[,2], type = "l", lwd=2.5, col=ycols[2])
  axis(1, las=1, tcl=-0.2, lwd=1.8, cex.axis=1.2, labels = labxlim$labs, at = labxlim$labs)
  axis(2, las=1, tcl=-0.2, lwd=1.8, cex.axis=1.2, labels = zlb, at = zat)
  title(xlab = xlab, mgp=c(2.2, 0.4, 0), cex.lab=1.2)
  title(ylab = zlab, mgp=c(3, 0.4, 0), cex.lab=1.2)
  par(xpd=TRUE)
  ylim <- format(format(labylim$lim), width = nchar(ylab))
  legend(x = "topright", legend = rev(ylim), col = rev(ycols), bty = "n", pt.cex = 1.2, 
         cex=0.8, title = ylab, title.adj = 0, pch = 15, inset = c(-0.25,0))
  par(op)
}

##------------------------------------------------------------------------------
.plot3D <- function(prmat, obdata, labxlim, labylim, labzlim, xlab, 
                    ylab, zlab, zcols, zlog=TRUE, showdata=FALSE, 
                    datacols="grey80",lab_pvalue){
  
  #-- set prmat and zlim
  prmat[prmat<labzlim$lim[1]] <- labzlim$lim[1]
  prmat[prmat>labzlim$lim[2]] <- labzlim$lim[2]
  
  #-- get labs and color palette
  bks <- labzlim$at
  ptlow <- colorRampPalette(c(zcols[1], "white"))(sum(bks<=labzlim$center))
  pthigh <- colorRampPalette(c("white", zcols[2]))(sum(bks>labzlim$center))[-1]
  pal <- c(ptlow,pthigh)
  
  #-- plot
  op <- par(no.readonly=TRUE)
  layout(matrix(c(1,2), ncol=2), widths = c(4,1))
  par(mar=c(4,4.3,2.8,0.5), mgp=c(2, 0.4, 0.2), cex=1)
  image(prmat, col = pal, breaks = bks, axes=F)
  axis(1,las=1, tcl=-0.2, lwd=1.8, at = labxlim$at, labels = labxlim$labs, cex.axis=1.2)
  axis(2,las=1, tcl=-0.2, lwd=1.8, at = labylim$at, labels = labylim$labs, cex.axis=1.2)
  
  #-- add obdata
  if(showdata){
    tp <- obdata
    tp[,1] <- tp[,1] + abs(labxlim$lim[1])
    tp[,1] <- tp[,1]/(labxlim$lim[2]*2)
    tp[,2] <- tp[,2] + abs(labylim$lim[1])
    tp[,2] <- tp[,2]/(labylim$lim[2]*2)
    idx <- rowSums(tp>1 | tp<0)>0
    lncols <- lighter(datacols, factor = -0.4)
    lncols <- adjustcolor(lncols, alpha.f = 0.6)
    bgcols <- lighter(datacols, factor = 1)
    bgcols <- adjustcolor(bgcols, alpha.f = 0.9)
    points(x=tp[,1],y=tp[,2], cex=0.7, pch=21, col=lncols, bg=bgcols, lwd=0.6)
  }
  
  #-- add labs
  title(ylab = ylab, adj = 0.5, cex.lab = 1.2, mgp = c(if(labylim$nchar<=2) 1.8 else 2.5, 0.4, 0))
  title(xlab = xlab, adj = 0.5, cex.lab = 1.2, mgp = c(1.6, 0.4, 0))

  #--- add leg
  par(xpd=TRUE)
  if(!is.na(lab_pvalue)){
    lab_pvalue <- paste("Interaction P: ", format(lab_pvalue, digits = 3, scientific = TRUE))
    legend("topright", cex = 0.8, legend = lab_pvalue, bty = "n", inset = c(0,-0.1))
  }
  par(mar=c(11,0,4,3.5), mgp=c(3, 0.4, 0), cex=1)
  image(matrix(bks, ncol=length(bks)), col = pal, breaks = bks, axes=F)
  at <- 0:(length(labzlim$at)-1)
  at <- at/(length(labzlim$at)-1)
  idx <- !is.na(labzlim$labs)
  axis(4,labels = labzlim$labs[idx], at = at[idx],lwd=2,las=1, tcl=-0.2, cex.axis = 0.8)
  mtext(zlab, side = 3, adj = 0, line = 0.3)
  par(op)
}

##------------------------------------------------------------------------------
#--- var1*var2 for the interaction matrix
.get.vars <- function(vars, xlim, ylim, n=64){
  var1 <- seq(xlim[1],xlim[2], length=n)
  var2 <- seq(ylim[1],ylim[2], length=n)
  dt <- matrix(data=NA, nrow=n*n, ncol=2)
  for (i in 1:n){
    for (j in 1:n){
      dt[j+(i-1)*n,] <- c(var1[j], var2[i])
    }
  }
  colnames(dt) <- vars
  return(list(dt=dt,var1=var1, var2=var2, n=n))
}

##------------------------------------------------------------------------------
#predictions from lm or glm
.get.predictions.glm <- function(model,dtvars,n){
  pm <- predict(model, newdata=data.frame(dtvars$dt), type="response")
  prmat <- array(pm,dim=c(dtvars$n,dtvars$n))
  return(prmat)
}

##------------------------------------------------------------------------------
#predictions from coxph
.get.predictions.coxph<-function(model,dtvars){
  pm <- predict(model, newdata=data.frame(dtvars$dt), type="risk")
  prmat <- array(pm,dim=c(dtvars$n,dtvars$n))
  return(prmat)
}

##------------------------------------------------------------------------------
.prettycenter <- function(lim, n=10, center=1){
  if(center<lim[1] || center>lim[2]) center <- mean(lim)
  if(center>=1)lim[2] <- ceiling(lim[2])
  lim <- signif(lim, digits = 2)
  xlabs <- sort(c(center, lim))
  xat <- c(seq(lim[1], center, length.out = n), 
           seq(center, lim[2], length.out = n)[-1])
  labs <- c(format(xat, scientific = FALSE, digits=1))
  labs[!xat %in% xlabs] <- NA
  lim <- range(xlabs)
  list(xlabs=xlabs, labs = labs, at = xat, lim=lim, center=1)
}
.prettynum <- function(lim, n=10){
  xlabs <- pretty(lim)
  center <- xlabs[ceiling(length(xlabs)/2)]
  lim <- range(xlabs)
  labs2 <- seq(lim[1],lim[2], length.out = n)
  xat <- labs2 <- sort(unique(c(xlabs,labs2)))
  labs <- format(labs2, scientific = FALSE, digits=1)
  idx <- labs2 %in% xlabs
  labs[!idx] <- NA
  list(xlabs=xlabs, labs = labs, at = xat, lim=lim, center=center)
}
.prettylabs <- function(lim){
  if(length(lim)==1)
    lim <- c(-lim,lim)
  labs <- pretty(lim)
  lim <- range(labs)
  at <- labs+abs(lim[1])
  at <- at/max(at)
  if(length(labs)>5 & max(labs)>1){
    idx <- labs == round(labs)
    if(sum(idx)>sum(!idx)){
      labs <- labs[idx]
      at <- at[idx]
    } else {
      labs <- labs[!idx]
      at <- at[!idx]
    }
  }
  nch <- max(nchar(labs))
  return(list(labs=labs, at=at, lim=lim, nchar=nch))
}

##------------------------------------------------------------------------------
#set lighter colors
lighter <- function(color, factor = 0.2) {
  x <- rgb2hsv(col2rgb(color))
  v <- pmax(pmin(x[3, ] + factor, 1), 0)
  hsv(h = x[1, ], s = x[2, ], v = v)
}

################################################################################
##########################    plot dual regulons    ############################
################################################################################

#' Plot shared targets between regulons.
#'
#' This function plots the shared targets for a regulon pair.
#'
#' @param object A processed object of class \linkS4class{MBR} evaluated by 
#' the method \code{\link{mbrAssociation}}.
#' @param dualreg A string indicating the name of a dual regulon.
#' @param filepath A string indicating the file path where the plot 
#' should be saved.
#' @param cols A vector of length 2 indicating a diverging color scheme 
#' for negative and positive correlations, respectively.
#' @return A plot showing targets of dual regulons.
#' @examples
#' 
#' ##--- load a dataset for demonstration
#' data("dt4rtn", package = "RTN")
#' gexp <- dt4rtn$gexp
#' annot <- dt4rtn$gexpIDs
#' tfs <- dt4rtn$tfs[c("IRF8","IRF1","PRDM1","AFF3","E2F3")]
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
#' \dontrun{
#'
#' ##--- get inferred duals and plot the shared cloud of targets
#' duals <- mbrGet(rmbr, what="dualRegulons")
#' mbrPlotDuals(rmbr, dualreg=duals[1])
#'
#'}
#' @importFrom grDevices adjustcolor dev.off pdf colorRampPalette col2rgb
#' @importFrom graphics abline axis par plot.new plot.window points title segments
#' legend
#' @export

##------------------------------------------------------------------------------
mbrPlotDuals <- function(object, dualreg, filepath=NULL, 
                         cols=c("#006400FF","#CD6600FF")){
  ##----checks
  mbr.checks(name="object", para=object)
  mbr.checks(name="dualreg", para=dualreg)
  mbr.checks(name="filepath", para=filepath)
  mbr.checks(name="cols2", para=cols)
  ##---
  mbr_para <- mbrGet(object,"para")
  mbr_estimator <- mbr_para$MBR$association$estimator
  mbr_nper <- mbr_para$MBR$association$nPermutations
  rtni <- mbrGet(object, "TNI")
  rtni_para <- tni.get(rtni, what="para")
  rtni_nper <- rtni_para$perm$nPermutations
  dualstb <- mbrGet(object, what="dualsCorrelation")
  if(!dualreg%in%rownames(dualstb)){
    tp <- unlist(strsplit(dualreg, split = "~", fixed=TRUE))
    dualreg <- paste(tp[2:1], collapse = "~")
    if(!dualreg%in%rownames(dualstb)){
      stop("'dualreg' should be listed in 'dualsCorrelation' table!\nsee 'mbrGet' function. \n", call.=FALSE)
    }
    dualstb <- dualstb[dualreg,]
    reg1 <- dualstb$Regulon2
    reg2 <- dualstb$Regulon1
  } else {
    dualstb <- dualstb[dualreg,]
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
                     filename=filename, cols=cols, mbr_estimator=mbr_estimator,
                     mbr_nper=mbr_nper, rtni_nper=rtni_nper)
}

##------------------------------------------------------------------------------
##subfunction for 'mbrPlotDuals'
.tni.plot.greement<-function(rtni, duals, r_val, r_adjpval, filename, cols, 
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
  bgcols <- cols
  bgcols[1]<-colorRampPalette(c(cols[1],"white"))(30)[15]
  bgcols[2]<-colorRampPalette(c(cols[2],"white"))(30)[15]
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
    points(tpp, col=cols[1], pch=21, cex=cexleg1, bg=bgcols[1], lwd=0.70)
    #compete
    tpp<-xy[sign(xy[, 1])==-1 & sign(xy[, 2])==1,]
    points(tpp,col=cols[1],pch=21,cex=cexleg1,bg="white", lwd=0.70)
    #agree
    tpp<-xy[rowSums(sign(xy))!=0, ]
    points(tpp,col="grey",pch=21,cex=cexleg1,bg="grey", lwd=0.70)
    legend("topright", c(leg1), bty="n", cex = cexleg2, text.font = c(2,1,1,1))
  } else {
    ##---positive Dual
    #co-activation
    tpp<-xy[rowSums(sign(xy))==2, ]
    points(tpp,col=cols[2],pch=21,cex=cexleg1,bg="white", lwd=0.70)
    #co-repression
    tpp<-xy[rowSums(sign(xy))==-2, ]
    points(tpp,col=cols[2],pch=21,cex=cexleg1,bg=bgcols[2], lwd=0.70)
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
