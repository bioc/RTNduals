
##------------------------------------------------------------------------------
##This function is used for argument checking
mbr.checks <- function(name, para, paraSuppl){
  if(name == "gexp") {
    if(!is.matrix(para) || !is.numeric(para)){
      tp <- paste("NOTE: 'gexp' should be a numeric matrix with, genes on rows",
                  "and samples on cols!")
      stop(tp, call.=FALSE)
    }
    else if(is.null(rownames(para)) || is.null(colnames(para)) || 
       length(unique(rownames(para))) < length(rownames(para)) || 
       length(unique(colnames(para))) < length(colnames(para))){
      tp <- paste("NOTE: 'gexp' matrix should be named on rows and cols",
                  "(with unique names)")
      stop(tp, call.=FALSE)
    }
  }
  else if(name == "regulatoryElements"){
    if(!is.null(para)){
      if( !all.characterValues(para) || any(para == "") ){
        tp <- paste("NOTE: 'regulatoryElements' should be a character vector,",
                    "without 'NA' or empty names!")
        stop(tp, call.=FALSE)
      }
      if( any(duplicated(para)) )
        stop("NOTE: 'regulatoryElements' should have unique identifiers!", 
             call.=FALSE)
    }
  }
  else if(name=="object"){
    if(class(para) != 'MBR')
      stop("NOTE: 'object' should be a 'MBR' class object", call.=FALSE)
  }
  else if(name=="doSizeFilter"){
    if(!is.singleLogical(para))
      stop("NOTE: 'doSizeFilter' should be a logical value!", call.=FALSE)
  }
  else if(name=="zlog"){
    if(!is.singleLogical(para))
      stop("NOTE: 'zlog' should be a logical value!", call.=FALSE)
  }
  else if(name=="model"){
    tp <- c("lm", "glm", "coxph")
    if( sum(tp %in% class(para))==0 )
      stop("NOTE: 'model' object should be of class 'lm', 'glm' or 'coxph'!", 
           call.=FALSE)
  }
  else if(name=="showdata"){
    if(!is.singleLogical(para))
      stop("NOTE: 'showdata' should be a logical value!", call.=FALSE)
  }
  else if(name=="verbose"){
    if(!is.singleLogical(para))
      stop("NOTE: 'verbose' should be a logical value!", call.=FALSE)
  }
  else if(name=="miFilter"){
    if(!is.singleLogical(para))
      stop("NOTE: 'miFilter' should be a logical value!", call.=FALSE)
  }
  else if(name == "minRegulonSize"){
    if(!is.singleNumber(para) || !para>0)
      stop("NOTE: 'minRegulonSize' should be a numeric value > 0", call.=FALSE)
  }
  else if(name=="nPermutations") {
    if(!is.singleInteger(para) || length(para)!=1 || para<1)
      stop("'nPermutations' should be an integer >=1 !",call.=FALSE)
  }
  else if(name == "pCutoff"){
    if(!is.singleNumber(para) || (!para>=0) && (!para<=1))
      stop("NOTE: 'pCutoff' should be a numeric value >= 0 and <= 1!", call.=FALSE)
  }
  else if(name == "estimator"){
    if(!is.singleString(para) || !para %in% c("spearman", "kendall", "pearson"))
      stop("NOTE: 'estimator' should be one of 'spearman', 'kendall', 'pearson'!", 
           call.=FALSE)
  }
  else if(name == "dualreg"){
    if(!is.singleString(para))
      stop("NOTE: 'dualreg' should be a single string!", 
           call.=FALSE)
  }
  else if(name == "filepath") {
    if(!is.null(para)){
      if(!is.singleString(para) || !dir.exists(para))
        stop("NOTE: 'filepath' should be a valid single path name!", call.=FALSE)
    }
  }
  else if(name == "fpath") {
    if(!is.null(para)){
      if(!is.singleString(para) || !dir.exists(para))
        stop("NOTE: 'fpath' should be a valid single path name!", call.=FALSE)
    }
  }
  else if (name == "fname"){
    if (!is.singleString(para)) 
      stop("'fname' must be a single character.", call. = FALSE)
    #---check name
    validname <- gsub("[^0-9A-Za-z\\.]", '_',para)
    if(validname!=para){
      stop("NOTE: please provide 'fname' without special charaters or path information!",
           call. = FALSE)
    }
  }
  else if (name == "width"){
    if (!is.singleNumber(para)) 
      stop("'width' must be a single numeric values.", call. = FALSE)
  }
  else if (name == "height"){
    if (!is.singleNumber(para)) 
      stop("'height' must be a single numeric values.", call. = FALSE)
  }
  else if(name == "cols"){
    if(!is.color(para))
      stop("NOTE: 'cols' should be a vector with valid colors!", 
           call.=FALSE)
  }
  else if(name == "datacols"){
    if(!is.color(para))
      stop("NOTE: 'datacols' should be a vector with valid colors!", 
           call.=FALSE)
  }
  else if(name == "zcols") {
    if(!is.color(para) || length(para)!=2)
      stop("NOTE: 'zcols' should be a vector (length = 2) with valid colors!", 
           call.=FALSE)
  }
  else if(name == "ycols") {
    if(!is.color(para) || length(para)!=2)
      stop("NOTE: 'ycols' should be a vector (length = 2) with valid colors!", 
           call.=FALSE)
  }
  else if(name == "regulatoryElements"){
    if(!all.characterValues(para) || any(duplicated(para)) ){
      stop("NOTE: 'regulatoryElements' should be a character vector with unique values !", call. = FALSE)
    }
  } 
  else if(name == "vars"){
    if(!all.characterValues(para) || any(duplicated(para)) || length(para)!=2 ){
      stop("NOTE: 'vars' should be a character vector (length = 2) with unique values!", call. = FALSE)
    }
  }
  else if (name == "xlim"){
    if (!is.numeric(para) || length(para) != 2) 
      stop("'xlim' must be a numeric vector of length 2.", call. = FALSE)
  } 
  else if (name == "ylim"){
    if (!is.numeric(para) || length(para) != 2) 
      stop("'ylim' must be a numeric vector of length 2.", call. = FALSE)
  } 
  else if (name == "zlim"){
    if (!is.numeric(para) || length(para) != 2) 
      stop("'zlim' must be a numeric vector of length 2.", call. = FALSE)
  } 
  else if (name == "zcenter"){
    if (!is.singleNumber(para)) 
      stop("'zcenter' must be a single numeric value.", call. = FALSE)
  } 
  else if (name == "ylab"){
    if (!is.singleString(para)) 
      stop("'ylab' must be a single character.", call. = FALSE)
  } 
  else if (name == "xlab"){
    if (!is.singleString(para)) 
      stop("'xlab' must be a single character.", call. = FALSE)
  }
  else if (name == "zlab"){
    if (!is.singleString(para)) 
      stop("'zlab' must be a single character.", call. = FALSE)
  }
  else if(name=="numberRegElements"){
    if(length(para)<2)
      stop("NOTE: at least 2 regulatory elements should be valid!", 
           call.=FALSE)
  }
  else if(name=="pValueCutoff") {
    if(!is.singleNumber(para) || length(para)!=1 || para>1 || para<0)
      stop("'pValueCutoff' should be an integer or numeric value >=0 and <=1  !",call.=FALSE)
  }
  else if(name=="pAdjustMethod"){
    tp <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
    if(!is.singleString(para) || !(para %in% tp))
      stop("NOTE: 'pAdjustMethod' should be any one of: ",
           paste(tp, collapse = ", "),call.=FALSE)
  }
  else if(name == "priorEvidenceTable"){
    #--- general checks
    if(!is.data.frame(para) || !ncol(para)==3 || !nrow(para)>=1 || 
       is.null(dim(para))){
      tp <- paste("NOTE: 'priorEvidenceTable' should be a 'data.frame' object",
                  "with 3 columns (Regulon1, Regulon2, Evidence)!")
      stop(tp, call.=FALSE)
    }
    if(is.null(colnames(para))){
      stop("NOTE: columns in 'priorEvidenceTable' should be named!", call.=FALSE)
    }
    bl <- paraSuppl %in% colnames(para)
    if(!is.character(paraSuppl) || is.null(paraSuppl) || 
       !is.singleString(paraSuppl) || !bl){
      tp <- paste("NOTE: 'evidenceColname' should be a character value listed in",
                  "the 'priorEvidenceTable' colnames!")
      stop(tp, call.=FALSE)
    }
    if(!is.numeric(para[[paraSuppl]])){
      tp <- paste("NOTE: evidences in the 'priorEvidenceTable' should be",
                  "of numeric or integer type!")
      stop(tp, call.=FALSE)
    }
    #--- set ordering
    col_3 <- which(colnames(para)%in%paraSuppl)
    col_1_2 <- !c(1,2,3) %in% col_3
    col123 <- c(c(1,2,3)[col_1_2],col_3)
    para <- para[,col123]
    colnames(para) <- c("Regulon1","Regulon2",paraSuppl)
    para$Regulon1 <- as.character(para$Regulon1)
    para$Regulon2 <- as.character(para$Regulon2)
    #--- check duplications
    duplnms <- c(paste(para[,1],para[,2],sep="~"), paste(para[,2], para[,1],sep="~"))
    if(sum(duplicated(duplnms))>0) {
      duplnms <- duplnms[duplicated(duplnms)]
      duplnms <- unique(duplnms)
      tp <- paste("NOTE: all pairs of regulators should be unique in 'priorEvidenceTable'.\n", 
                  "The following pairs seem duplicated: ")
      stop(c(paste(tp, paste(duplnms,collapse=", ")), call.=FALSE))
    }
    return(para)
  }
  else if(name=="mbrGet"){
    opts <- c("summary", "status", "results", 
              "dualRegulons", "dualsCorrelation", "dualsOverlap",
              "TNI", "regulatoryElements","para","dualsCorMatrix")
    if(!is.character(para) || length(para)!=1 || !(para %in% opts))
      stop(paste("NOTE: 'what' should be any one of the options:", 
                 paste(opts,collapse = ", ") ), call.=FALSE)
  }
}

##------------------------------------------------------------------------------
.txtcollapse<-function(vec){
  paste("'",paste(vec[-length(vec)], collapse = "', '"),
        "'"," and '",vec[length(vec)],"'!", sep=""
  )
}

##------------------------------------------------------------------------------
is.singleNumber <- function(para){
  (is.integer(para) || is.numeric(para)) && length(para)==1L && !is.na(para)
}
is.singleInteger <- function(para){
  lg <- (is.integer(para) || is.numeric(para)) && length(para)==1L && 
    !is.na(para)
  if(lg) lg <- (para / ceiling(para)) == 1
  return(lg)
}
is.singleString <- function(para){
  is.character(para) && length(para) == 1L && !is.na(para)
}
is.singleLogical <- function(para){
  is.logical(para) && length(para) == 1L && !is.na(para)
}
all.binaryValues <- function(para){
  all( para %in% c(0, 1, NA) )
}
all.integerValues <- function(para){
  lg <- ( all(is.integer(para)) || all(is.numeric(para)) ) && 
    !any(is.na(para))
  if(lg) lg <- all ( (para / ceiling(para)) == 1 )
  return(lg)
}
all.characterValues <- function(para){
  all(is.character(para)) && !any(is.na(para))
}
is.color <- function(x){
  res <- try(col2rgb(x),silent=TRUE)
  return(!"try-error"%in%class(res))
}

