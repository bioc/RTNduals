% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MBR-methods.R
\docType{methods}
\name{mbrAssociation,MBR-method}
\alias{mbrAssociation,MBR-method}
\alias{mbrAssociation}
\title{Motifs analysis and inference of 'dual regulons'.}
\usage{
\S4method{mbrAssociation}{MBR}(
  object,
  regulatoryElements = NULL,
  minRegulonSize = 15,
  doSizeFilter = FALSE,
  pValueCutoff = 0.001,
  pAdjustMethod = "bonferroni",
  estimator = "spearman",
  nPermutations = 1000,
  miFilter = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{object}{A processed object of class \linkS4class{MBR}}

\item{regulatoryElements}{An optional character vector specifying which 
'TNI' regulatory elements should be evaluated. If 'NULL' all regulatory 
elements will be evaluated.}

\item{minRegulonSize}{A single integer or numeric value specifying the 
minimum number of elements in a regulon. Gene sets with fewer than this 
number are removed from the analysis.}

\item{doSizeFilter}{a logical value. If TRUE, negative and positive targets are 
independently verified by the 'minRegulonSize' argument.}

\item{pValueCutoff}{a single numeric value specifying the cutoff for p-values 
considered significant.}

\item{pAdjustMethod}{A single character value specifying the p-value 
adjustment method to be used (see 'p.adjust' function for details).}

\item{estimator}{A character value specifying the estimator used in the 
association analysis. One of "spearman" (default), "kendall", or "pearson".}

\item{nPermutations}{A single integer value specifying the number of 
permutations for deriving p-values associating regulon pairs.}

\item{miFilter}{A single logical value specifying to apply the 'miFilter' 
between two regulators.}

\item{verbose}{A single logical value specifying to display detailed 
messages (when verbose=TRUE) or not (when verbose=FALSE).}
}
\value{
An \linkS4class{MBR} object with two data.frames in the slot 
'results' listing the inferred 'dual regulons' and correspoding statistics.
}
\description{
This function takes an MBR object and compares the shared regulon 
targets in order to test whether regulon pairs agree on the predicted 
downstream effects.
}
\examples{
##--- load a dataset for demonstration
data("tniData", package = "RTN")
gexp <- tniData$expData
annot <- tniData$rowAnnotation
tfs <- c("IRF8","IRF1","PRDM1","E2F3","STAT4","LMO4","ZNF552")

##--- construct a tni object
rtni <- tni.constructor(gexp, regulatoryElements = tfs, rowAnnotation=annot)

##--- compute regulons 
## set nPermutations>=1000
rtni <- tni.permutation(rtni, nPermutations=30)
## set nBootstrap>=100
rtni <- tni.bootstrap(rtni, nBootstrap=30)
## 'eps=NA' estimates threshold from empirical null
rtni <- tni.dpi.filter(rtni, eps=NA)

##--- construct a mbr object
rmbr <- tni2mbrPreprocess(rtni)

##--- run mbrAssociation 
## set nPermutations>=1000
rmbr <- mbrAssociation(rmbr, pValueCutoff = 0.05, nPermutations=30)

}
