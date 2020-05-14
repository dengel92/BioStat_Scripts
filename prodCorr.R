prodCorr <- function(pr, 
                     do.plot = FALSE,
                     types = c("pearson", "spearman"),
                     dist = c("original",
                              "AOE", "FPC", "hclust", "alphabet")){
  if (xor(is.null(pr), is.null(types)))
    stop("Must have Data and Specify type of correlation")
  if (pr != is.null(pr)){
    tlist <- lapply(pr, typeof)
    lngth <- length(tlist)
    for (i in 1:lngth){
      if (any(tlist[i] == "Character"))
        stop("Data must only have values")
    }
  }
  corm <- rcorr(as.matrix(pr), type = types)
  
  cormrm <- corm$r
  cormp <- corm$P
  
  flattenCorrMatrix <- function(cormat, pmat) {
    ut <- lower.tri(cormat)
    data.frame(
      row = rownames(cormat)[row(cormat)[ut]],
      column = rownames(cormat)[col(cormat)[ut]],
      cor  =(cormat)[ut],
      p = pmat[ut]
    )
  }
  
  corr <- flattenCorrMatrix(cormat = cormrm, pmat = cormp)
  
  if (do.plot==TRUE){
    corrplot(cormrm, type = "lower", order = dist, method = "square")
  }
  return(corr)
}