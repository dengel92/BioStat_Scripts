#################################
#                               #
#         Description           #
#                               #
#################################

# dx: is the expression set or provided data being used. Input must be a matrix or dataframe with numerical data.

# sampleAnnot: determines of sample information will be provided to determine phenotypes being tested and associated samples

# sampleAnnot_file: if sampleAnnot is true an annotation source must be provided.

# comp: reference to covariates being used to create model matrix

# gender: determines if you are doing differential analysis for each gender

# gender




#' @import gdata 
#' @import comprehenr 
#' @export

# options(error = recover)

diffMetab <- function (dx,
                       class_id,
                       sampleAnnot_file = NA,
                       gender = F,
                       gender_id,
                       otherCtrl = F,
                       samplID,
                       specCTRL,
                       specTRTMNt) {
  
  if(gender && is.null(gender_id))
    stop("Since comparing within genders need column identifying sample gender")
  if(otherCtrl && is.null(specCTRL) && is.null(specTRTMNt))
    stop("Atypical controls are prefent. Please specify specialized comparisons")
  if(is.null(samplID) && is.null(class_id) && is.null(sampleAnnot_file))
    stop("Need to identify samples, specify treatments, and provide sample annotation")
  if(xor(xor(!is.character(class_id), !is.character(gender_id)), !is.character(samplID)))
    stop("ID inputs must be characters")
    
  sampleAnnot_file[,gender_id] <- factor(sampleAnnot_file[,gender_id])
  sampleAnnot_file[,class_id] <- factor(sampleAnnot_file[,class_id])
  if(gender != F){
    contrast.matrix <- list()
    dx.fit.con <- list()
    samplelst <- list()
    dxlst <- list()
    dx.design <- list()
    dx.fit <- list()
    indices <- list()
    cntrsts <- list()
    res.GenderCompare <- list()
    dir.Gender <- list()
    for (i in levels(sampleAnnot_file[,gender_id])){
      
      samplelst[[i]] <- sampleAnnot_file[startsWith(sampleAnnot_file[,gender_id], i,
                                                    ignore.case = T), ]
      testi <- match(names(dx), samplelst[[i]][, samplID])
      samplelst[[i]]<- samplelst[[i]][na.omit(testi), ]
      
      samplelst[[i]] <- samplelst[[i]][order(samplelst[[i]][,class_id]), ]
      
      indices[[i]] <- match(samplelst[[i]][, samplID], names(dx))
      dxlst[[i]] <- dx[ ,na.omit(indices[[i]])]
      rownames(dxlst[[i]]) <- metAnnot$HMDB_id
      dx.design[[i]] <- model.matrix(~0+samplelst[[i]][,class_id], data = samplelst[[i]])
      colnames(dx.design[[i]]) <- levels(samplelst[[i]][,class_id])
      dx.design[[i]] <- dx.design[[i]][,!apply(dx.design[[i]] == 0, 2, all)]
      dx.fit[[i]] <- lmFit(dxlst[[i]], dx.design[[i]])
      cntrsts[[i]] <- c(to_vec(for(j in levels(samplelst[[i]][,class_id]))
        if(startsWith(j, "co", ignore.case = T))
          for(k in levels(samplelst[[i]][,class_id])[-c(2, 4, 5, 7, 8)])
            paste0(j,"-", k)), to_vec(for (l in 1:length(specCTRL))
              if (any(startsWith(names(as.data.frame(dx.design[[i]])), "gh", ignore.case = T)) && otherCtrl) paste0(specCTRL[l], "-", specTRTMNt[l])))
      contrast.matrix[[i]] <- makeContrasts(contrasts=cntrsts[[i]], levels = dx.design[[i]])
      dx.fit.con[[i]] <- contrasts.fit(dx.fit[[i]], contrast.matrix[[i]])
      dx.fit.con[[i]] <- eBayes(dx.fit.con[[i]])
      res.GenderCompare[[i]] <- topTableF(dx.fit.con[[i]], adjust="BH", sort.by = "none", number=Inf)
      dir.Gender[[i]] <- decideTests(dx.fit.con[[i]], adjust.method = "BH", method="global", lfc = 0.002, p.value = 0.01)
      res <- list(res.GenderCompare, dir.Gender)
    }
    return(res)
    }
  else{
    indices <- match(sampleAnnot_file[, samplID], names(dx))
    data <- dx[ ,na.omit(indices)]
    indices2 <- match(names(dx), sampleAnnot_file[, samplID])
    sampleInfo <- sampleAnnot_file[na.omit(indices2), ]
    dx.design.all <- model.matrix(~0+sampleInfo[,class_id]+sampleInfo[,gender_id], data = sampleInfo)
    colnames(dx.design.all) <- c(levels(sampleInfo[,class_id]), gender_id)
    dx.fit.all <- lmFit(dx[,-1], dx.design.all)
    contrsts <- c(to_vec(for(j in levels(sampleAnnot_file[,class_id]))
      if(startsWith(j, "co", ignore.case = T))
        for(k in levels(sampleAnnot_file[,class_id])[-c(2, 4, 5, 7, 8)])
          paste0(j,"-", k)), to_vec(for (l in 1:length(specCTRL))
            if (any(startsWith(names(as.data.frame(dx.design.all)), "gh", ignore.case = T)) && otherCtrl) paste0(specCTRL[l], "-", specTRTMNt[l])))
    contrast.matrix.all <- makeContrasts(contrasts=contrsts, levels = dx.design.all)
    dx.fit.con.all <- contrasts.fit(dx.fit.all, contrast.matrix.all)
    dx.fit.con.all <- eBayes(dx.fit.con.all)
    res.GenderCompare.all <- topTableF(dx.fit.con.all, adjust="BH", sort.by = "none", number=Inf)
    dir.Gender.all <- decideTests(dx.fit.con.all, adjust.method = "BH", method="global", lfc = 0.002, p.value = 0.01)
    res.all <- list(res.GenderCompare.all, dir.Gender.all)
    return(res.all)
  }
   }