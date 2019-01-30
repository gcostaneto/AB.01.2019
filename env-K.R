#' ================================================================================================== #'
#' Environmental kernels and assist tools for environmental covariate quality and genomic prediciton
#' ================================================================================================== #'
#'
#' Author: Germano Martins F. Costa-Neto  germano.cneto@usp.br
#' Version: 0.7 Piracicaba, Jan 28th
#' ================================================================================================== #'


# provide environmental kerneles for reaction norm models
envKernel<-function(df.cov,Zenv){
  #' df.cov: dataframe of environmental covariates
  #' Zenv  : incidence matrix of environments
  require(psych)
  O <- cor(t(df.cov))
  H <- cor(df.cov)
#  O <- tcrossprod(as.matrix(df.cov))
 # H <- crossprod(as.matrix(df.cov))
 # O <- O/(tr(O)/ncol(O))
#  H <- H/(tr(H)/nrow(H))
  O <- O[match(colnames(Zenv),rownames(O)),match(colnames(Zenv),colnames(O))]
  H <- H[match(colnames(df.cov),rownames(H)),match(colnames(df.cov),colnames(H))]
  O <- O +diag(5E-5,nrow=nrow(O),ncol=ncol(O))
  H <- H +diag(5E-5,nrow=nrow(H),ncol=ncol(H))
  return(list(H=H,O=O))
}


# provide environmental covariates formats for genomic prediction and phenotypic GxE diagnosis
envMarker <- function(df.cov,digits=0, format=c("012","rel-01","PCA"), ncp=2){
  #' df.cov : dataframe of environmental covariates
  #' digits: number of output digits
  #' format: environmental marker (012), 
  #'         centering and scaled to mean 0 and variance 1 (rel-01)
  #'         principal component analysis (PCA)
  #' ncp   : number of principal components used (valid only for format="PCA")
  
  require(scales)
  require(FactoMineR)
  
  COVARIATES <- df.cov
  
  if(format == "rel-01"){
    COVARIATES <- as.matrix(round(scale(COVARIATES,scale = T,center = T),digits))
    return(W=COVARIATES)
  }
  
  if(format == "012"){
    COVARIATES<-round(apply(COVARIATES,2,function(X) scales::rescale(x=X,to=c(0,2))),digits)
    return(W=COVARIATES)}
  
  if(format == "PCA"){
    res.pca <- PCA(COVARIATES,graph=FALSE)
    COVARIATESp <- res.pca$svd$U[,1:ncp]
    rownames(COVARIATESp) <- row.names(COVARIATES)
    return(list(W=COVARIATESp,eig=res.pca$eig))}
}


# provide environmental covariate matrix for genomic prediction
envK <-function(df.cov,df.p,skip=3){
  df.p <-data.frame(df.p)
  df.cov <-data.frame(df.cov)
  df.cov$env <- as.factor(rownames(df.cov))
  W <- as.matrix(merge(df.p,df.cov, by="env")[,-c(1:skip)])
  return(W)
}

# provide variance-covariance matrices for genomic prediciton using BGLR
var.BGLR <- function(Z,K,cov=c("I","K")){
  if(cov == "I"){
    I <- diag(1,nrow=ncol(Z))
    ZKZt <-  tcrossprod(Z %*% I)
    return(ZKZt)
  }
  if(cov == "K"){
    ZKZt <-  tcrossprod(Z %*% K)
    return(ZKZt)
  }
  
}

# provide covariate otimization for genomic prediciton and GxE diagnosis
OtimCov <- function(df.cov,size=1,res.hcpc=T){
  require(FactoMineR)
  oticv<-c()
  Ws.red <- PCA(t(df.cov))
  Ws.red <- HCPC(Ws.red,consol = T,graph = F)
  A<-Ws.red$call$X$clust
  for(i in 1:nlevels(A)){
    line<-sample(1:dim(Ws.red$call$X[Ws.red$call$X$clust == levels(A)[i],])[1], size = size, replace = FALSE)
    oticv<-rbind(oticv,rownames(Ws.red$call$X[Ws.red$call$X$clust == levels(A)[i],])[line])
  }
  COV <- data.frame(cov=oticv)
  if(res.hcpc == TRUE){
    return(list(res.hcpc = Ws.red, newCov=df.cov[,match(COV$cov,colnames(df.cov))]))
  }
  if(res.hcpc == FALSE){
    return(df.cov[,match(COV$cov,colnames(df.cov))])
  }
  
}

