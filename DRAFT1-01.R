#'-------------------------------------------------------------------------#
#' Script 1: preparo das matrizes de variancia e covariancia para USP e HEL
#' 
#'-------------------------------------------------------------------------#
# USP data set
#'-------------------------------------------------------------------------#

home.dir <- getwd()
source("env-K.R")

# directorie of data set
setwd("C:/Users/germano/Documents/Doutorado/2019-1/Draft1/USP data") # USP
#setwd("C:/Users/germano/Documents/Doutorado/2019-1/Draft1/HEL data") # HEL

# Phenotypic and molecular data
#'----------------------------------------------------------------------
GY.df <-readRDS("GY570.rds")
Ze <- model.matrix(~ env  - 1, data = GY.df)
Zg <- model.matrix(~ gid  - 1, data = GY.df)
colnames(Ze)<-gsub("env","",colnames(Ze))
colnames(Zg)<-gsub("gid","",colnames(Zg))
LA <- t(chol(readRDS("Ga570H")))
#LD <- t(chol(readRDS("Gd570H")))

# Environmental markers                                                                                         
#----------------------------------------------------------------------
W <- data.frame(as.matrix(readRDS("W.full")))
Wp <- envMarker(W,format = "PCA",ncp=2)$W
colnames(Wp) <- c("PC1","PC2")
Ws <- envMarker(W,format = "012")
Wc <- envMarker(W,format = "rel-01",digits = 5)

# Environmental effects                                             
#----------------------------------------------------------------------
Op <- t(chol(diag(1,nrow=8)));#rownames(Op) = colnames(Op) = colnames(O)
Os <-t(chol(envKernel(df.cov = Ws, Zenv = Ze)$O))
Oc <-t(chol(envKernel(df.cov = Wc, Zenv = Ze)$O))

## Covariate effects                                                 
#----------------------------------------------------------------------

Hs <-t(chol(envKernel(df.cov = Ws, Zenv = Ze)$H))
Hc <-t(chol(envKernel(df.cov = Wc, Zenv = Ze)$H))

Ws <- envK(df.cov = Ws,  df.p =  GY.df,  skip=3)
Wp <- envK(df.cov = Wp,  df.p =  GY.df,  skip=3)
Wc <- envK(df.cov = Wc,  df.p =  GY.df,  skip=3)

# Var-cov for MET BGLR (reproducing GBLUP using spherical blups)                           
#'----------------------------------------------------------------------
# environment
ZEZt  <- var.BGLR(Z = Ze,          cov = "I")
ZE1Zt  <- var.BGLR(Z = Ze, K = Os, cov = "K")
ZE2Zt  <- var.BGLR(Z = Ze, K = Oc, cov = "K")
# genotype
ZGZt  <- var.BGLR(Z = Zg,          cov = "I")
ZAZt  <- var.BGLR(Z = Zg,  K = LA, cov = "K")
# GxE
ZH1Zt <- var.BGLR(Z = Ws,  K = Hs, cov = "K")
ZH2Zt <- var.BGLR(Z = Wc,  K = Hc, cov = "K")
ZH3Zt <- var.BGLR(Z = Wp,          cov = "I")



