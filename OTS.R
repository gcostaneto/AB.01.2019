#'----------------------------------------------------------------------------------------------
# Optmizing training set                                                                          -----
#'----------------------------------------------------------------------------------------------


OTS <- function(G,Ne,mc.cores=4){
  library(STPGA)
  cat(paste("#----------------------------------------------------------------------------#","\n",sep=""))
  cat(paste("Starting:",Sys.time(),"\n",sep=" "))
  svdG <- svd(G, nu = nrow(G), nv = ncol(G))
  pcsG <- G %*% svdG$v
  rownames(pcsG) <- rownames(G)
  
  TS <- GenAlgForSubsetSelectionNoTest(P = as.matrix(pcsG), ntoselect = Ne,plotiters = F, 
                                       lambda = 1e-5, errorstat = "PEVMEAN", mc.cores = mc.cores)[[1]]
  cat(paste("Optmization Complete:",Sys.time(),"\n",sep=" "))
  return(TS)
}


## proximo passo:

criar dataframe contendo o Ne
rodar OST usando dlply

svdG <- svd(Ga, nu = 100, nv = 100)
pcsG <- Ga %*% svdG$v
rownames(pcsG) <- rownames(G)