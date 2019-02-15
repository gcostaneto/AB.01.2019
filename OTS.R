#'----------------------------------------------------------------------------------------------
# Optmizing training set                                                                          -----
#'----------------------------------------------------------------------------------------------


OTS <- function(G,Ne,mc.cores=4,make.svd=TRUE){
  library(STPGA)
  cat(paste("#----------------------------------------------------------------------------#","\n",sep=""))
  cat(paste("Starting:",Sys.time(),"\n",sep=" "))
  if(make.svd==TRUE){
  svdG <- svd(G, nu = nrow(G), nv = ncol(G))
  pcsG <- G %*% svdG$v
  rownames(pcsG) <- rownames(G)}
  if(make.svd==FALSE){pcsG <- G }
  
  TS <- GenAlgForSubsetSelectionNoTest(P = as.matrix(pcsG), ntoselect = Ne,plotiters = F, 
                                       lambda = 1e-5, errorstat = "PEVMEAN", mc.cores = mc.cores)[[1]]
  cat(paste("Optmization Complete:",Sys.time(),"\n",sep=" "))
  return(TS)
}
