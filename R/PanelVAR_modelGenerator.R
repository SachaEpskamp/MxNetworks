# Generate stationary PanelVAR model (internally):
panelVAR_modelGen_stat <- function(
  # data, # Data
  covMat,
  means,
  sampleSize,
  designMatrix,
  nNode,
  nTime,
  kappa_mu,
  sigma_mu,
  kappa_zeta,
  sigma_zeta,
  beta,
  mu,
  temporal = TRUE,
  contemporaneous = TRUE,
  betweenSubjects = TRUE,
  startValues = list(),
  name = "PanelVAR",
  groupEqual = "none",
  group = 1
){
  # Some defaults:
  # if (missing(beta)){
  #   beta <- matrix(NA,nNode,nNode)
  # }
  # if (missing(kappa_zeta) && missing(sigma_zeta)){
  #   kappa_zeta <- matrix(NA,nNode,nNode)
  # }
  # if (missing(kappa_mu) && missing(sigma_mu)){
  #   kappa_mu <- matrix(NA,nNode,nNode)
  # }
  # if (missing(mu)){
  #   mu <- rep(NA,nNode)
  # }
  # if (!is.matrix(mu)){
  #   mu <- matrix(mu)
  # }
  
  ### check Design matrix ###
  if (missing(designMatrix)){
    stop("'designMatrix' argument may not be missing. See ?panelVAR for more details.")
  }
  
  # check for duplicates:
  if (any(duplicated(na.omit(c(designMatrix))))){
    stop("Multiple variables with the same name found in the designMatrix")
  }

  # Estraxt info:
  nNode <- nrow(designMatrix)
  nTime <- ncol(designMatrix)
  
  # Check row and column names:
  if (is.null(rownames(designMatrix))){
    rownames(designMatrix) <- paste0("V",seq_len(nNode))
  }
  if (is.null(colnames(designMatrix))){
    colnames(designMatrix) <- paste0("t",seq_len(nTime))
  }
  
  # Obtain labels:
  nodeLabels <- rownames(designMatrix)
  allLabels <- c(designMatrix)

  ### Setup the filtering matrix ###
  isIncluded <- !is.na(allLabels)
  nIncluded <- sum(isIncluded)
  includedLabels <- allLabels[isIncluded]
  FilterMat <- matrix(0,nrow=nIncluded,ncol=length(allLabels))
  for (i in 1:nrow(designMatrix)){
    for (j in 1:ncol(designMatrix)){
      FilterMat[includedLabels == designMatrix[i,j], allLabels == designMatrix[i,j]]  <- 1
    }
  }
  
  
  ### Setup model matrices ###
  # filter matrix:
  MxFilter <- mxMatrix("Full",
                       nrow = nIncluded,
                       ncol = nNode * nTime,
                       values = FilterMat,
                       free = FALSE,
                       name = "F")
  
  ## Temporal effects:

  if ("temporal" %in% groupEqual){
    labs <- toLabel(beta[,,group],"beta",group=0)
  } else {
    labs <- toLabel(beta[,,group],"beta",group=group)
  }
  
  if (temporal){
    MxBeta <- mxMatrix("Full",
                       nrow=nNode,
                       ncol=nNode,
                       free = isFree(beta[,,group]),
                       labels = labs,
                       values = start("beta",startValues,ifelse(isFree(beta[,,group]),0,toFree(beta[,,group]))),
                       name = "Beta")
  } else {
    MxBeta <- mxMatrix("Full",
                       nrow=nNode,
                       ncol=nNode,
                       free = FALSE,
                       labels = labs,
                       values = 0,
                       name = "Beta")
  }

  
  ## Contemporaneous effects:
  # Model via kappa_zeta:
  if (!contemporaneous){
    MxKappa_zeta <-  mxMatrix("Symm",
                              nrow=nNode,
                              ncol=nNode,
                              free = FALSE,
                              values = 0,
                              name = "Kappa_zeta")
    
    
    MxSigma_zeta <-  mxMatrix("Symm",
                              nrow=nNode,
                              ncol=nNode,
                              free = FALSE,
                              values = 0,
                              name = "Sigma_zeta")  
    
  } else {
    if (!missing(kappa_zeta)){
      if ("contemporaneous" %in% groupEqual){
        labs <- toLabel(kappa_zeta[,,group],"kappa_zeta",symmetric = TRUE,group=0)
      } else {
        labs <- toLabel(kappa_zeta[,,group],"kappa_zeta",symmetric = TRUE,group=group)
      }
      
      MxKappa_zeta <-  mxMatrix("Symm",
                                nrow=nNode,
                                ncol=nNode,
                                free = isFree(kappa_zeta[,,group]),
                                labels = labs,
                                values = start("kappa_zeta",startValues,ifelse(isFree(kappa_zeta[,,group]),diag(nNode),toFree(kappa_zeta[,,group]))),
                                name = "Kappa_zeta")
      
      
      MxSigma_zeta <-  mxAlgebraFromString("solve(Kappa_zeta)",name = "Sigma_zeta")
    } else {
      if ("contemporaneous" %in% groupEqual){
        labs <- toLabel(sigma_zeta[,,group],"sigma_zeta",symmetric = TRUE,group=0)
      } else {
        labs <- toLabel(sigma_zeta[,,group],"sigma_zeta",symmetric = TRUE,group=group)
      }
      
      # Model via sigma_zeta
      MxSigma_zeta <-  mxMatrix("Symm",
                                nrow=nNode,
                                ncol=nNode,
                                free = isFree(sigma_zeta[,,group]),
                                labels = labs,
                                values = start("sigma_zeta",startValues,ifelse(isFree(sigma_zeta[,,group]),diag(nNode),toFree(sigma_zeta[,,group]))),
                                name = "Sigma_zeta")
      
      
      MxKappa_zeta <-  mxAlgebraFromString("solve(Sigma_zeta)",name = "Kappa_zeta")
    }
  }


  ## Between-subject effects:
  if (!betweenSubjects){
    MxKappa_mu <-  mxMatrix("Symm",
                            nrow=nNode,
                            ncol=nNode,
                            free = FALSE,
                            values = 0,
                            name = "Kappa_mu")    
    
    MxSigma_mu <-  mxMatrix("Symm",
                            nrow=nNode,
                            ncol=nNode,
                            free = FALSE,
                            values = 0,
                            name = "Sigma_mu")    
    
  } else {
    
    # Model via kappa_zeta:
    if (!missing(kappa_mu)){
      if ("betweenSubjects" %in% groupEqual){
        labs <- toLabel(kappa_mu[,,group],"kappa_mu",symmetric = TRUE,group=0)
      } else {
        labs <- toLabel(kappa_mu[,,group],"kappa_mu",symmetric = TRUE,group=group)
      }
      
      MxKappa_mu <-  mxMatrix("Symm",
                              nrow=nNode,
                              ncol=nNode,
                              free = isFree(kappa_mu[,,group]),
                              labels = labs,
                              values = start("kappa_mu",startValues,ifelse(isFree(kappa_mu[,,group]),diag(nNode),toFree(kappa_mu[,,group]))),
                              name = "Kappa_mu")
      
      
      MxSigma_mu <-  mxAlgebraFromString("solve(Kappa_mu)",name = "Sigma_mu")
    } else {
      if ("betweenSubjects" %in% groupEqual){
        labs <- toLabel(sigma_mu[,,group],"sigma_mu",symmetric = TRUE,group=0)
      } else {
        labs <- toLabel(sigma_mu[,,group],"sigma_mu",symmetric = TRUE,group=group)
      }
      
      # Model via sigma_mu
      MxSigma_mu <-  mxMatrix("Symm",
                              nrow=nNode,
                              ncol=nNode,
                              free = isFree(sigma_mu[,,group]),
                              labels = labs,
                              values = start("sigma_mu",startValues,ifelse(isFree(sigma_mu[,,group]),diag(nNode),toFree(sigma_mu[,,group]))),
                              name = "Sigma_mu")
      
      
      MxKappa_mu <-  mxAlgebraFromString("solve(Sigma_mu)",name = "Kappa_mu")
    }
    
  }

  if ("means" %in% groupEqual){
    labs <- toLabel(mu[,group,drop=FALSE],"mu",symmetric = TRUE,group=0)
  } else {
    labs <- toLabel(mu[,group,drop=FALSE],"mu",symmetric = TRUE,group=group)
  }
  ## Mean structure:
  MxMu <- mxMatrix("Full",ncol=1,nrow=nNode,
                   free = isFree(mu[,group,drop=FALSE]),
                   labels = labs,
                   values = start("mu",startValues,ifelse(isFree(mu[,group,drop=FALSE]),diag(nNode),toFree(mu[,group,drop=FALSE]))),
                   name="mu")
  
  # Implied mu for entire dataset:
  MxMuFull <- mxAlgebraFromString(paste0("F %*% rbind(",paste0(rep("mu",nTime),collapse=","),")"), name = "MuFull",dimnames = list(includedLabels,"Mean"))
  
  # Dummy matrix structures:
  Iden <-  mxMatrix("Diag",nrow=nNode,ncol=nNode,free=FALSE,
                    
                    # Starting values:
                    values = diag(nNode),
                    name = "Iden")
  
  Iden2 <-  mxMatrix("Diag",nrow=nNode^2,ncol=nNode^2,free=FALSE,
                     
                     # Starting values:
                     values = diag(nNode^2),
                     name = "Iden2")
  
  # Expected var-cov structure:
  Sigma_y_vec <- OpenMx::mxAlgebraFromString("solve(Iden2 - Beta %x% Beta) %*% cvectorize(Sigma_mu - Beta %*% Sigma_mu %*% t(Beta) + Sigma_zeta)",name = "Sigma_y_vec")

  # Create silly index columns:
  Cols <- lapply(seq_len(nNode), function(x){
    mxMatrix(type = "Full", nrow = nNode^2, ncol = 1, values = 1*(seq_len(nNode^2) %in% (1 + ((x-1)*nNode):(x*nNode-1))), free = FALSE, name = paste0("col",x))
  })
  
  # vec -> matrix:
  Sigma_y <- OpenMx::mxAlgebraFromString(paste0("cbind(",paste0("omxSelectRows(Sigma_y_vec,col",1:nNode,")",collapse=","),")"),name = "Sigma_y")
  
  # Expected lag 1 structure:
  Sigma_lag1 <- OpenMx::mxAlgebraFromString("Sigma_mu %*% t(Iden - Beta) + Sigma_y %*% t(Beta)",name = "Sigma_lag1")
  
  # Generate lag 2 - k:
  if (nTime > 2){
    Lags <- lapply(2:(nTime-1),function(k) {
      OpenMx::mxAlgebraFromString(paste0("Sigma_mu %*% t(Iden - Beta) + Sigma_lag",k-1, "%*% t(Beta)"),name = paste0("Sigma_lag",k))
    } )
  } else {
    Lags <- list()
  }
  
  # Create the entire implied variance-covariance matrix:
  diag= rep("y",nTime)
  offdiag = paste0("lag",c(unlist(sapply(rev(seq_len(nTime-1)),function(k)1:k))))
  m <- matrix(NA, ncol = nTime, nrow = nTime)
  m[lower.tri(m)] <- offdiag
  m[upper.tri(m)] <- t(m)[upper.tri(t(m))]
  diag(m) <- diag
  sigmat <- matrix(paste0("Sigma_",m),nTime,nTime)
  sigmat[upper.tri(sigmat)] <- paste0("t(",sigmat[upper.tri(sigmat)],")")
  str <- paste("F %*% cbind(",paste(apply(sigmat,1,function(x)paste0("rbind(",paste(x,collapse=","),")")),collapse=","),") %*% t(F)")
  Sigma <- OpenMx::mxAlgebraFromString(str, 
                                       name = "Sigma",dimnames = list(includedLabels,includedLabels)
  )
  
  # Expectation:
  logLik <- OpenMx::mxExpectationNormal(covariance = "Sigma",means = "MuFull")
  
  # Fit function:
  fitFunction <- OpenMx::mxFitFunctionML()
  
  # Data: covMat,
  openMxData <- OpenMx::mxData(covMat, type = "cov",means=means,numObs=sampleSize)
  
  # Combine model:
  allArgs <- list(name = name,
                  MxBeta,
                  # MxSigma_mu,MxSigma_mu_lower,
                  MxKappa_mu,MxSigma_mu,
                  MxSigma_zeta,MxKappa_zeta,
                  Sigma,
                  Sigma_y_vec, Sigma_y,#@col1,col2,col3,col4,col5,col6,
                  Sigma_lag1,
                  # Sigma_lag2,
                  logLik,
                  fitFunction,
                  openMxData,
                  Iden,Iden2,
                  MxMu,MxMuFull,
                  MxFilter)
  
  allArgs <- c(allArgs,Lags,Cols)
  
  FullModel <- do.call(OpenMx::mxModel,allArgs)
  
  # Return model:
  return(FullModel)
}