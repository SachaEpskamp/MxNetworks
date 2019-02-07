# Generate stationary PanelVAR model (internally):
panelVAR_modelGen_saturated <- function(
  data, # Data
  designMatrix,
  nNode,
  nTime,
  name = "PanelVAR_saturated",
  group = 1,
  fixedSaturated = FALSE
){
 
  
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
  isIncluded <- !is.na(allLabels)
  nIncluded <- sum(isIncluded)
  includedLabels <- allLabels[isIncluded]
  
  ## Covariances:
  MxSigma <- mxMatrix("Symm",ncol=length(includedLabels),nrow=length(includedLabels),
                   free = !fixedSaturated, values = cov(data[,includedLabels],use="pairwise.complete.obs"),
                   name="Sigma", dimnames = list(includedLabels,includedLabels))
  

  ## Mean structure:
  MxMu <- mxMatrix("Full",ncol=1,nrow=length(includedLabels),
                   free = !fixedSaturated, values = colMeans(data[,includedLabels],na.rm = TRUE),
                   name="mu",dimnames = list(includedLabels,"mean"))
  
  # Expectation:
  logLik <- OpenMx::mxExpectationNormal(covariance = "Sigma",means = "mu")
  
  # Fit function:
  fitFunction <- OpenMx::mxFitFunctionML()
  
  # Data:
  openMxData <- OpenMx::mxData(data, type = "raw")
  
  # Combine model:
  allArgs <- list(name = name,
                  MxSigma,
                  logLik,
                  fitFunction,
                  openMxData,
                  MxMu)
  
  FullModel <- do.call(OpenMx::mxModel,allArgs)
  
  # Return model:
  return(FullModel)
}