# function to run a panelVAR model, or to stepwise search for one:
panelVAR <- function(
  data, # Data
  designMatrix, # this matrix takes the following form: rows stand for variables, columns for lags. Use NA if a var is missing ina  lag and rownames to indicate the variable labels
  groupVar, # Grouping variable, if missing it is added
  stepup = FALSE, # Should step-up model selection be used?
  MIcrit = 10, # Mod index criterion
  temporalModel = c("default","full","empty","manual"), # Models to use: default will pick, full = saturated, empty = empty, manual to take input form args below
  contemporaneousModel = c("default","full","empty","manual"),
  betweenSubjectsModel = c("default","full","empty","manual"),
  kappa_mu,
  sigma_mu,
  kappa_zeta,
  sigma_zeta,
  beta,
  mu,
  modIndices = stepup, # Set to FALSE if stepup = FALSE
  startValues = list(),
  name = "PanelVAR",
  searchMatrices = c("Beta","Kappa_mu","Kappa_zeta")
){
  if (stepUp & !modIndices){
    modIndices <- TRUE
  }
  
  ### Check input ###
  if (is.matrix(data)){
    data <- as.data.frame(data)
  }
  
  # models:
  temporalModel <- match.arg(temporalModel)
  if (temporalModel == "default"){
    if (!missing(beta)){
      temporalModel <- "manual"
    } else if (stepup){
      temporalModel <- "empty"
    } else {
      temporalModel <- "full"
    }
  }
  
  contemporaneousModel <- match.arg(contemporaneousModel)
  if (contemporaneousModel == "default"){
    if (!missing(kappa_zeta) || !missing(sigma_zeta)){
      contemporaneousModel <- "manual"
    } else if (stepup){
      contemporaneousModel <- "empty"
    } else {
      contemporaneousModel <- "full"
    }
  }
  
  betweenSubjectsModel <- match.arg(betweenSubjectsModel)
  if (betweenSubjectsModel == "default"){
    if (!missing(kappa_mu) || !missing(sigma_mu)){
      betweenSubjectsModel <- "manual"
    } else if (stepup){
      betweenSubjectsModel <- "empty"
    } else {
      betweenSubjectsModel <- "full"
    }
  }
  
  ### Setup groups ###
  if(missing(groupVar)){
    data$GROUP <- "singleGroup"
    groupVar <- "GROUP"
  }
  
  # check if group is number:
  data[[groupVar]] <- ifelse(is.na(as.numeric(data[[groupVar]])),data[[groupVar]],paste0("g",data[[groupVar]]))
  
  allGroups <- unique(data[[groupVar]])
  nGroup <- length(allGroups)
  
  ### Prepare the model matrices
  if (temporalModel == "empty"){
    beta <- array(0, c(nNode, nNode,nGroup))
  } else if (temporalModel == "full"){
    beta <- array(NA, c(nNode, nNode,nGroup))
  }
  
  if (contemporaneousModel == "empty"){
    kappa_zeta <- array(diag(NA, nNode),c(nNode,nNode,nGroup))
  } else if (contemporaneousModel == "full"){
    kappa_zeta <- array(NA,c(nNode,nNode,nGroup))
  }
  
  if (betweenSubjectsModel == "empty"){
    kappa_mu <- array(diag(NA, nNode),c(nNode,nNode,nGroup))
  } else if (betweenSubjectsModel == "full"){
    kappa_mu <- array(NA, c(nNode,nNode,nGroup))
  }
  
  if (missing(mu)){
    mu <- matrix(NA,nNode,1)
  }
  
  # Make arrays if needed:
  if (!missing(beta) && is.na(dim(beta)[3])){
    beta <- array(beta,c(dim(beta),nGroup))
  }
  if (!missing(kappa_mu) && is.na(dim(kappa_mu)[3])){
    kappa_mu <- array(kappa_mu,c(dim(kappa_mu),nGroup))
  }
  if (!missing(sigma_mu) && is.na(dim(sigma_mu)[3])){
    sigma_mu <- array(sigma_mu,c(dim(sigma_mu),nGroup))
  }
  if (!missing(kappa_zeta) && is.na(dim(kappa_zeta)[3])){
    kappa_zeta <- array(kappa_zeta,c(dim(kappa_zeta),nGroup))
  }
  if (!missing(sigma_zeta) && is.na(dim(sigma_zeta)[3])){
    sigma_zeta <- array(sigma_zeta,c(dim(sigma_zeta),nGroup))
  }
  if (!missing(mu) && (is.na(dim(mu)[2])) || (ncol(mu) == 1)){
    mu <- matrix(mu,nNode,nGroup)
  }

  ### Obtain the model ###
  GroupModels <- lapply(seq_len(nGroup),function(i){
    Model <- panelVAR_modelGen_stat(data[data[[groupVar]] == allGroups[i],], 
                                    nNode = nNode, nTime = nTime, nodeLabels = nodeLabels,allLabels = allLabels,
                                    designMatrix = designMatrix,
                                    kappa_mu = kappa_mu,
                                    sigma_mu = sigma_mu,
                                    kappa_zeta = kappa_zeta,
                                    sigma_zeta = sigma_zeta,
                                    mu = mu,
                                    beta = beta,
                                    startValues = list(),
                                    name = allGroups[i],
                                    group = i)
  })
  
  
  # Combine group models:
  allArgs <- c(
    list(name = name,mxFitFunctionMultigroup(allGroups)),
    GroupModels
  )
  FullModel <- do.call(mxModel,allArgs)
  
  # Run model:
  Fit <- mxRun(FullModel)
  
  # Compute MIs:
  if (modIndices){
    miRes <- mxMI(Fit$g1, matrices = searchMatrices)
  }
  
  # Model selection
  if (stepup){
    
  }
  browser()
  
  return(Fit)
}