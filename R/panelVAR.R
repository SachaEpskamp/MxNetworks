# function to run a panelVAR model, or to stepwise search for one:
panelVAR <- function(
  data, # Data
  designMatrix, # this matrix takes the following form: rows stand for variables, columns for lags. Use NA if a var is missing ina  lag and rownames to indicate the variable labels
  groupVar, # Grouping variable, if missing it is added
  # fixedSaturated = FALSE, # Just sets saturated model means and covs to fixed from sample covs (pairwise deletion). Not optimal... But fast!
  stepup = FALSE, # Should step-up model selection be used?
  speedStart = FALSE,
  jointStructure = FALSE, # Should the same structure (but not nessisarily parameters) be estimated?
  groupEqual = "none", # options: "networks", "means", "all", "temporal", "contemporaneous", "between-subjects"
  MIcrit = 10, # Mod index criterion
  temporalModel = c("default","diag","full","empty","manual"), # Models to use: default will pick, full = saturated, empty = empty, manual to take input form args below
  contemporaneousModel = c("default","full","empty","manual"),
  contemporaneousPartial, # Default to estimating GGM instead of cov matrix?
  betweenSubjectsModel = c("default","full","empty","manual"),
  betweenSubjectsPartial, # Default to estimating GGM instead of cov matrix?
  kappa_mu,
  sigma_mu,
  kappa_zeta,
  sigma_zeta,
  beta,
  mu,
  temporal = TRUE,
  contemporaneous = TRUE,
  betweenSubjects = TRUE,
  modIndices = stepup, # Set to FALSE if stepup = FALSE
  startValues = list(),
  name = "PanelVAR",
  searchMatrices = c("Beta","Kappa_mu","Kappa_zeta"),
  alpha = 0.01, # alpha used for modification
  verbose = TRUE,
  # saturatedModel,
  baselineModel,
  missing = "fiml",
  optimizeBIC = TRUE
){
  # optimizeBIC = FALSE # If TRUE, stepup search stops when BIC is not improved
  # Set now to FALSE for dummy
  
  # obtain from designMatrix:
  if (missing(designMatrix)){
    stop("'designMatrix' argument may not be missing. See ?panelVAR for more details.")
  }
  nNode <- nrow(designMatrix)
  nTime <- ncol(designMatrix)
  
  if (stepup & !modIndices){
    modIndices <- TRUE
  }
  
  # Defaults:
  if (missing(contemporaneousPartial)){
    if (missing(kappa_zeta) && missing(sigma_zeta)){
      contemporaneousPartial <- TRUE
    } else if (!missing(kappa_zeta) && missing(sigma_zeta)){
      contemporaneousPartial <- TRUE
    } else if (missing(kappa_zeta) && !missing(sigma_zeta)){
      contemporaneousPartial <- FALSE
    } else {
      stop("'kappa_zeta' and 'sigma_zeta' cannot both be missing.")
    }
  }
  
  if (missing(betweenSubjectsPartial)){
    if (missing(kappa_mu) && missing(sigma_mu)){
      betweenSubjectsPartial <- TRUE
    } else if (!missing(kappa_mu) && missing(sigma_mu)){
      betweenSubjectsPartial <- TRUE
    } else if (missing(kappa_mu) && !missing(sigma_mu)){
      betweenSubjectsPartial <- FALSE
    } else {
      stop("'kappa_mu' and 'sigma_mu' cannot both be missing.")
    }
  }
  
  ### Check input ###
  if (is.matrix(data)){
    data <- as.data.frame(data)
  }
  
  if (!all(groupEqual %in% c("none","networks", "means", "all", "temporal", "contemporaneous", "between-subjects"))){
    stop("'groupEqual' must be one or more of the following: 'none', 'all', 'networks', 'means', 'contemporaneous', 'temporal', 'between-subjects'")
  } 
  if (any(groupEqual == "all")){
    groupEqual <- c("means","temporal","contemporaneous","between-subjects")
  }
  if (any(groupEqual == "networks")){
    groupEqual <- c(groupEqual[groupEqual!="networks"],"temporal","contemporaneous","between-subjects")
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
  } else  if(temporalModel == "diag"){
    beta <- array(diag(NA,nNode), c(nNode, nNode,nGroup))
  }
  
  if (contemporaneousPartial){
    if (contemporaneousModel == "empty"){
      kappa_zeta <- array(diag(NA, nNode),c(nNode,nNode,nGroup))
    } else if (contemporaneousModel == "full"){
      kappa_zeta <- array(NA,c(nNode,nNode,nGroup))
    }
  } else {
    if (contemporaneousModel == "empty"){
      sigma_zeta <- array(diag(NA, nNode),c(nNode,nNode,nGroup))
    } else if (contemporaneousModel == "full"){
      sigma_zeta <- array(NA,c(nNode,nNode,nGroup))
    }
  }
  
  if (betweenSubjectsPartial){
    if (betweenSubjectsModel == "empty"){
      kappa_mu <- array(diag(NA, nNode),c(nNode,nNode,nGroup))
    } else if (betweenSubjectsModel == "full"){
      kappa_mu <- array(NA, c(nNode,nNode,nGroup))
    }
  } else {
    if (betweenSubjectsModel == "empty"){
      sigma_mu <- array(diag(NA, nNode),c(nNode,nNode,nGroup))
    } else if (betweenSubjectsModel == "full"){
      sigma_mu <- array(NA, c(nNode,nNode,nGroup))
    }
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
  
  # Estimate saturated model:
  ### This is taking too long... So I'll just use lavaan instead for this.
    # if (missing(saturatedModel)){
  #   GroupModels <- lapply(seq_len(nGroup),function(i){
  #     panelVAR_modelGen_saturated(data[data[[groupVar]] == allGroups[i],], 
  #      nNode = nNode, nTime = nTime, 
  #      designMatrix = designMatrix,
  #      name = paste0(allGroups[i],"_saturated"),
  #     group = i)
  #   })
  #   allArgs <- c(
  #     list(name = name,mxFitFunctionMultigroup(paste0(allGroups,"_saturated"))),
  #     GroupModels
  #   )
  #   FullModel <- do.call(mxModel,allArgs)
  #   
  #     
  #     if (verbose){
  #       message("Estimating saturated model...")
  #     }
  #   # Run model:
  #   saturatedModel <- mxRun(FullModel, silent = TRUE)
  # }
  
  # Estimate baseline model:
  # baseline model is simply a normal model with beta = 0 and sigma_zeta and sigma_mu diagonal.
  # equality constraints ARE included in this model. 
  if (verbose){
    message("Estimating saturated model...")
  }
  lavRes <- lavCor(data[,c(na.omit(c(designMatrix)), groupVar)], missing = missing, output = "fit", group = groupVar)
  sat_covs <- lavInspect(lavRes,"Sigma")
  sat_means <- lavInspect(lavRes,"Mu") 
  for (i in 1:nGroup){
    class(sat_covs[[i]]) <- "matrix"
    sat_means[[i]] <- as.vector(sat_means[[i]])
    names(sat_means[[i]]) <- colnames(sat_covs[[i]])
  }
  nObs <- lavInspect(lavRes, "nObs")
  
  
  if (missing(baselineModel)){
    GroupModels <- lapply(seq_len(nGroup),function(i){
      panelVAR_modelGen_stat(covMat = sat_covs[[i]],means = sat_means[[i]],sampleSize = nObs[i],
                             nNode = nNode, nTime = nTime, 
                             designMatrix = designMatrix,
                             sigma_mu = array(diag(NA,nNode),dim = c(nNode,nNode,nGroup)),
                             sigma_zeta = array(diag(NA,nNode),dim = c(nNode,nNode,nGroup)),
                             mu = matrix(NA,nNode,nGroup),
                             beta = array(0, dim=c(nNode, nNode,nGroup)),
                             startValues = list(),
                             name = paste0(allGroups[i],"_baseline"),
                             groupEqual = groupEqual,
                             group = i,
                             temporal = FALSE,
                             contemporaneous = TRUE,
                             betweenSubjects = FALSE)
    })
    allArgs <- c(
      list(name = name,mxFitFunctionMultigroup(paste0(allGroups,"_baseline"))),
      GroupModels
    )
    FullModel <- do.call(mxModel,allArgs)
    
    
    if (verbose){
      message("Estimating baseline model...")
    }
    # Run model:
    baselineModel <- mxRun(FullModel, silent = TRUE)
  }
  

  
  # Start of main loop (simply breaks if stepup = FALSE)
  run <- 1
  repeat{
    ### Obtain the model ###
    GroupModels <- lapply(seq_len(nGroup),function(i){
      panelVAR_modelGen_stat(covMat = sat_covs[[i]],means = sat_means[[i]],sampleSize = nObs[i],
                                      nNode = nNode, nTime = nTime, 
                                      designMatrix = designMatrix,
                                      kappa_mu = kappa_mu,
                                      sigma_mu = sigma_mu,
                                      kappa_zeta = kappa_zeta,
                                      sigma_zeta = sigma_zeta,
                                      mu = mu,
                                      beta = beta,
                                      startValues = list(),
                                      name = allGroups[i],
                                      groupEqual = groupEqual,
                                      group = i,
                                       temporal = temporal,
                                       contemporaneous = contemporaneous,
                                       betweenSubjects = betweenSubjects)
    })
    
    if (verbose){
      if (stepup){
        if (run == 1){
          message("Estimating starting model...")          
        } else {
          message("Estimating new model...")
        }
        
      } else {
        message("Estimating model...")
      }
    }
    # Combine group models:
    allArgs <- c(
      list(name = name,mxFitFunctionMultigroup(allGroups)),
      GroupModels
    )
    FullModel <- do.call(mxModel,allArgs)
    
    # Run model:
    Fit <- mxRun(FullModel, silent = TRUE)
    
    
    # Fit inds:
    
    # Model means:
    model_means <- lapply(allGroups,function(g)Fit[[g]][['MuFull']]$result)
    
    # Model covs:
    model_covs <- lapply(allGroups,function(g)Fit[[g]][['Sigma']]$result)
    
    # Baseline means:
    bas_means <- lapply(allGroups,function(g)baselineModel[[paste0(g,"_baseline")]][['MuFull']]$result)
    
    # Baseline Covs:
    bas_covs <- lapply(allGroups,function(g)baselineModel[[paste0(g,"_baseline")]][['Sigma']]$result)
    
    # Descriptives:
    nPar <- summary(Fit)[['estimatedParameters']]
    nSample <- summary(Fit)$numObs
    
    FitInds <- mxNetworkFit(
      model_means = model_means,
      model_covs = model_covs,
      sat_means = sat_means,
      sat_covs = sat_covs,
      bas_means = bas_means,
      bas_covs = bas_covs,
      nPar = nPar, # Number of parameters (total)
      bas_nPar = summary(baselineModel)[['estimatedParameters']],
      sampleSize = nSample, # Total sample size
      verbose = verbose,
      nGroup=nGroup
    )
    
    # Prune after speedstart?
    
    # if (speedStart && run == 2){
    #   message("Pruning all edges no longer significant at alpha / nGroup")
    #   
    #   names(summary(Fit))
    #   summary(Fit)$parameters
    #   
    # }
    
    # Check for BIC:
    if (optimizeBIC && run > 1){
      if (!FitInds$bic < oldFitInds$bic){
        if (verbose){
          message("New model did not improve BIC, returning previous model.")
          Fit <- oldFit
          FitInds <- oldFitInds
          break
        }
      }
    }
    
    # Compute MIs:
    if (modIndices){
      if (verbose){
        message("Computing modification indices...")
      }
      suppressMessages(miRes <- lapply(allGroups, function(g) mxMI(Fit[[g]], matrices = searchMatrices)))
      
      # Make list of MIs:
      miList <- lapply(miRes,"[[","MI")
      
 
    }
    
    # Stepup search:
    if (!stepup){
      break
    } else {
      
      # If any MI not significant, break:
      if (any(is.na(unlist(miList)))){
        if (verbose){
          message("NA modification index found. Stopping.")
        }
        break
      }
      
      if (!any(unlist(miList) > qchisq(alpha,1,lower.tail = FALSE))){
        if (verbose){
          message("No parameter can be added at given alpha level")
        }
        break
      }
      
      miPars <-  lapply(miRes,function(x)names(x$MI))
      
      if (speedStart && run == 1){
        message("Speedstart! Adding all MIs that are significant at alpha / nTest")
        }
      
      # Modify models:
      for (g in 1:nGroup){
        if (jointStructure){
          # Compute mean per MI (needed if jointStructure = TRUE):
          meanMIs <- rowMeans(do.call(cbind,miList))
          
          if (speedStart && run == 1){
            nTests <- length(unlist(miList))
            optimalMIs <- miPars[[g]][meanMIs > qchisq(alpha/nTests,1,lower.tail=FALSE)] 
            if (length(optimalMIs) == 0){
              optimalMIs <- miPars[[g]][which.max(meanMIs)]
            }
          } else {
            optimalMIs <- miPars[[g]][which.max(meanMIs)]
          }
        } else {
          # Test if any improves fit:
          if (any(is.na(miList[[g]]))){
            if (verbose){
              message("NA modification index found. Stopping.")
            }
            break
          }
          if (!any(miList[[g]] > qchisq(alpha,1,lower.tail = FALSE))){
            next
          }
          # Optimal parameter to add:
          if (speedStart && run == 1){

            nTests <- length(unlist(miList))
            optimalMIs <- miPars[[g]][miList[[g]] > qchisq(alpha/nTests,1,lower.tail=FALSE)] 
            if (length(optimalMIs) == 0){
              optimalMIs <- miPars[[g]][which.max(miList[[g]])] 
            }
          } else {
            optimalMIs <- miPars[[g]][which.max(miList[[g]])] 
          }
        }
        
        # Name of matrix:
        for (m in seq_along(optimalMIs)){
          optimalMI <- optimalMIs[m]
          matName <- gsub("\\_\\d*\\_\\d*\\_\\d*","",optimalMI)
          
          # Indices:
          inds <- as.numeric(unlist(regmatches(optimalMI,gregexpr("\\d+",optimalMI))))
          eval(parse(text=paste0(matName,"[",inds[1],",",inds[2],",",g,"] <- NA")))
          # Symmetric?
          if (grepl("(kappa)|(sigma)",optimalMI)){
            # Replace first two digits:
            eval(parse(text=paste0(matName,"[",inds[2],",",inds[1],",",g,"] <- NA")))
          }
        }

      }
    }
    
    run <- run + 1
    oldFit <- Fit
    oldFitInds <- FitInds
  }
  
  
  ### Compute fit indices  as per Lavaan ###
  # Saturated means:
  # sat_means <- lapply(allGroups,function(g)saturatedModel[[paste0(g,"_saturated")]][['mu']][['values']])
  
  # Saturated Covs:
  # sat_covs <- lapply(allGroups,function(g)saturatedModel[[paste0(g,"_saturated")]][['Sigma']][['values']])
  
  # Model means:
  model_means <- lapply(allGroups,function(g)Fit[[g]][['MuFull']]$result)
  
  # Model covs:
  model_covs <- lapply(allGroups,function(g)Fit[[g]][['Sigma']]$result)
  
  # Baseline means:
  bas_means <- lapply(allGroups,function(g)baselineModel[[paste0(g,"_baseline")]][['MuFull']]$result)
  
  # Baseline Covs:
  bas_covs <- lapply(allGroups,function(g)baselineModel[[paste0(g,"_baseline")]][['Sigma']]$result)
  
  # Descriptives:
  nPar <- summary(Fit)[['estimatedParameters']]
  nSample <- summary(Fit)$numObs
  
  FitInds <- mxNetworkFit(
    model_means = model_means,
    model_covs = model_covs,
    sat_means = sat_means,
    sat_covs = sat_covs,
    bas_means = bas_means,
    bas_covs = bas_covs,
    nPar = nPar, # Number of parameters (total)
    bas_nPar = summary(baselineModel)[['estimatedParameters']],
    sampleSize = nSample, # Total sample size
    verbose = verbose,
    nGroup=nGroup
  )
  

  # Store model matrices:
  modelMatrices <- list()
  if (!missing(beta)){
    modelMatrices[["beta"]] <- beta
  }
  if (!missing(kappa_zeta)){
    modelMatrices[["kappa_zeta"]] <- kappa_zeta
  }
  if (!missing(sigma_zeta)){
    modelMatrices[["sigma_zeta"]] <- sigma_zeta
  }
  if (!missing(kappa_mu)){
    modelMatrices[["kappa_mu"]] <- kappa_mu
  }
  if (!missing(sigma_mu)){
    modelMatrices[["sigma_mu"]] <- sigma_mu
  }
  
  # Store networks:
  # results <- list()
  PDC <- array(NA,c(nNode,nNode,nGroup)) # Partial directed correlations
  # DC  <- array(NA,c(nNode,nNode,nGroup)) # Directed correlations
  PCC <- array(NA,c(nNode,nNode,nGroup)) # partial contemporaneous correlations
  CC <-  array(NA,c(nNode,nNode,nGroup)) # Contemporaneous correlations
  PBC <- array(NA,c(nNode,nNode,nGroup)) # Partial between-subjects correlations
  BC <-  array(NA,c(nNode,nNode,nGroup)) # Between subjects correlations
  
  for (i in 1:nGroup){
    B <- Fit[[allGroups[i]]]$Beta$values
    if (contemporaneousPartial){
      Kz <- Fit[[allGroups[i]]]$Kappa_zeta$values
      Sz <- Fit[[allGroups[i]]]$Sigma_zeta$result
    } else {
      Kz <-Fit[[allGroups[i]]]$Kappa_zeta$results
      Sz <- Fit[[allGroups[i]]]$Sigma_zeta$values
    }
    if (betweenSubjectsPartial){
      Km <- Fit[[allGroups[i]]]$Kappa_mu$values
      Sm <- Fit[[allGroups[i]]]$Sigma_mu$result
    } else {
      Km <-Fit[[allGroups[i]]]$Kappa_mu$results
      Sm <- Fit[[allGroups[i]]]$Sigma_mu$values
    }

    PDC[,,i] <- computePDC(B,Kz)
    PCC[,,i] <- computePCC(Kz)
    PBC[,,i] <- computePCC(Km)
    CC[,,i] <- cov2cor(Sz)
    BC[,,i] <- cov2cor(Sm)
  }
  
  results <- list(
    PDC = PDC,
    PCC = PCC,
    PBC = PBC,
    CC = CC,
    BC = BC
  )
  
  Res <- list(
    modelFit = Fit,
    # saturated = saturatedModel,
    baseline = baselineModel,
    fitMeasures = FitInds,
    modelMatrices = modelMatrices,
    results = results
  )
  return(Res)
}
