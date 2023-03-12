# devtools::document()

## Algorithm 1: Alternate Maximization
OrMIG_alternate_maximization <- function(Xmis, group, types, q,
                                         epsLogLike=1e-5, maxIter=10,
                                         verbose=FALSE, parallel=FALSE, lambda=1e-5){

  tmpList <- preprocess_data(Xmis, group, types)
  Xmis <- tmpList$Xmis
  group <- tmpList$group
  types <- tmpList$types

  resList <- gfmImpute(Xmis, group, types, q,
                       epsLogLike=epsLogLike, maxIter=maxIter,
                       verbose=verbose, parallel=parallel, lambda=lambda)

  # out <- try({resList <- gfmImpute(Xmis, group, types, q,
  #                                  epsLogLike=epsLogLike, maxIter=maxIter,
  #                                  verbose=verbose, parallel=parallel, lambda=lambda)},
  #            silent = TRUE)
  # if(class(out) == 'try-error'){
  #   message("For this number of factors q, MIG does not converge")
  #   message("The linear factor model-based method is used!")
  #
  #   resList <- LFM_impute(Xmis, q, group, types)
  #
  # }


  ## if include binary variable, use the imputation results of imputeFAMD
  hX_tmp <- resList$hX
  n <- nrow(hX_tmp)
  types_scale <- tmpList$type_scale
  ng <- length(types)
  hX <- tmpList$Xmiss
  hX[is.na(hX)] <- hX_tmp[is.na(hX)]
  for(s in 1:ng){
    switch (types[s],
            poisson = {
              id <- types_scale[[types[s]]][,1]
              hX[,id] <- hX_tmp[,id] * matrix(types_scale[[types[s]]][,2],
                                              n, length(id), byrow = TRUE)

            },
            gaussian = {
              id <- types_scale[[types[s]]][,1]
              hX[,id] <- hX_tmp[,id] * matrix(types_scale[[types[s]]][,2],
                                              n, length(id), byrow = TRUE)
            }
    )

  }

  return(return(list(hX=hX, fitList=resList)))
}
## Algotihm 2: Variational EM algorithm
OrMIG_VEM <- function(Xmis, group, types, q, offset=FALSE, epsELBO=1e-5, maxIter=30, verbose=TRUE, seed = 1){


  XmisList <- transferMat2List(Xmis, types,group)

  reslist <- OrMIG_vb.fit(XmisList, types= types, q=q, offset=offset, epsELBO=epsELBO, maxIter=maxIter,
                          verbose=verbose, seed = seed)

  return(reslist)
}

OrMIG <- function(...) UseMethod("OrMIG")

OrMIG.default <- function(Xmis, group, types, q, algorithm= c("VEM", "AM"),
                             epsLogLike=1e-5, maxIter=30, offset = FALSE,
                             verbose=FALSE, parallel=FALSE, lambda=1e-5, seed=1){

  algorithm <- match.arg(algorithm)

  if(algorithm == "VEM"){
    message('Starting the varitional EM(VEM) algorithm...\n')
    reslist <- OrMIG_VEM(Xmis, group, types= types, q=q, offset=offset, epsELBO=epsLogLike, maxIter=maxIter,
                         verbose=verbose, seed = seed)
    message('Finish the varitional EM(VEM) algorithm...\n')

  }else if(algorithm == "AM"){

    message('Starting the  alternate maximization(AM) algorithm...\n')
    reslist <- OrMIG_alternate_maximization(Xmis, group, types, q,
                                            epsLogLike=epsLogLike, maxIter=round(maxIter/2),
                                            verbose=verbose, parallel=parallel, lambda=lambda)
    message('Finish the alternate maximization(AM) algorithm...\n')

  }else{
    stop("OrMIG.default: unsupported algorithm!")
  }

  return(reslist)
}
OrMIG.matrix <- OrMIG.default

OrMIG.list <- function(XmisList, types, q, algorithm= c("VEM", "AM"),
                       epsLogLike=1e-5, maxIter=30, offset = FALSE,
                       verbose=FALSE, parallel=FALSE, lambda=1e-5, seed=1){


  if(length(XmisList) != length(types)) stop("OrMIG.list: The legnth of XmisList must be the same as that of types!")
  tmpList <- transferList2Mat(XmisList, types)

  Xmis <- tmpList$Xmis; group <- tmpList$group
  rm(tmpList)
  reslist <- OrMIG.default(Xmis, group, types, q, algorithm, epsLogLike=epsLogLike, maxIter=maxIter, offset = offset,
                           verbose=verbose, parallel=parallel, lambda=lambda, seed=seed)


  return(reslist)
}

OrMIG.data.frame <- function(DFmis,types, q, ...){

  ## handle the numeric names
  ns <- names(DFmis)

  if(sum(duplicated(ns))>0){
    message("There exists repeated colname in DFmis that will added 'XX' in front of raw names!")
    ns[duplicated(ns)] <- paste0("XX",ns[duplicated(ns)])
  }

  names(DFmis) <- ns
  ns_numeric <- as.numeric(ns)
  if(sum(!is.na(ns_numeric))> 0){
    message("There exists colname in DFmis that are numerical value such as '1', and it will replace with 'X1'!")
    ns_numeric[is.na(ns_numeric)] <- -1e20
    names(DFmis)[ns_numeric>-1e20] <- paste0("X", ns[ns_numeric>-1e20] )
  }


  var_class <- unique(sapply(DFmis, class))
  if(length(var_class) ==1 && var_class=='factor' && length(types)>1)
    stop("OrMIG: all variables in DFmis are factors, please set types='binomial'!")
  if(is.element('factor', var_class) && (!is.element('binomial', types)))
    stop('OrMIG: There must be "binomial" types since there exists "factor" in DFmis!')


  ## transfer data frame into matrix
  if(is.element('binomial', types)){
    Xmis <- model.matrix.lm(object = ~.+1, data=DFmis, na.action = "na.pass")
    Xmis <- Xmis[,-1] # remove the first column
  }else{
    Xmis <- sapply(DFmis, as.numeric)
    DFmis <- as.data.frame(Xmis)
  }

  # intialize group
  p <- ncol(Xmis)
  group <- rep(1, p)

  # discuss by three cases
  if(length(types) == 1){ # can directly run

    if(types == 'binomial'){
      var_class <- unique(sapply(DFmis, class))
      if("numeric" %in% var_class)
        stop("OrMIG: there are numeric variable in DFmis, please transfer it into factor if types='binomial'!")
      misList <- OrMIG.default(Xmis=Xmis, group=group, types=types, q=q, ...)
      levelList <- lapply(DFmis, levels)
      for(jname in names(levelList)){
        # message("jname = ", jname, '\n')
        DFmis[jname] <- imputeFactor(DFmis[jname], misList$hX)
      }
    }else if(types %in% c('gaussian','poisson') ){
      misList <- OrMIG.default(Xmis=Xmis, group=group, types=types, q=q, ...)
      ns <- names(DFmis)
      for(jname in ns){
        # message("jname = ", jname, '\n')
        DFmis[jname] <- imputeNumeric(DFmis[jname], misList$hX[,jname])
      }
    }else{
      stop(paste0("OrMIG: unsupported types: ", types, "!\n"))
    }

  }else if(length(types)==2){

    var_types <- apply(Xmis, 2, FindVartypes)
    ## transfor errorous count into continuous var.
    if(all(c("gaussian", "binomial") %in% types) && (! 'gaussian' %in% var_types)){
      var_types[var_types=="poisson"] <- "gaussian"
    }

    var_types <- factor(var_types)
    group <- as.numeric(var_types)
    types <- levels(var_types)
    misList <- OrMIG(Xmis, group = group, types=types, q=q, ...)

    levelList <- lapply(DFmis, levels)
    ns <- names(levelList)
    for(jname in ns){
      # message("jname = ", jname, '\n')
      if(is.null(levelList[[jname]])){ ## numeric
        DFmis[jname] <- imputeNumeric(DFmis[jname], misList$hX[,jname])
      }else{ ## factor
        DFmis[jname] <- imputeFactor(DFmis[jname], misList$hX) # each step can remove used columns
      }

    }
  }else if(length(types)==3){
    # obtain the variable types for each variable
    var_types <- apply(Xmis, 2, FindVartypes)
    var_types <- factor(var_types)
    group <- as.numeric(var_types)
    types <- levels(var_types)
    misList <- OrMIG(Xmis, group = group, types=types, q=q, ...)
    # misList <- MIG(Xmis, group = group, types=types, q=q, maxIter=3)

    levelList <- lapply(DFmis, levels)
    ns <- names(levelList)
    for(jname in ns){
      # message("jname = ", jname, '\n')
      if(is.null(levelList[[jname]])){ ## numeric
        DFmis[jname] <- imputeNumeric(DFmis[jname], misList$hX[,jname])
      }else{ ## factor
        DFmis[jname] <- imputeFactor(DFmis[jname], misList$hX) # each step can remove used columns
      }

    }
  }else{
    stop("The types has a length greater than 3!")
  }


 # misList$hX <- NULL
 # misList$DFimp <- DFmis
 return(list(DFimp=DFmis, used_types=types, gfmList = misList))
}

imputeFactor <- function(x, hX){
  if(!inherits(x[,1], 'factor')){stop("x must be a factor!")}
  id_na <- which(is.na(x))
  if(!length(id_na)) return(x)
  level <-  levels(x[,1])
  var_names <- paste0(names(x),  level)
  levelMat <- matrix(level[-1], nrow=length(id_na), ncol=length(level)-1, byrow=T)
  x[id_na,1] <- level[1]
  index_logical <- hX[id_na,var_names[-1]] == 1
  if(sum(index_logical)>0 && length(index_logical)<= length(id_na))
    x[id_na,1] <- levelMat[hX[id_na,var_names[-1]] == 1]
  return(x)
}

imputeNumeric <- function(x, hx){
  na_logical <- is.na(x[,1])
  x[na_logical,1] <- hx[na_logical]
  return(x)
}


FindVartypes <- function(x){
  xu <- unique(x)
  id_logical <- (!is.na(xu))
  if(sum(id_logical) <= 2){ # only two class
    y <- 'binomial'
  }else if(sum(id_logical) > 2 && all(xu[id_logical]>=0)&& all(abs(xu[id_logical]- round(xu[id_logical])) < 1e-20)){
    y <- 'poisson'
  }else{
    y <- 'gaussian'
  }
  return(y)
}
