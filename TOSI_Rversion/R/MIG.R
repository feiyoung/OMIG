# devtools::document()


OrMIG <- function(...) UseMethod("OrMIG")
OrMIG.default <- function(Xmis, group, type, q,
                             epsLogLike=1e-5, maxIter=10,
                             verbose=FALSE, parallel=FALSE, lambda=1e-5, ...){
  tmpList <- preprocess_data(Xmis, group, type)
  Xmis <- tmpList$Xmis
  group <- tmpList$group
  type <- tmpList$type
  out <- try({resList <- gfmImpute(Xmis, group, type, q,
                                   epsLogLike=epsLogLike, maxIter=maxIter,
                                   verbose=verbose, parallel=parallel, lambda=lambda)},
             silent = TRUE)
  if(class(out) == 'try-error'){
    message("For this number of factors q, MIG does not converge")
    message("The linear factor model-based method is used!")

    resList <- LFM_impute(Xmis, q, group, type)

  }


  ## if include binary variable, use the imputation results of imputeFAMD
  hX_tmp <- resList$hX
  n <- nrow(hX_tmp)
  type_scale <- tmpList$type_scale
  ng <- length(type)
  hX <- tmpList$Xmiss
  hX[is.na(hX)] <- hX_tmp[is.na(hX)]
  for(s in 1:ng){
    switch (type[s],
            poisson = {
              id <- type_scale[[type[s]]][,1]
              hX[,id] <- hX_tmp[,id] * matrix(type_scale[[type[s]]][,2],
                                                  n, length(id), byrow = TRUE)

            },
            gaussian = {
              id <- type_scale[[type[s]]][,1]
              hX[,id] <- hX_tmp[,id] * matrix(type_scale[[type[s]]][,2],
                                                  n, length(id), byrow = TRUE)
            }
    )

  }


 return(return(list(hX=hX, gfmList=resList)))
}
OrMIG.matrix <- OrMIG.default
OrMIG.data.frame <- function(DFmis,type, q, ...){

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
  if(length(var_class) ==1 && var_class=='factor' && length(type)>1)
    stop("OrMIG: all variables in DFmis are factors, please set type='binomial'!")
  if(is.element('factor', var_class) && (!is.element('binomial', type)))
    stop('OrMIG: There must be "binomial" type since there exists "factor" in DFmis!')


  ## transfer data frame into matrix
  if(is.element('binomial', type)){
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
  if(length(type) == 1){ # can directly run

    if(type == 'binomial'){
      var_class <- unique(sapply(DFmis, class))
      if("numeric" %in% var_class)
        stop("OrMIG: there are numeric variable in DFmis, please transfer it into factor if type='binomial'!")
      misList <- OrMIG.default(Xmis=Xmis, group=group, type=type, q=q, ...)
      levelList <- lapply(DFmis, levels)
      for(jname in names(levelList)){
        # message("jname = ", jname, '\n')
        DFmis[jname] <- imputeFactor(DFmis[jname], misList$hX)
      }
    }else if(type %in% c('gaussian','poisson') ){
      misList <- OrMIG.default(Xmis=Xmis, group=group, type=type, q=q, ...)
      ns <- names(DFmis)
      for(jname in ns){
        # message("jname = ", jname, '\n')
        DFmis[jname] <- imputeNumeric(DFmis[jname], misList$hX[,jname])
      }
    }else{
      stop(paste0("OrMIG: unsupported type: ", type, "!\n"))
    }

  }else if(length(type)==2){

    var_type <- apply(Xmis, 2, FindVarType)
    ## transfor errorous count into continuous var.
    if(all(c("gaussian", "binomial") %in% type) && (! 'gaussian' %in% var_type)){
      var_type[var_type=="poisson"] <- "gaussian"
    }

    var_type <- factor(var_type)
    group <- as.numeric(var_type)
    type <- levels(var_type)
    misList <- OrMIG(Xmis, group = group, type=type, q=q, ...)

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
  }else if(length(type)==3){
    # obtain the variable type for each variable
    var_type <- apply(Xmis, 2, FindVarType)
    var_type <- factor(var_type)
    group <- as.numeric(var_type)
    type <- levels(var_type)
    misList <- OrMIG(Xmis, group = group, type=type, q=q, ...)
    # misList <- MIG(Xmis, group = group, type=type, q=q, maxIter=3)

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
    stop("The type has a length greater than 3!")
  }


 # misList$hX <- NULL
 # misList$DFimp <- DFmis
 return(list(DFimp=DFmis, used_type=type, gfmList = misList))
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


FindVarType <- function(x){
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
