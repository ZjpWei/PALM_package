palm.get.multi.summary <- function(null.obj,
                                   covariate.interest,
                                   cluster = NULL,
                                   parallel.core = NULL,
                                   verbose = FALSE) {

  #=== Check input data ===#
  study.ID <- names(null.obj)
  feature.ID <- NULL
  cov.names <- NULL
  for(d in names(null.obj)){
    feature.ID <- c(feature.ID, colnames(null.obj[[d]]$Y_I))
    if(is.null(colnames(covariate.interest[[d]]))){
      stop("covariate.interest is not a matrix with column names. Please check your data.")
    }else{
      covariate.interest[[d]] <- data.frame(covariate.interest[[d]])
      cov.names <- c(cov.names, colnames(covariate.interest[[d]]))
    }
  }
  feature.ID <- unique(feature.ID)
  cov.names <- unique(cov.names)
  K <- length(feature.ID)
  if(verbose){
    message.num <- 0
  }

  ## Match rel.data, covariate.interest and covariate.adjust.
  SUB.id <- list()
  for(d in study.ID){
    Z <- cbind(null.obj[[d]]$Z[,-ncol(null.obj[[d]]$Z),drop=FALSE], as.matrix(covariate.interest[[d]]))
    if(all(rownames(Z) == rownames(covariate.interest[[d]]))){
      null.obj[[d]]$Z <- Z
    }else{
      stop("The rownames of covariate of inyterests don't match the input null.obj, please check your input data.")
    }

    if(is.null(covariate.interest[[d]])){
      stop("Study IDs in rel.data and covariate.interest don't match, please check the input data.")
    }
    if(!is.data.frame(covariate.interest[[d]])){
      stop("covariate.interest is not a list of data frames.\n")
    }
    if(nrow(covariate.interest[[d]]) != (nrow(null.obj[[d]]$Y_I) + length(null.obj[[d]]$rm.sample.idx))){
      stop("The sample size of covariate.interest is not correct, please check the input data.")
    }else{
      if(length(null.obj[[d]]$rm.sample.idx) > 0){
        covariate.interest[[d]] <- covariate.interest[[d]] %>% dplyr::slice(-null.obj[[d]]$rm.sample.idx)
      }
      rownames(covariate.interest[[d]]) <- rownames(null.obj[[d]]$Y_I)
    }
    if(is.null(cluster[[d]])){
      cluster.nm <- rownames(null.obj[[d]]$Y_I)
      names(cluster.nm) <- rownames(null.obj[[d]]$Y_I)
      SUB.id[[d]] <- cluster.nm
    }else{
      if(!is.vector(cluster[[d]])){
        stop("cluster is not a list of vectors. \n")
      }
      if(length(cluster[[d]]) != (nrow(null.obj[[d]]$Y_I) + length(null.obj[[d]]$rm.sample.idx))){
        stop("The sample size of cluster is not correct, please check the input data. \n")
      }else{
        if(length(null.obj[[d]]$rm.sample.idx) > 0){
          cluster[[d]] <- (cluster[[d]])[-null.obj[[d]]$rm.sample.idx]
        }
      }
      cluster.nm <- cluster[[d]]
      names(cluster.nm) <- rownames(null.obj[[d]]$Y_I)
      SUB.id[[d]] <- cluster.nm
    }
  }

  Sample.info <- list()
  Cov.int.info <- list()
  for(d in study.ID){
    Cov.int.info[[d]] <- colnames(covariate.interest[[d]])
    Sample.info[[d]] <- !is.na(rowSums(null.obj[[d]]$Z))
  }
  ## add comments
  ## Get summary statistics by palm_rcpp
  summary.stat.study <- palm_multi_rcpp(null_obj = null.obj, SUB_id = SUB.id,
                                        study_ID = study.ID, feature_ID = feature.ID,
                                        Cov_int_info = Cov.int.info,
                                        Sample_info = Sample.info)

  for(d in study.ID){
    Cov.int.info[[d]] <- colnames(covariate.interest[[d]])
    Sample.info[[d]] <- !is.na(rowSums(null.obj[[d]]$Z))
  }

  return(summary.stat.study)
}
