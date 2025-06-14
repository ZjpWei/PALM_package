#' @title Generate palm summary statistics for individual studies.
#'
#' @description This function directly takes the output object from the "palm.null.model" function and
#' constructs summary statistics for covariates of interest in each study.
#' The output can be directly used by the function "plam.meta.summary" to perform meta-analysis.
#'
#' @param null.obj The output of function "palm.null.model".
#' @param ... See function palm. In "covariate.interest" and "cluster", the order of samples should be matched
#' with the order in "rel.abd" used in the "palm.null.model" function.
#'
#' @return Output a list with each component for a study. The component includes the following elements.
#' \item{est}{A matrix contains the association effect estimates for the study. The row names are microbial feature IDs
#' and the column names are the covariates of interest IDs.}
#' \item{stderr}{A matrix contains the standard errors of the association effect estimates for the study. The row names
#' are microbial feature IDs and the column names are the covariates of interest IDs.}
#' \item{n}{Sample size for the study.}
#'
#' @seealso \code{\link{palm.null.model}},
#' \code{\link{palm.meta.summary}},
#' \code{\link{palm}}
#'
#' @import dplyr Rcpp
#' @export
#'
#' @examples
#' \donttest{
#' library("PALM")
#' data("CRC_data", package = "PALM")
#' CRC_abd <- CRC_data$CRC_abd
#' CRC_meta <- CRC_data$CRC_meta
#'
#' ########## Generate summary statistics ##########
#' rel.abd <- list()
#' covariate.interest <- list()
#' for(d in unique(CRC_meta$Study)){
#'   rel.abd[[d]] <- CRC_abd[CRC_meta$Sample_ID[CRC_meta$Study == d],]
#'   disease <- as.numeric(CRC_meta$Group[CRC_meta$Study == d] == "CRC")
#'   names(disease) <- CRC_meta$Sample_ID[CRC_meta$Study == d]
#'   covariate.interest[[d]] <- matrix(disease, ncol = 1, dimnames = list(names(disease), "disease"))
#' }
#'
#'  null.obj <- palm.null.model(rel.abd = rel.abd)
#'
#'  summary.stats <- palm.get.summary(null.obj = null.obj, covariate.interest = covariate.interest)
#' }
#'

palm.get.summary <- function(null.obj,
                             covariate.interest,
                             cluster = NULL,
                             correct = TRUE) {

  #=== Check input data ===#
  if(length(null.obj) == 1){
    if(!is.matrix(covariate.interest)){
      stop("null.obj only includes one study, covariate.interest should be a matrix.")
    }
    if(!is.null(cluster)){
      if(!is.vector(cluster)){
        stop("null.obj only includes one study, cluster should be a vector.")
      }
    }
    covariate.interest <- list("Study" = covariate.interest)
    cluster <- list("Study" = cluster)
  }

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

  ## Match rel.data, covariate.interest and covariate.adjust.
  SUB.id <- list()
  for(d in study.ID){
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
      SUB.id[[d]] <- as.numeric(factor(cluster.nm))
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
      SUB.id[[d]] <- as.numeric(factor(cluster.nm))
    }
  }

  Sample.info <- list()
  Cov.int.info <- list()
  for(d in study.ID){
    cov.int.nm <- colnames(covariate.interest[[d]])
    sample.ids <- rownames(null.obj[[d]]$Y_I)
    kp.sample.mat <- matrix(FALSE, nrow = length(sample.ids), ncol = length(cov.int.nm), dimnames = list(sample.ids, cov.int.nm))
    for(cov.name in cov.int.nm){
      kp.sample.mat[rownames(covariate.interest[[d]]),cov.name] <- TRUE
    }
    Sample.info[[d]] <- kp.sample.mat
    Cov.int.info[[d]] <- cov.int.nm
  }

  ## Get summary statistics by palm_rcpp
  summary.stat.study <- palm_rcpp(null_obj = null.obj, covariate_interest = covariate.interest,
                                  SUB_id = SUB.id,
                                  study_ID = study.ID, feature_ID = feature.ID, Cov_int_info = Cov.int.info,
                                  Sample_info = Sample.info)

  for(d in study.ID){
    summary.stat.study[[d]]$n <- nrow(covariate.interest[[d]])
  }

  ## Correct summary statistics
  if(correct){
    covariate.interest <- unique(unlist(lapply(summary.stat.study, function(d){colnames(d$est)})))
    for(cov.int in covariate.interest){
      for(d in study.ID){
        if(cov.int %in% colnames(summary.stat.study[[d]]$est)){
          ## adjust estimates
          min.delta <- median( - summary.stat.study[[d]]$est[,cov.int], na.rm = TRUE)
          non.na <- !is.na(summary.stat.study[[d]]$est[,cov.int])
          summary.stat.study[[d]]$est[non.na,cov.int] <- summary.stat.study[[d]]$est[non.na,cov.int] + min.delta

          ## adjust variances
          adjust.part <- sum(non.na) / (2 * sum(1 / sqrt(2*pi) / summary.stat.study[[d]]$stderr[non.na,cov.int]))^2
          summary.stat.study[[d]]$stderr[non.na,cov.int] <- sqrt((summary.stat.study[[d]]$stderr[non.na,cov.int])^2 + adjust.part)
        }
      }
    }
  }

  return(summary.stat.study)
}
