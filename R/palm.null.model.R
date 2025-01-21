#' @title Get information from the null model.
#'
#' @description The "palm.null.model" function computes the estimated mean microbial feature proportions and residuals
#' of the relative abundance under the null model of no associations for each study (so this function does not need the
#' covariate.interest information). The output of this function will be fed into the "palm.get.summary" function to
#' construct summary statistics for covariates of interest in each study.
#'
#' @param ... See function palm
#'
#' @return Output a list with each component for a study. The component includes the following elements.
#' \item{Y_I}{A matrix of predicted abundance.}
#' \item{Y_R}{A matrix of Firth's bias corrected abundance residual.}
#' \item{Z}{Cleaned ad formatted covariate.adjust.}
#' \item{rm.sample.idx}{The index of the removed samples for the study.}
#'
#' @seealso \code{\link{palm.get.summary}},
#' \code{\link{palm.meta.summary}},
#' \code{\link{palm}}
#'
#' @import dplyr brglm2
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
#' for(d in unique(CRC_meta$Study)){
#'   rel.abd[[d]] <- CRC_abd[CRC_meta$Sample_ID[CRC_meta$Study == d],]
#' }
#'
#' null.obj <- palm.null.model(rel.abd = rel.abd)
#' }
#'

palm.null.model <- function(rel.abd,
                            covariate.adjust = NULL,
                            depth = NULL,
                            depth.filter = 0,
                            prev.filter = 0.1){

  #=== Check input data ===#
  if(is.matrix(rel.abd)){
    rel.abd <- list("Study" = rel.abd)
    if(!is.null(covariate.adjust)){
      if(!is.data.frame(covariate.adjust)){
        stop("rel.abd is a matrix, covariate.adjust should be a data frame.")
      }
    }
    if(!is.null(depth)){
      if(!is.vector(depth)){
        stop("rel.abd is a matrix, depth should be a vector.")
      }
    }
    covariate.adjust <- list("Study" = covariate.adjust)
    depth <- list("Study" = depth)
  }

  study.ID <- names(rel.abd)
  if(is.null(study.ID)){
    stop("Please check the study name in rel.data.\n")
  }else{
    if(length(study.ID) != length(unique(study.ID))){
      stop("Multiple studies has the same name in rel.data, please check the input data.\n")
    }
  }

  #=== match rel.data and covariate.adjust ===#
  for(d in study.ID){
    if(!all(!is.na(rel.abd[[d]]))){
      stop("Detect NA in relative abundant counts.\n")
    }
    if(min(rel.abd[[d]]) < 0){
      stop("Detect negative value in relative abundant counts.\n")
    }
    if(!is.null(covariate.adjust[[d]])){
      if(!is.data.frame(covariate.adjust[[d]])){
        stop("covariate.adjust is not a list of data frames.\n")
      }
      if(nrow(covariate.adjust[[d]]) != nrow(rel.abd[[d]])){
        stop("The sample size doesn't match between rel.abd and covariate.adjust, please check the input data.\n")
      }else{
        rownames(covariate.adjust[[d]]) <- rownames(rel.abd[[d]])
      }
      if(!all(!is.na(covariate.adjust[[d]]))){
        stop("Detect NA in covariate.adjust.\n")
      }
    }
    if(is.null(depth[[d]])){
      depth[[d]] <- rowSums(rel.abd[[d]])
    }else{
      if(!is.vector(depth[[d]])){
        stop(paste0("depth in study ", d, " is not a vector, please check the input data.\n"))
      }
      if(length(depth[[d]]) != nrow(rel.abd[[d]])){
        stop(paste0("The length of depth in study ", d, " doesn't match the dimension of input data.\n"))
      }
    }
  }

  #=== Filter the samples using depth.filter ===#
  rm.sample.idx <- list()
  for(d in study.ID){
    depth.kp <- which(depth[[d]] > depth.filter)
    rm.sample.idx[[d]] <- which(depth[[d]] <= depth.filter)
    rel.abd[[d]] <- as.matrix(rel.abd[[d]][depth.kp,,drop=FALSE])
    depth[[d]] <- depth[[d]][depth.kp]
    if(!is.null(covariate.adjust[[d]])){
      covariate.adjust[[d]] <- covariate.adjust[[d]] %>% dplyr::slice(depth.kp)
    }
  }

  #=== Match samples in relative abundant counts and sample data ===#
  data.relative <- list()
  for(d in study.ID){
    Y.pool <- rel.abd[[d]]

    if(nrow(Y.pool) < 20){
      warning(paste0("Less than 20 samples in study ", d,
                     ", the summary statistics maybe not stable.\n"))
    }

    X.pool <- matrix(NA, nrow = nrow(Y.pool))
    if(!is.null(covariate.adjust)){
      cov.adjust <- covariate.adjust[[d]]
      for(cov_name in colnames(cov.adjust)){
        if(all(!is.na(cov.adjust[[cov_name]]))){
          if(is.factor(cov.adjust[[cov_name]]) | is.character(cov.adjust[[cov_name]])){
            dummys <- as.data.frame(model.matrix(formula(paste("~", cov_name)), data = cov.adjust[cov_name]))
            X.pool <- as.matrix(cbind(dummys[,-1], X.pool))
          }else{
            X.pool <- cbind(cov.adjust[[cov_name]], X.pool)
          }
        }else{
          stop(paste0("NA presents in `", cov_name, "`, please check the input data.\n"))
        }
      }
    }
    data.relative[[d]] <- list(Y = Y.pool, X = X.pool, N = depth[[d]])
  }

  #=== Generate null model ===#
  reg.fit <- list()
  for(d in study.ID){
    Y.sub <- data.relative[[d]]$Y
    X.sub <- cbind(1, data.relative[[d]]$X)
    N.sub <- data.relative[[d]]$N
    colnames(X.sub) <- c("Intercept", paste0("V_", as.character(1:(ncol(X.sub)-1))))
    rownames(X.sub) <- rownames(Y.sub)

    ## Apply prevalence filter for each study.
    feature.set.filter <- colMeans(Y.sub != 0) > prev.filter
    Y.sub <- Y.sub[,feature.set.filter,drop=FALSE]
    feature.ids <- colnames(Y.sub)

    if(ncol(X.sub) == 2){
      est.single <- matrix(NA, ncol = 2, nrow = ncol(Y.sub),
                           dimnames = list(feature.ids, c(":(Intercept)", ":X")))
      est.single[,":X"] <- 0
    }else if(ncol(X.sub) == 3){
      est.single <- matrix(NA, ncol = 3, nrow = ncol(Y.sub),
                           dimnames = list(feature.ids, c(":(Intercept)", ":X", ":V")))
      est.single[,":V"] <- 0
    }else{
      est.single <- matrix(NA, ncol = ncol(X.sub), nrow = ncol(Y.sub),
                           dimnames = list(feature.ids, c(":(Intercept)", paste0(":XV_", 1:(ncol(X.sub)-1)))))
      est.single[,paste0(":XV_", ncol(X.sub)-1)] <- 0
    }

    Y_b <- matrix(0, nrow = nrow(Y.sub), ncol = ncol(Y.sub), dimnames = list(rownames(Y.sub), colnames(Y.sub)))
    for(k in 1:ncol(Y.sub)){
      # Try brglmFit model
      suppressWarnings(
        if(ncol(X.sub) == 2){
          input.data.tmp = list(Y=Y.sub[,k])
          glm.out.tmp =  try(glm(Y ~ 1, data = input.data.tmp,
                                 family = poisson(link = "log"),
                                 offset = log(N.sub),
                                 method = brglm2::brglmFit,
                                 type = "AS_mean"),
                             silent = TRUE)
        }else{
          input.data.tmp = list(Y=Y.sub[,k], X = X.sub[,-c(1, ncol(X.sub))])
          glm.out.tmp =  try(glm(Y ~ X, data = input.data.tmp,
                                 family = poisson(link = "log"),
                                 offset = log(N.sub),
                                 method = brglm2::brglmFit,
                                 type = "AS_mean"),
                             silent = TRUE)
        }
      )
      if(class(glm.out.tmp)[1] != "try-error"){
        if(glm.out.tmp$converged){
          names(glm.out.tmp$coefficients) <- paste0(":", names(glm.out.tmp$coefficients))
          est.single[k,names(glm.out.tmp$coefficients)] <- glm.out.tmp$coefficients

          ## compute the bias part
          Qmat <- qr.Q(glm.out.tmp$qr)
          Y_b[,k] <- 0.5 * rowSums(Qmat * Qmat)
        }else{
          warning("Cannot converge for feature ", feature.ids[k],
                  ", remove this taxa in nul model.\n")
          est.single[k,] <- NA
        }
      }else{
        warning("Cannot generate score estimate for feature ", feature.ids[k],
                ", remove this taxa in null model.\n")
        est.single[k,] <- NA
      }
    }
    ## Remove NA estimates
    non_na_taxa <- which(!is.na(rowSums(est.single)))
    est.single <- est.single[non_na_taxa,,drop=FALSE]
    Y.sub <- Y.sub[,non_na_taxa,drop=FALSE]
    Y_b <- Y_b[,non_na_taxa,drop=FALSE]

    #=== Summarize null model output ===#
    Y_R <- NULL
    Y_I <- NULL
    n.rows <- nrow(est.single)
    n.cols <- ncol(X.sub)
    for(i in 1:length(N.sub)){
      dd <- colSums(matrix(rep(X.sub[i,], n.rows), nrow = n.cols) * t(est.single),
                    na.rm = TRUE)
      Y_I.i <- exp(dd) * N.sub[i]
      Y_I <- rbind(Y_I, Y_I.i)
      Y_R <- rbind(Y_R, Y.sub[i,] - Y_I.i)
    }
    Y_R <- Y_R + Y_b
    rownames(Y_I) <- rownames(Y.sub)
    rownames(Y_R) <- rownames(Y.sub)
    reg.fit.one <- list(Y_I = Y_I, Y_R = Y_R, Z = X.sub, rm.sample.idx = rm.sample.idx[[d]])

    #=== output ===#
    reg.fit[[d]] <- reg.fit.one
  }

  return(reg.fit)
}
