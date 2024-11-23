#' @title Meta-analyze summary statistics across studies
#'
#' @description This function directly takes the summary statistics output from the "palm.get.summary"
#' function and combines the summary statistics across studies for selecting microbial signatures associated with each covariate of interest.
#'
#' @param summary.stats The output of function "palm.get.summary".
#' @param ... See function palm.
#'
#' @return Same output as the function "palm".
#'
#' @seealso \code{\link{palm.null.model}},
#' \code{\link{palm.get.summary}},
#' \code{\link{palm}}
#'
#' @import dplyr
#' @export
#'
#' @examples
#' \donttest{
#' library("PALM")
#' data("CRC_data", package = "miMeta")
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
#' null.obj <- palm.null.model(rel.abd = rel.abd, ref = "Coprococcus catus [ref_mOTU_v2_4874]")
#'
#' summary.stats <- palm.get.summary(null.obj = null.obj, covariate.interest = covariate.interest)
#'
#' ########## Meta-analysis ##########
#' meta.result <- palm.test(summary.stats = summary.stats)
#' }
#'

palm.test <- function(summary.stats = summary.stats,
                      p.adjust.method = "fdr"){


  data.type <- "RA"
  cor.method = "adjust"

  palm.meta <- list()
  covariate.interest <- unique(unlist(lapply(summary.stats, function(d){colnames(d$est)})))
  for(cov.int in covariate.interest){
    palm_fits <- list()
    study.ID <- NULL
    feature.ID <- NULL
    for(d in names(summary.stats)){
      if(cov.int %in% colnames(summary.stats[[d]]$est)){
        study.ID <- c(study.ID, d)
        feature.ID <- c(feature.ID, rownames(summary.stats[[d]]$est))
      }
    }
    feature.ID <- unique(feature.ID)

    AA.est <- matrix(NA, nrow = length(feature.ID), ncol = length(study.ID),
                     dimnames = list(feature.ID, as.character(study.ID)))
    AA.var <- matrix(NA, nrow = length(feature.ID), ncol = length(study.ID),
                     dimnames = list(feature.ID, as.character(study.ID)))

    for(d in study.ID){
      if(data.type == "AA"){
        min.delta <- 0
      }else if(data.type == "RA"){
        min.delta <- median(- summary.stats[[d]]$est[,cov.int], na.rm = TRUE)
      }else{
        stop("The data type should only be `AA` or `RA`, please check your input data.\n")
      }
      non.na <- !is.na(summary.stats[[d]]$est[,cov.int])
      median.var <- 0
      beta.coef <- summary.stats[[d]]$est[non.na,cov.int] + min.delta

      if(cor.method == "original"){
        pval <- 1 - pchisq(beta.coef^2 / summary.stats[[d]]$var[non.na,cov.int], df = 1)
        qval <- p.adjust(p = pval, method = p.adjust.method)
        palm_fits[[d]] <- data.frame(feature = names(beta.coef),
                                     coef = beta.coef,
                                     stderr = sqrt(summary.stats[[d]]$var[non.na,cov.int]),
                                     pval = pval,
                                     qval = qval)
      }else{
        pval <- 1 - pchisq(beta.coef^2 / (summary.stats[[d]]$var[non.na,cov.int] + sum(summary.stats[[d]]$var[non.na,cov.int])/sum(non.na)^2), df = 1)
        qval <- p.adjust(p = pval, method = p.adjust.method)
        palm_fits[[d]] <- data.frame(feature = names(beta.coef),
                                     coef = beta.coef,
                                     stderr = sqrt(summary.stats[[d]]$var[non.na,cov.int] + sum(summary.stats[[d]]$var[non.na,cov.int])/sum(non.na)^2),
                                     pval = pval,
                                     qval = qval)
      }


      AA.est[names(beta.coef),d] <- beta.coef
      if(cor.method == "original"){
        AA.var[names(beta.coef),d] <- summary.stats[[d]]$var[non.na,cov.int]
      }else{
        AA.var[names(beta.coef),d] <- summary.stats[[d]]$var[non.na,cov.int] + sum(summary.stats[[d]]$var[non.na,cov.int])/sum(non.na)^2
      }
    }

    if(length(study.ID) > 1){
      ## Meta statistics
      est.statics <- rowSums(AA.est / AA.var, na.rm = TRUE)
      var.statics <- rowSums(1 / AA.var, na.rm = TRUE)
      meta.coef <- est.statics / var.statics
      meta.var <- 1 / var.statics
      q.coef <- (meta.coef)^2 / meta.var

      pval.sin <- 1 - pchisq(q.coef, df = 1)
      qval.sin <- p.adjust(pval.sin, method = p.adjust.method)

      meta_fits <- data.frame(feature = rownames(AA.est),
                              coef = meta.coef,
                              stderr = sqrt(meta.var),
                              pval = pval.sin,
                              qval = qval.sin)

      palm.meta[[cov.int]] <- list(meta_fits = meta_fits, palm_fits = palm_fits)
    }else{
      palm.meta[[cov.int]] <- palm_fits[[study.ID]]
    }
  }
  return(palm.meta)
}
