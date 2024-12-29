#' @title Meta-analyze summary statistics across studies
#'
#' @description This function directly takes the summary statistics output from the "palm.get.summary"
#' function and combines the summary statistics across studies for feature-level association testing for each covariate of interest.
#'
#' @param summary.stats The output of function "palm.get.summary".
#' @param ... See function palm.
#'
#' @return Same output as the function "palm" when performing meta-analysis.
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
#' null.obj <- palm.null.model(rel.abd = rel.abd)
#'
#' summary.stats <- palm.get.summary(null.obj = null.obj, covariate.interest = covariate.interest)
#'
#' ########## Meta-analysis ##########
#' meta.result <- palm.meta.summary(summary.stats = summary.stats)
#' }
#'

palm.meta.summary <- function(summary.stats = summary.stats, p.adjust.method = "fdr"){

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
      non.na <- !is.na(summary.stats[[d]]$est[,cov.int])
      beta.coef <- summary.stats[[d]]$est[non.na,cov.int]
      var.coef <- (summary.stats[[d]]$stderr[non.na,cov.int])^2
      pval <- 1 - pchisq(beta.coef^2 / (var.coef + sum(var.coef)/sum(non.na)^2), df = 1)
      qval <- p.adjust(p = pval, method = p.adjust.method)
      palm_fits[[d]] <- data.frame(feature = names(beta.coef),
                                   coef = beta.coef,
                                   stderr = sqrt(var.coef + sum(var.coef)/sum(non.na)^2),
                                   pval = pval,
                                   qval = qval)

      AA.est[names(beta.coef),d] <- beta.coef
      AA.var[names(beta.coef),d] <- var.coef + sum(var.coef)/sum(non.na)^2
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

      ## Calculate Heterogeneity
      pval.het <- NULL
      for(k in 1:nrow(AA.est)){
        nonna.id <- !is.na(AA.est[k,])
        m <- try(metafor::rma(yi = AA.est[k,nonna.id], vi = AA.var[k,nonna.id],
                              control=list(maxiter=2000), method = "EE"), silent = TRUE)
        if(class(m)[1] != "try-error"){
          pval.het <- c(pval.het, m$QEp)
        }else{
          pval.het <- c(pval.het, NA)
        }
      }
      qval.het <- p.adjust(pval.het, method = p.adjust.method)

      ## Summary
      meta_fits <- data.frame(feature = rownames(AA.est),
                              coef = meta.coef,
                              stderr = sqrt(meta.var),
                              pval = pval.sin,
                              qval = qval.sin,
                              pval.het = pval.het,
                              qval.het = qval.het)

      for(d in study.ID){
        non.na <- !is.na(summary.stats[[d]]$est[,cov.int])
        meta_fits[[paste0(d, "_effect")]] <- summary.stats[[d]]$est[,cov.int]
        meta_fits[[paste0(d, "_stderr")]] <- sqrt((summary.stats[[d]]$stderr[,cov.int])^2 +
                                                    sum((summary.stats[[d]]$stderr[non.na,cov.int])^2)/sum(non.na)^2)
      }
      rownames(meta_fits) <- NULL
      palm.meta[[cov.int]] <- meta_fits
    }else{
      palm.meta[[cov.int]] <- palm_fits[[study.ID]]
    }
  }
  return(palm.meta)
}
