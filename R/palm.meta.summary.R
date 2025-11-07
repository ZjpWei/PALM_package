#' @title Meta-analyze summary statistics across studies
#'
#' @description
#' This function takes the summary statistics output from the \code{palm.get.summary} function
#' and combines the results across studies to perform feature-level association testing
#' for each covariate of interest.
#'
#' @param summary.stats The output object from the \code{palm.get.summary} function.
#' @param ... Additional arguments passed to \code{palm}.
#'
#' @return
#' An object with the same structure as the output of the \code{palm} function when performing meta-analysis.
#'
#' @seealso
#' \code{\link{palm.null.model}},
#' \code{\link{palm.get.summary}},
#' \code{\link{palm}}
#'
#' @import dplyr abess metafor
#' @export
#'
#' @examples
#' \donttest{
#' library(PALM)
#' data("CRC_data", package = "PALM")
#' CRC_abd <- CRC_data$CRC_abd
#' CRC_meta <- CRC_data$CRC_meta
#'
#' ########## Generate summary statistics ##########
#' rel.abd <- list()
#' covariate.interest <- list()
#' for (d in unique(CRC_meta$Study)) {
#'   rel.abd[[d]] <- CRC_abd[CRC_meta$Sample_ID[CRC_meta$Study == d], ]
#'   disease <- as.numeric(CRC_meta$Group[CRC_meta$Study == d] == "CRC")
#'   names(disease) <- CRC_meta$Sample_ID[CRC_meta$Study == d]
#'   covariate.interest[[d]] <- matrix(disease, ncol = 1,
#'                                     dimnames = list(names(disease), "disease"))
#' }
#'
#' null.obj <- palm.null.model(rel.abd = rel.abd)
#' summary.stats <- palm.get.summary(null.obj = null.obj,
#'                                   covariate.interest = covariate.interest)
#'
#' ########## Meta-analysis ##########
#' meta.result <- palm.meta.summary(summary.stats = summary.stats)
#' }

palm.meta.summary <- function(summary.stats,
                              p.adjust.method = "fdr",
                              meta.method = "EE"){

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
      AA.est[rownames(summary.stats[[d]]$est),d] <- summary.stats[[d]]$est[,cov.int]
      AA.var[rownames(summary.stats[[d]]$stderr),d] <- (summary.stats[[d]]$stderr[,cov.int])^2
    }

    if(length(study.ID) > 1){

      ## Calculate Heterogeneity
      meta_fits <- sapply(seq_len(nrow(AA.est)), function(i) {
        non.id <- !is.na(AA.est[i,])
        m <- try(metafor::rma(yi = AA.est[i,non.id], vi = AA.var[i,non.id], method = meta.method),
                 silent = TRUE)
        if(class(m)[1] != "try-error"){
          return(c(coef = m$beta, stderr = m$se, pval = m$QMp, pval.het = m$QEp))
        }else{
          return(c(coef = NA, stderr = NA, pval = NA, pval.het = NA))
        }
      })

      meta_fits <- data.frame(t(meta_fits)) %>% dplyr::transmute(
        feature = feature.ID,
        coef, stderr, pval,
        qval = p.adjust(pval, method = p.adjust.method),
        pval.het,
        qval.het = p.adjust(pval.het, method = p.adjust.method)
      )

      ## Summary
      if (sum(meta_fits$qval <= 0.05, na.rm = TRUE) >= nrow(AA.est) * 0.2) {
        paste0(
          "Over 20% of features are significant (FDR ≤ 0.05) for the covariate of interest '",
          cov.int,
          "'. Consider setting `correct = 'tune'` before performing the meta-analysis if you have not already done so."
        )
      } else if (sum(meta_fits$qval <= 0.05, na.rm = TRUE) >= nrow(AA.est) * 0.5) {
        paste0(
          "Over 50% of features are significant (FDR ≤ 0.05) for the covariate of interest '",
          cov.int,
          "'. This may violate model assumptions, and the results may not be reliable."
        )
      }

      for(d in study.ID){
        meta_fits[[paste0(d, "_effect")]] <- summary.stats[[d]]$est[,cov.int]
        meta_fits[[paste0(d, "_stderr")]] <- summary.stats[[d]]$stderr[,cov.int]
      }
      palm.meta[[cov.int]] <- meta_fits
    }else{
      beta.coef <- summary.stats[[study.ID]]$est[,cov.int]
      std.coef <- summary.stats[[study.ID]]$stderr[,cov.int]
      pval <- 1 - pchisq((beta.coef / std.coef)^2, df = 1)
      qval <- p.adjust(p = pval, method = p.adjust.method)
      palm_fits <- data.frame(feature = names(beta.coef),
                              coef = beta.coef,
                              stderr = std.coef,
                              pval = pval,
                              qval = qval)

      ## Summary
      if (sum(palm_fits$qval <= 0.05, na.rm = TRUE) >= nrow(AA.est) * 0.2) {
        paste0(
          "Over 20% of features are significant (FDR ≤ 0.05) for the covariate of interest '",
          cov.int,
          "'. Consider setting `correct = 'tune'` before performing the meta-analysis if you have not already done so."
        )
      } else if (sum(palm_fits$qval <= 0.05, na.rm = TRUE) >= nrow(AA.est) * 0.5) {
        paste0(
          "Over 50% of features are significant (FDR ≤ 0.05) for the covariate of interest '",
          cov.int,
          "'. This may violate model assumptions, and the results may not be reliable."
        )
      }
      palm.meta[[cov.int]] <- palm_fits
    }
  }
  return(palm.meta)
}
