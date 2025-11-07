#' @title Meta-analysis for feature-level association testing for each covariate of interest
#'
#' @description
#' Implements the PALM framework for both single-study association analysis and
#' meta-analysis across multiple studies. This function conducts feature-level
#' association testing, evaluating one covariate of interest at a time.
#'
#' @author Zhoujingpeng Wei <zwei74@wisc.edu>
#'
#' @param rel.abd For a single association study, provide a matrix of relative abundance counts,
#' with sample IDs as row names and microbial feature IDs as column names.
#' The matrix must not contain any missing values.
#' For a meta-analysis, provide a list of such matrices, where each element corresponds
#' to one study, and the name of each element represents the study ID.
#'
#' @param covariate.interest For a single association study, provide a matrix where rows represent samples
#' and columns represent numeric covariates of interest. The order of samples must match the order in \code{rel.abd},
#' and the matrix must not contain any missing values.
#' For a meta-analysis, provide a list of such matrices, where each element corresponds to one study,
#' and the name of each element represents the study ID. The set of covariates may differ across studies.
#'
#' @param covariate.adjust For a single association study, provide a data frame where rows represent samples
#' and columns represent covariates to be adjusted for in the model. Covariates must be either factors or numeric vectors.
#' The order of samples must match the order in \code{rel.abd}, and the data frame must not contain missing values.
#' For a meta-analysis, provide a list of such data frames, where each element corresponds to one study,
#' and the name of each element represents the study ID. The set of covariates may differ across studies.
#' The default is \code{NULL}, meaning no covariates are adjusted.
#'
#' @param cluster For a single association study with correlated samples, provide a vector defining
#' the sample clusters. The order of samples in the vector must match the order in \code{rel.abd}.
#' For example, the values may represent subject IDs if each subject has multiple correlated samples
#' (e.g., in a longitudinal study).
#' For a meta-analysis, provide a list of such vectors, where each element pertains to one study with correlated samples,
#' and the name of each element represents the study ID. The order of samples in each vector must match the order in
#' \code{rel.abd} for the corresponding study.
#' The default is \code{NULL}, assuming all samples across all studies are independent.
#'
#' @param depth For a single association study, provide a vector representing the sequencing depth for each sample.
#' The length of the vector must match the number of samples, and the order must align with \code{rel.abd}.
#' For a meta-analysis, provide a list of such vectors, where each element corresponds to one study,
#' and the name of each element represents the study ID.
#' The default is \code{NULL}, in which case the sequencing depth is calculated as the row sums of \code{rel.abd}.
#'
#' @param depth.filter A cutoff value used to remove samples with sequencing depth less than or equal to the cutoff.
#' Default is \code{0}.
#'
#' @param prev.filter A cutoff value to remove microbial features with prevalence (proportion of nonzero observations)
#' less than or equal to the cutoff. This cutoff is applied to each study, so a feature may be removed in a subset of studies.
#' The cutoff value must be in the range \code{0–1}. Default is \code{0.1}.
#'
#' @param p.adjust.method Character string specifying the multiple testing correction method.
#' Options include \code{"holm"}, \code{"hochberg"}, \code{"hommel"}, \code{"bonferroni"}, \code{"BH"},
#' \code{"BY"}, \code{"fdr"}, and \code{"none"}. See \code{?stats::p.adjust} for more details.
#'
#' @param meta.method Character string specifying whether an equal-effects or random-effects model should be fitted.
#' An equal-effects model is fitted when \code{method = "EE"}.
#' A random-effects model is fitted by setting \code{method} to one of the following:
#' \code{"DL"}, \code{"HE"}, \code{"HS"}, \code{"HSk"}, \code{"SJ"}, \code{"ML"}, \code{"REML"},
#' \code{"EB"}, \code{"PM"}, \code{"GENQ"}, \code{"PMM"}, or \code{"GENQM"}.
#' The default is \code{"EE"}. See \code{?metafor::rma} for more details.
#'
#' @param correct Compositional effect correction method. The default is \code{"median"}.
#' If \code{"median"}, the compositional effect is corrected using the median of the estimates.
#' If \code{"tune"}, the compositional effect is corrected using a data-driven tuned value.
#' If \code{NULL}, no compositional correction is applied (i.e., the analysis uses relative-abundance summary
#' statistics directly and outputs RA-level summary statistics).
#' Recommendation: For sparse signals (<= 20% differential AA features), use \code{"median"}.
#' For moderately dense signals (> 20% and <= 50%), \code{"tune"} tends to be more accurate.
#' Both methods assume that the proportion of differential AA features is <= 50%.
#'
#' @return
#' A list containing one component for each covariate of interest.
#'
#' For single-study analysis, the component includes:
#' \item{palm_single}{A data frame containing association effect estimates, standard errors, p-values,
#' q-values, and summary statistics.}
#'
#' For meta-analysis across multiple studies, the component includes:
#' \item{palm_meta}{A data frame containing overall association effect estimates, standard errors,
#' p-values, and q-values for testing overall effects; p-values and q-values for testing cross-study
#' heterogeneity; and summary statistics (association effect estimates and standard errors) for individual studies.}
#'
#' @seealso
#' \code{\link{palm.get.summary}},
#' \code{\link{palm.null.model}},
#' \code{\link{palm.meta.summary}}
#'
#' @import dplyr
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
#' meta.result <- palm(rel.abd = rel.abd, covariate.interest = covariate.interest)
#' }

palm <- function(rel.abd,
                 covariate.interest,
                 covariate.adjust = NULL,
                 cluster = NULL,
                 depth = NULL,
                 depth.filter = 0,
                 prev.filter = 0.1,
                 p.adjust.method = "fdr",
                 meta.method = "EE",
                 correct = "median"
) {

  #=== Generate summary statistics ===#
  null.obj <- palm.null.model(rel.abd = rel.abd,
                              covariate.adjust = covariate.adjust,
                              depth = depth,
                              depth.filter = depth.filter,
                              prev.filter = prev.filter)

  summary.stats <- palm.get.summary(null.obj = null.obj,
                                    covariate.interest = covariate.interest,
                                    cluster = cluster,
                                    correct = correct)

  # Meta-analysis
  palm.model <- palm.meta.summary(summary.stats = summary.stats,
                                  p.adjust.method = p.adjust.method,
                                  meta.method = meta.method)

  return(palm.model)
}
