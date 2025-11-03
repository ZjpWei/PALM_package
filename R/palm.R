#' @title Meta-analysis for feature-level association testing for each covariate of interest.
#'
#' @description Implementation of PALM for single-study association analysis and meta-analysis across
#' multiple studies. This function conducts feature-level association testing, evaluating one covariate
#' of interest at a time.
#'
#' @author Zhoujingpeng Wei <zwei74@wisc.edu>.
#'
#' @param rel.abd For a single association study, provide a matrix of relative abundance counts, with sample IDs as row names and microbial feature IDs as column names. The matrix must not contain any missing values.
#' For a meta-analysis, provide a list of such matrices, where each element corresponds to one study, and the name of each element represents the study ID.
#' @param covariate.interest For a single association study, provide a matrix where rows represent samples, and columns represent numeric covariates of interest. The order of samples must match the order in "rel.abd",
#' and the matrix must not contain any missing values. For a meta-analysis, provide a list of such matrices, where each element corresponds to one study, and the name of each element represents the study ID.
#' The set of covariates can differ across studies.
#' @param covariate.adjust For a single association study, provide a data frame where rows represent samples, and columns represent covariates to be adjusted in the model. Covariates must be either factors or numeric vectors.
#' The order of samples must match the order in rel.abd, and the data frame must not contain any missing values. For a meta-analysis, provide a list of such data frames, where each element corresponds to one study,
#' and the name of each element represents the study ID. The set of covariates can differ across studies. The default is NULL, meaning no covariates are adjusted.
#' @param cluster For a single association study with correlated samples, provide a vector that defines the sample clusters. The order of samples in the vector must match the order in "rel.abd". For example,
#' the values of this variable can be subject IDs if each subject has multiple correlated samples (e.g., in a longitudinal study). For a meta-analysis, provide a list of such vectors,
#' where each element pertains to one study with correlated samples, and the name of each element must be the study ID. The order of samples in each vector must match the order in "rel.abd" for the corresponding study.
#' The default is NULL, assuming all samples across all studies are independent, in which case this input is not required.
#' @param depth For a single association study, provide a vector representing the sequencing depth for each sample. The length of the vector must match the number of samples, and the order must align with "rel.abd".
#' For a meta-analysis, provide a list of such vectors, where each element corresponds to one study, and the name of each element must be the study ID. The length of each vector must match the sample size of the corresponding study.
#' The default is NULL, in which case the sequencing depth is calculated as the row sum of "rel.abd".
#' @param depth.filter A cutoff value to remove samples with sequencing depth less than or equal to the cutoff. Default is 0.
#' @param prev.filter A cutoff value remove microbial features with prevalence (proportion of nonzero observations)
#' less than or equal to the cutoff. This cutoff is applied to each study. So, a feature could be removed in a subset of the studies. The cutoff value must be in the range of 0-1. Default is 0.1.
#' @param p.adjust.method p-value correction method, a character string, should be one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none". Default is "fdr".
#' @param correct This argument takes the default value of "median". If "median", it will correct the compositional effect by median of estimates; If "tune", it will correct the compositional effect by tuned value;
#' If NULL, it will not correct the compositional effect (i.e., directly use RA summary statistic to do analysis and output RA-level summary statistics).
#'
#' @return Output a list with each component for a covariate of interest.
#'
#' For the single-study analysis, the component includes the following elements:
#' \item{palm_single}{A data framework containing association effect estimates, standard errors, p-values, and q-values
#' and summary statistics (association effect estimates and standard errors).}
#'
#' For the met-analysis of multiple association studies, the component includes the following elements:
#' \item{palm_meta}{A data framework containing overall association effect estimates, standard errors, p-values, and q-values for testing the overall effect;
#' p-values and q-values for testing cross-study effect heterogeneity; and summary statistics (association effect estimates and standard errors) for individual studies.}
#'
#' @seealso \code{\link{palm.get.summary}},
#' \code{\link{palm.null.model}},
#' \code{\link{palm.meta.summary}}
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
#' meta.result <- palm(rel.abd = rel.abd, covariate.interest = covariate.interest)
#' }
#'

palm <- function(rel.abd,
                 covariate.interest,
                 covariate.adjust = NULL,
                 cluster = NULL,
                 depth = NULL,
                 depth.filter = 0,
                 prev.filter = 0.1,
                 p.adjust.method = "fdr",
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
  palm.model <- palm.meta.summary(summary.stats = summary.stats, p.adjust.method = p.adjust.method)

  return(palm.model)
}
