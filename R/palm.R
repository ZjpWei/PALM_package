#' @title Meta-analysis for selecting microbial signatures associated with covariate of interests.
#'
#' @description PALM is a meta-analysis method designed to account for the unique features of compositional microbiome data
#' in selecting microbial signatures.
#'
#' @author Zhoujingpeng Wei, Guanhua Chen, Zheng-Zheng Tang
#'
#' @param rel.abd A list of matrices for relative abundance counts. Each element of the list pertains to one study,
#' and the name of each element must be the study ID. Each matrix must have sample IDs as row names and microbial feature IDs as column names.
#' No missing value is allowed.
#' @param covariate.interest A list of matrices for single or multiple covariates of interest. Each element of the list pertains to one study, and the name of each element must be
#' the study ID. In a matrix of a given study, each row is a sample, and each column is a covariate of interest and can only be a numeric vector.
#' The set of covariates can be different among studies. In a given study, the order of samples should be matched with the order in "rel.abd". No missing value is allowed.
#' @param covariate.adjust A list of data frames for covariates that will be adjusted for in the model. Each element of the list pertains to one study,
#' and the name of each element must be the study ID. In a data frame of a given study, each row is a sample, and each column is a covariate to be adjusted and can only be a factor or numeric vector.
#' The set of covariates can be different among studies. In a given study, the order of samples should be matched with the order in "rel.abd". No missing value is allowed.
#' Default is NULL (not adjusting for any covariates)
#' @param cluster A list of vectors of variable that define the sample cluster for studies with correlated samples. Each element of the list pertains to one study with correlated samples,
#' and the name of each element must be the study ID. The order of samples should be matched with the order in "rel.abd" in the corresponding study.
#' For example, the values of this variable are subject IDs if each subject has multiple correlated samples
#' (e.g., measured in a longitudinal study). Default is NULL (all samples in all studies are independent, so no need to provide this input.).
#' @param depth A list of vectors for sequence depth. Each element of the list pertains to one study, the length of each vector must match the sample size of each study. For each element of vector, the value must be equal with
#' or larger than the sequence depth that are calculated from `covariate.interest`. Default is NULL, the calculated sequence depth will be used in this case.
#' @param depth.filter A cutoff value to remove samples with sequencing depth less than or equal to the cutoff. Default is 0.
#' @param prev.filter A cutoff value remove microbial features with prevalence (proportion of nonzero observations)
#' less than or equal to the cutoff. This cutoff is applied to each study. So, a feature could be removed in a subset of the studies. The cutoff value must be in the range of 0-1. Default is 0.1.
#' @param p.adjust.method p-value correction method, a character string, should be one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none". Default is "fdr".
#' @param verbose: Whether to print verbose information. Default is FALSE. (see details in Value)
#'
#' @return Output a list with each component for a covariate of interest.
#'
#' If only one study is in `rel.abd`, the component includes the following elements to deliver the single-study analysis.
#' \item{palm_fits}{A list of one data frame includes features, coefficients, standardized error, p-value and q-value.}
#'
#' If only more than one study are in `rel.abd`, the component includes the following elements to deliver the meta-analysis results and single-study results for each study.
#' \item{palm_fits}{A list of multiple data frames includes features, coefficients, standardized error, p-value and q-value.}
#' \item{palm_meta}{A data frame includes meta-analysis results with features, coefficients, standardized error, p-value and q-value.}
#'
#' @seealso \code{\link{palm.get.summary}},
#' \code{\link{palm.null.model}},
#' \code{\link{meta.summary}}
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
                 verbose = FALSE
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
                                    verbose = verbose)

  # Meta-analysis
  palm.model <- palm.meta.summary(summary.stats = summary.stats, p.adjust.method = p.adjust.method)

  return(palm.model)
}
