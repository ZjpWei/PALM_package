// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include <unordered_map>
#include <string>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]

using Rcpp::List;
using Rcpp::NumericMatrix;
using Rcpp::NumericVector;
using Rcpp::CharacterVector;
using Rcpp::LogicalMatrix;
using Rcpp::LogicalVector;
using Rcpp::IntegerVector;
using Rcpp::Named;
using Rcpp::as;
using Rcpp::wrap;
using Rcpp::stop;

// ---- Small utilities ---------------------------------------------------------

// Convert a LogicalVector mask to arma::uvec of 0-based row indices where TRUE.
static inline arma::uvec mask_to_uvec(const LogicalVector &mask) {
  std::vector<arma::uword> idx;
  idx.reserve(mask.size());
  for (int i = 0; i < mask.size(); ++i) if (mask[i]) idx.push_back(static_cast<arma::uword>(i));
  return arma::uvec(idx);
}

// Build a fast string->row index map for feature_ID.
static inline std::unordered_map<std::string, int>
  build_feature_index_map(const CharacterVector &feature_ID) {
    std::unordered_map<std::string, int> mp;
    mp.reserve(feature_ID.size() * 2);
    for (int i = 0; i < feature_ID.size(); ++i) {
      mp.emplace(Rcpp::as<std::string>(feature_ID[i]), i); // 0-based row index
    }
    return mp;
  }

// Fill a NumericMatrix with NA
static inline void fill_na(NumericMatrix &M) {
  std::fill(M.begin(), M.end(), Rcpp::NumericVector::get_na());
}

// Safe solve with SPD hint; falls back to generic solve if needed.
static inline arma::mat solve_spd(const arma::mat &A, const arma::mat &B) {
  arma::mat X;
  bool ok = arma::solve(X, A, B, arma::solve_opts::fast + arma::solve_opts::likely_sympd);
  if (!ok) ok = arma::solve(X, A, B);
  if (!ok) stop("Linear solve failed (matrix may be singular/ill-conditioned).");
  return X;
}

// ---- Core per-feature computation -------------------------------------------
static inline void estimate_one_feature(
    const arma::mat &X,
    const arma::vec &yR,
    const arma::vec &yI,
    const arma::uword G,
    const arma::uvec &gidx,
    double &est_out,
    double &var_out
) {
  const arma::uword n = X.n_rows;
  const arma::uword p = X.n_cols;

  // Scores per subject: S_i = x_i * yR_i   (n x p)
  arma::mat S = X.each_col() % yR;

  // Hessian: H = X' diag(yI) X = (X .* yI)^T * X   (p x p)
  arma::mat H = (X.each_col() % yI).t() * X;

  if (p == 0) stop("X has zero columns.");
  const arma::uword beta_col = p - 1;

  arma::uvec idx_beta(1); idx_beta[0] = beta_col;
  arma::uvec idx_gamm;
  if (p > 1) idx_gamm = arma::regspace<arma::uvec>(0, beta_col - 1);

  const arma::mat Hbb = H(idx_beta, idx_beta);                   // 1x1
  arma::mat I1;
  if (p == 1) {
    I1 = Hbb;
  } else {
    const arma::mat Hbg = H(idx_beta, idx_gamm);                // 1 x (p-1)
    const arma::mat Hgg = H(idx_gamm, idx_gamm);                // (p-1) x (p-1)
    I1 = Hbb - Hbg * solve_spd(Hgg, Hbg.t());                   // Schur complement
  }

  const arma::rowvec ssum = arma::sum(S, 0);                    // 1 x p
  const double ssum_beta = ssum[beta_col];

  arma::mat I1inv = solve_spd(I1, arma::eye(I1.n_rows, I1.n_cols)); // 1x1 inverse via solve
  const double est_beta = (I1inv(0,0) * ssum_beta);

  arma::vec U_beta(n, arma::fill::zeros);                       // per-observation adj. score

  if (p == 1) {
    U_beta = S.col(beta_col);
  } else {
    const arma::mat S_gamm = S.cols(idx_gamm);                  // n x (p-1)
    const arma::mat Hgg     = H(idx_gamm, idx_gamm);

    const arma::mat Hgg_inv_SgT = solve_spd(Hgg, S_gamm.t());   // (p-1) x n
    const arma::mat Hbg         = H(idx_beta, idx_gamm);        // 1 x (p-1)

    const arma::vec adjust = (Hbg * Hgg_inv_SgT).t();           // n x 1
    U_beta = S.col(beta_col) - adjust;
  }

  arma::vec U_combined(G, arma::fill::zeros);
  for (arma::uword i = 0; i < n; ++i) {
    const arma::uword g = gidx[i];
    if (g >= G) stop("Cluster index out of range.");
    U_combined[g] += U_beta[i];
  }

  const double meat = arma::as_scalar(U_combined.t() * U_combined);
  const double var_beta = I1inv(0,0) * meat * I1inv(0,0);

  est_out = est_beta;
  var_out = var_beta;
}

// ---- Main exported function --------------------------------------------------

// [[Rcpp::export]]
Rcpp::List palm_rcpp(
    Rcpp::List null_obj,
    Rcpp::List covariate_interest,
    Rcpp::List SUB_id,
    Rcpp::CharacterVector study_ID,
    Rcpp::CharacterVector feature_ID,
    Rcpp::List Cov_int_info,
    Rcpp::List Sample_info
) {
  const int K = feature_ID.size();
  auto feat_map = build_feature_index_map(feature_ID); // name -> row index

  Rcpp::List summary_stat_study;

  for (int idx_d = 0; idx_d < study_ID.size(); ++idx_d) {
    const Rcpp::String d = study_ID[idx_d];

    const CharacterVector cov_int_nm = Cov_int_info[d];
    const LogicalMatrix   kp_sample_id = Sample_info[d];

    NumericMatrix est_mat(K, cov_int_nm.size());
    NumericMatrix std_err(K, cov_int_nm.size());
    fill_na(est_mat);
    fill_na(std_err);

    Rcpp::colnames(est_mat) = cov_int_nm;
    Rcpp::rownames(est_mat) = feature_ID;
    Rcpp::colnames(std_err) = cov_int_nm;
    Rcpp::rownames(std_err) = feature_ID;

    const List null_obj_d = null_obj[d];
    const List cov_int_d  = covariate_interest[d];

    arma::mat Y_R = as<arma::mat>(null_obj_d["Y_R"]); // n x m
    arma::mat Y_I = as<arma::mat>(null_obj_d["Y_I"]); // n x m
    arma::mat X    = as<arma::mat>(null_obj_d["Z"]);  // n x p
    arma::vec SUB  = as<arma::vec>(SUB_id[d]);        // n x 1 (1..G integers)

    if (Y_R.n_rows != X.n_rows || Y_I.n_rows != X.n_rows)
      stop("Row mismatch among Y_R, Y_I, and X in study %s.", std::string(d).c_str());

    if (Y_R.n_rows != static_cast<arma::uword>(kp_sample_id.nrow()))
      stop("Sample_info row count does not match data rows in study %s.", std::string(d).c_str());

    if (kp_sample_id.ncol() != cov_int_nm.size())
      stop("Sample_info columns do not match number of covariates in study %s.", std::string(d).c_str());

    for (int cov_name_idx = 0; cov_name_idx < cov_int_nm.size(); ++cov_name_idx) {
      const std::string cov_name = Rcpp::as<std::string>(cov_int_nm[cov_name_idx]);

      const LogicalVector keep_mask = kp_sample_id(Rcpp::_, cov_name_idx);
      const arma::uvec rows = mask_to_uvec(keep_mask);
      const arma::uword n_sub = rows.n_elem;
      if (n_sub == 0) continue;

      arma::mat X_sub  = X.rows(rows);            // n_sub x p
      arma::mat YR_sub = Y_R.rows(rows);          // n_sub x m
      arma::mat YI_sub = Y_I.rows(rows);          // n_sub x m
      arma::vec SUB_sub = SUB.elem(rows);         // n_sub x 1

      NumericVector cov_col_R = as<NumericVector>(cov_int_d[cov_name]);
      if (cov_col_R.size() != static_cast<int>(n_sub))
        stop("Length of covariate '%s' does not match subset size in study %s.",
             cov_name.c_str(), std::string(d).c_str());
      X_sub.col(X_sub.n_cols - 1) = as<arma::vec>(cov_col_R);

      arma::uvec gidx(n_sub);
      arma::uword G = 0;
      for (arma::uword i = 0; i < n_sub; ++i) {
        const double v = SUB_sub[i];
        if (v <= 0 || std::floor(v) != v) stop("SUB_id must be positive integers.");
        const arma::uword gi = static_cast<arma::uword>(v) - 1; // 0-based
        gidx[i] = gi;
        if (gi + 1 > G) G = gi + 1;
      }
      if (G == 0) G = 1;

      const arma::uword m = YR_sub.n_cols;
      if (YI_sub.n_cols != m) stop("Y_R and Y_I must have same number of columns.");

      // Get feature names from the (sub)set—same column order as originals
      Rcpp::CharacterVector feat_names = Rcpp::colnames(Rcpp::as<Rcpp::NumericMatrix>(null_obj_d["Y_I"]));
      if (feat_names.size() == 0)
        stop("Y_I must have column names for feature mapping.");

      for (arma::uword k = 0; k < m; ++k) {
        double est = Rcpp::NumericVector::get_na();
        double var = Rcpp::NumericVector::get_na();

        arma::vec yR_k = YR_sub.col(k);
        arma::vec yI_k = YI_sub.col(k);

        estimate_one_feature(X_sub, yR_k, yI_k, G, gidx, est, var);

        const std::string fname = Rcpp::as<std::string>(feat_names[k]);
        auto it = feat_map.find(fname);
        if (it == feat_map.end())
          stop("Feature '%s' not found in feature_ID.", fname.c_str());
        const int row = it->second; // 0-based

        est_mat(row, cov_name_idx) = est;
        std_err(row, cov_name_idx) = std::sqrt(var);
      }
    }

    List summary_stat_study_one = List::create(
      Named("est")    = est_mat,
      Named("stderr") = std_err
    );
    summary_stat_study[d] = summary_stat_study_one;
  }

  return summary_stat_study;
}
