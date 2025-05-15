#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
arma::cube setSlice(arma::cube M, const NumericMatrix& M1, int i) {
  // Retrieve dimensions of M
  int dim1 = M.n_rows;
  int dim2 = M.n_cols;
  int dim3 = M.n_slices;

  // Validate slice index
  if(i < 1 || i > dim1){
    stop("Index i is out of bounds. It must be between 1 and %d.", dim1);
  }

  // Validate dimensions of M1
  if(M1.nrow() != dim2 || M1.ncol() != dim3){
    stop("Dimensions of M1 must match (dim2 x dim3). Provided M1 has dimensions (%d x %d), but expected (%d x %d).",
         M1.nrow(), M1.ncol(), dim2, dim3);
  }

  // Convert NumericMatrix M1 to arma::mat
  arma::mat M1_arma = as<arma::mat>(M1);


  // Adjust for 0-based indexing in C++
  int i_idx = i - 1;

  // Assign M1_arma to the i-th slice of M
  // Perform element-wise assignment to the i-th slice
  for(int j = 0; j < dim2; j++){
    for(int k = 0; k < dim3; k++){
      M(i_idx, j, k) = M1_arma(j, k);
    }
  }

  return M;
}

NumericMatrix invert_matrix(NumericMatrix mat) {
  // Convert R matrix to Armadillo matrix
  arma::mat arma_mat = as<arma::mat>(mat);

  // Invert the matrix using Armadillo
  arma::mat arma_inv = arma::inv(arma_mat);

  // Convert back to R matrix
  NumericMatrix inv_mat = wrap(arma_inv);

  return inv_mat;
}

NumericMatrix subset_matrix(const NumericMatrix& mat,
                            const IntegerVector& row_indices,
                            const IntegerVector& col_indices) {
  int nrows = row_indices.size();
  int ncols = col_indices.size();

  NumericMatrix result(nrows, ncols);

  // Copy values
  for (int i = 0; i < nrows; ++i) {
    for (int j = 0; j < ncols; ++j) {
      result(i, j) = mat(row_indices[i], col_indices[j]);
    }
  }

  // Set row names if they exist
  CharacterVector original_rownames = rownames(mat);
  if (original_rownames.size() > 0) {
    CharacterVector new_rownames(nrows);
    for (int i = 0; i < nrows; ++i) {
      new_rownames[i] = original_rownames[row_indices[i]];
    }
    rownames(result) = new_rownames;
  }

  // Set column names if they exist
  CharacterVector original_colnames = colnames(mat);
  if (original_colnames.size() > 0) {
    CharacterVector new_colnames(ncols);
    for (int j = 0; j < ncols; ++j) {
      new_colnames[j] = original_colnames[col_indices[j]];
    }
    colnames(result) = new_colnames;
  }

  return result;
}

NumericMatrix na_matrix(int n, int k){
  NumericMatrix m(n,k) ;
  std::fill( m.begin(), m.end(), NumericVector::get_na() ) ;
  return m ;
}

NumericMatrix create_matrix(NumericVector X) {
  int n = X.size();
  NumericMatrix mat(n, n);

  for(int l = 0; l < n; ++l) {
    for(int ll = 0; ll < n; ++ll) {
      mat(l, ll) = X[l] * X[ll];
    }
  }

  return mat;
}

int which_character(StringVector X, String a) {
  for(int i = 0; i < X.size(); ++i) {
    if(X[i] == a) {
      return i; // R is 1-indexed
    }
  }
  return NA_INTEGER; // Return NA if no match is found
}

NumericMatrix matrixMultiply(NumericMatrix M1, NumericMatrix M2) {
  int n1_rows = M1.nrow();
  int n1_cols = M1.ncol();
  int n2_rows = M2.nrow();
  int n2_cols = M2.ncol();

  // Check if the number of columns in M1 equals the number of rows in M2
  if (n1_cols != n2_rows) {
    stop("Number of columns in M1 must equal number of rows in M2 for multiplication.");
  }

  // Initialize the output matrix with dimensions n1_rows x n2_cols
  NumericMatrix M_output(n1_rows, n2_cols);

  // Perform matrix multiplication
  for(int i = 0; i < n1_rows; ++i) {
    for(int j = 0; j < n2_cols; ++j) {
      double sum = 0.0;
      for(int k = 0; k < n1_cols; ++k) {
        sum += M1(i, k) * M2(k, j);
      }
      M_output(i, j) = sum;
    }
  }

  return M_output;
}

NumericMatrix matrixSubtract(NumericMatrix M1, NumericMatrix M2) {
  int n1_rows = M1.nrow();
  int n1_cols = M1.ncol();
  int n2_rows = M2.nrow();
  int n2_cols = M2.ncol();

  // Check if dimensions match
  if (n1_rows != n2_rows || n1_cols != n2_cols) {
    stop("Matrices M1 and M2 must have the same dimensions for subtraction.");
  }

  // Initialize the output matrix with the same dimensions
  NumericMatrix M_sub(n1_rows, n1_cols);

  // Perform matrix subtraction
  for(int i = 0; i < n1_rows; ++i) {
    for(int j = 0; j < n1_cols; ++j) {
      M_sub(i, j) = M1(i, j) - M2(i, j);
    }
  }

  return M_sub;
}

NumericMatrix matrixAdd(NumericMatrix M1, NumericMatrix M2) {
  int n1_rows = M1.nrow();
  int n1_cols = M1.ncol();
  int n2_rows = M2.nrow();
  int n2_cols = M2.ncol();

  // Check if dimensions match
  if (n1_rows != n2_rows || n1_cols != n2_cols) {
    stop("Matrices M1 and M2 must have the same dimensions for adding.");
  }

  // Initialize the output matrix with the same dimensions
  NumericMatrix M_sub(n1_rows, n1_cols);

  // Perform matrix subtraction
  for(int i = 0; i < n1_rows; ++i) {
    for(int j = 0; j < n1_cols; ++j) {
      M_sub(i, j) = M1(i, j) + M2(i, j);
    }
  }

  return M_sub;
}

NumericMatrix columnSumsMatrixCpp(NumericMatrix M) {
  int n_rows = M.nrow();
  int n_cols = M.ncol();

  // Initialize the output matrix with n_cols rows and 1 column
  NumericMatrix M_sub(n_cols, 1);

  for(int j = 0; j < n_cols; ++j) {          // Iterate over columns
    double sum = 0.0;
    for(int i = 0; i < n_rows; ++i) {      // Iterate over rows
      sum += M(i, j);
    }
    M_sub(j, 0) = sum;
  }

  // Optionally, assign a column name
  M_sub.attr("dimnames") = List::create(
    R_NilValue,                // No row names
    CharacterVector::create("Sum") // Column name
  );

  return M_sub;
}

NumericMatrix transposeMatrix(NumericMatrix M) {
  int n_rows = M.nrow();
  int n_cols = M.ncol();

  // Initialize the transposed matrix with dimensions n_cols x n_rows
  NumericMatrix M_transposed(n_cols, n_rows);

  // Perform the transposition
  for(int i = 0; i < n_rows; ++i){
    for(int j = 0; j < n_cols; ++j){
      M_transposed(j, i) = M(i, j);
    }
  }

  return M_transposed;
}

// [[Rcpp::export]]
List palm_rcpp(
    List null_obj,
    List covariate_interest, // assuming covariate_interest is also passed as input
    List SUB_id,
    CharacterVector study_ID,
    CharacterVector feature_ID,
    List Cov_int_info,
    List Sample_info
) {
  int K = feature_ID.size();

  List summary_stat_study;

  for (int idx_d = 0; idx_d < study_ID.size(); ++idx_d) {
    String d = study_ID[idx_d];
    CharacterVector cov_int_nm = Cov_int_info[d];
    LogicalMatrix kp_sample_id = Sample_info[d];
    NumericMatrix est_mat = na_matrix(K, cov_int_nm.size());
    NumericMatrix std_err = na_matrix(K, cov_int_nm.size());

    colnames(est_mat) = cov_int_nm;
    rownames(est_mat) = feature_ID;
    colnames(std_err) = cov_int_nm;
    rownames(std_err) = feature_ID;

    List null_obj_d = null_obj[d];
    List cov_int_d = covariate_interest[d];
    NumericMatrix Y_R = as<NumericMatrix>(null_obj_d["Y_R"]);
    NumericMatrix Y_I = as<NumericMatrix>(null_obj_d["Y_I"]);
    NumericMatrix X = as<NumericMatrix>(null_obj_d["Z"]);
    NumericVector SUBid = as<NumericVector>(SUB_id[d]);

    for (int cov_name_idx = 0; cov_name_idx < cov_int_nm.size(); ++cov_name_idx) {
      String cov_name = cov_int_nm[cov_name_idx];

      LogicalVector kp_indices = kp_sample_id(_, cov_name_idx);
      IntegerVector kp_indices_int = seq_along(kp_indices) - 1; // indices for subsetting
      IntegerVector kp_indices_int_sub = kp_indices_int[kp_indices];

      NumericMatrix Y_R_sub = subset_matrix(Y_R, kp_indices_int_sub, Range(0,Y_R.ncol()-1));
      NumericMatrix Y_I_sub = subset_matrix(Y_I, kp_indices_int_sub, Range(0,Y_I.ncol()-1));
      NumericMatrix X_sub = subset_matrix(X,kp_indices_int_sub, Range(0,X.ncol()-1));
      NumericVector SUBid_sub = SUBid[kp_indices_int_sub];

      NumericVector cov_col = as<NumericVector>(cov_int_d[cov_name]);

      for (int i = 0; i < kp_indices_int_sub.size(); ++i) {
        X_sub(i, X_sub.ncol() - 1) = cov_col[i];
      }
      NumericVector ests, covs;

      List XX_lst(SUBid_sub.size());
      for (int l = 0; l < SUBid_sub.size(); ++l) {
        NumericMatrix XX = create_matrix(X_sub(l, _));
        XX_lst[l] = XX;
      }

      for (int k = 0; k < Y_I_sub.ncol(); ++k) {
        NumericMatrix s_i_lst(max(SUBid_sub), X_sub.ncol());
        NumericMatrix I_mat(X_sub.ncol(), X_sub.ncol());

        for (int l = 0; l < SUBid_sub.size(); ++l) {
          NumericVector temp_row = s_i_lst(SUBid_sub(l) - 1, _); // Extract the row
          temp_row += Y_R_sub(l, k) * X_sub(l, _);               // Perform element-wise addition
          s_i_lst(SUBid_sub(l) - 1, _) = temp_row;

          NumericMatrix XX = XX_lst[l];
          I_mat += Y_I_sub(l, k) * XX;
        }

        int beta_name = X_sub.ncol() - 1;
        IntegerVector gamma_name = seq(0, X_sub.ncol() - 2);

        NumericMatrix I_gamma = subset_matrix(I_mat, gamma_name, gamma_name);
        NumericMatrix I_beta(1, 1);
        I_beta(0, 0) = I_mat(beta_name, beta_name);
        NumericMatrix Inv_I_gamma = invert_matrix(I_gamma);

        for (int l1 = 0; l1 < gamma_name.size(); ++l1) {
          for (int l2 = 0; l2 < gamma_name.size(); ++l2) {
            I_beta(0, 0) -= I_mat(beta_name, l1) * Inv_I_gamma(l1, l2) * I_mat(beta_name, l2);
          }
        }

        double est = sum(s_i_lst(_, beta_name)) / I_beta(0, 0);
        ests.push_back(est);

        double core_U = 0;
        for (int i = 0; i < s_i_lst.nrow(); ++i) {
          double tmp_U = s_i_lst(i, beta_name);
          for (int l1 = 0; l1 < gamma_name.size(); ++l1) {
            for (int l2 = 0; l2 < gamma_name.size(); ++l2) {
              tmp_U -= I_mat(beta_name, l1) * Inv_I_gamma(l1, l2) * s_i_lst(i, l2);
            }
          }
          core_U += tmp_U * tmp_U;
        }

        double cov = core_U / (I_beta(0, 0) * I_beta(0, 0));
        covs.push_back(cov);
      }

      CharacterVector feature_idd = colnames(Y_I_sub);
      for (int j = 0; j < Y_I_sub.ncol(); ++j) {
        String a = feature_idd[j];
        int tmp_id = which_character(feature_ID, a);
        est_mat(tmp_id, cov_name_idx) = ests[j];
        std_err(tmp_id, cov_name_idx) = sqrt(covs[j]);
      }
    }

    List summary_stat_study_one = List::create(
      Named("est") = est_mat,
      Named("stderr") = std_err
    );

    summary_stat_study[d] = summary_stat_study_one;
  }

  return summary_stat_study;
}
