#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
arma::mat matMult(const arma::mat& A, const arma::mat& B) {
    return A * B;
}

// [[Rcpp::export]]
arma::mat matMultThree(const arma::mat& A, const arma::mat& B, const arma::mat& C) {
    if (A.n_cols != B.n_rows || B.n_cols != C.n_rows) {
        Rcpp::stop("Matrix dimensions are not compatible for multiplication.");
    }
    arma::mat AB = A * B;
    arma::mat ABC = AB * C;
    return ABC;
}

// [[Rcpp::export]]
arma::mat choleskyDecomp(const arma::mat& X) {
    arma::mat L = arma::chol(X, "lower");
    return L;
}

// [[Rcpp::export]]
arma::vec solveArma(const arma::mat& A, const arma::vec& B) {
    arma::vec X = arma::solve(A, B);
    return X;
}

// [[Rcpp::export]]
arma::mat solveArma2(const arma::mat& A, const arma::mat& B) {
    arma::mat X = arma::solve(A, B);
    return X;
}

// [[Rcpp::export]]
arma::mat matInverse(const arma::mat& A) {
    arma::mat A_inv = arma::inv(A);
    return A_inv;
}

// [[Rcpp::export]]
arma::mat compute_L(const arma::mat& X, const arma::mat& Sig_x, const arma::mat& Sig_eps, double nugget) {
    int n = Sig_eps.n_rows;
    arma::mat result = X * Sig_x * X.t();
    result += Sig_eps;
    result += nugget * arma::eye(n, n);
    arma::mat chol_result = arma::chol(result);
    return chol_result.t();
}

// [[Rcpp::export]]
arma::vec compute_w(const arma::mat& L, const arma::vec& y, const arma::mat& X, const arma::vec& mu_x) {
    arma::vec X_mu_x = X * mu_x;
    arma::vec diff = y - X_mu_x;
    arma::vec result = arma::solve(L, diff);
    return result;
}

// [[Rcpp::export]]
arma::mat compute_W(const arma::mat& L, const arma::mat& X, const arma::mat& Sig_x) {
    arma::mat XSig = X * Sig_x;
    arma::mat W = arma::solve(arma::trimatl(L), XSig, arma::solve_opts::fast);
    return W;
}

// [[Rcpp::export]]
arma::vec compute_x_u(const arma::vec& mu_x, const arma::mat& Sig_x, int n) {
    arma::vec z = arma::randn(n);
    arma::mat chol_Sig_x = arma::chol(Sig_x);
    return mu_x + chol_Sig_x.t() * z;
}

// [[Rcpp::export]]
arma::mat compute_x(const arma::mat& x_u, const arma::mat& W, const arma::mat& L, const arma::vec& sss) {
    arma::mat L_inv = arma::inv(L);
    arma::mat Wt_L_inv = W.t() * L_inv;
    arma::mat result = x_u - Wt_L_inv * sss;
    return result;
}

// [[Rcpp::export]]
DataFrame get_MCMC_summary_cpp(const arma::mat& X, const Nullable<NumericVector>& h = R_NilValue) {
    int n = X.n_rows;
    mat summary_dat(n, 3);
    for(int i = 0; i < n; ++i) {
        rowvec row = X.row(i);
        row = sort(row);
        int size = row.n_elem;
        summary_dat(i, 0) = row(std::floor(0.025 * (size - 1)));
        summary_dat(i, 1) = row(std::floor(0.5 * (size - 1)));
        summary_dat(i, 2) = row(std::floor(0.975 * (size - 1)));
    }
    if (h.isNotNull()) {
        NumericVector h_vec(h);
        if (h_vec.length() != n) {
            stop("Length of h must match the number of rows in X");
        }
        return DataFrame::create(
            Named("h") = h_vec,
            Named("lower") = summary_dat.col(0),
            Named("median") = summary_dat.col(1),
            Named("upper") = summary_dat.col(2)
        );
    } else {
        return DataFrame::create(
            Named("lower") = summary_dat.col(0),
            Named("median") = summary_dat.col(1),
            Named("upper") = summary_dat.col(2)
        );
    }
}




