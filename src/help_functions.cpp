// help_functions.cpp
// This file contains helper functions used across different model implementations.
// These functions support various mathematical operations, probability calculations,
// and utility tasks required by the PLM, PLM0, GPLM, and GPLM0 models.

#include "cppFunctions.h"

using namespace Rcpp;



// [[Rcpp::export]]
arma::mat matMult(const arma::mat& A, const arma::mat& B) {
    return A * B;
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
DataFrame get_MCMC_summary_cpp(const arma::mat& X, const Nullable<NumericVector>& h ) {
    int n = X.n_rows;
    arma::mat summary_dat(n, 3);
    for(int i = 0; i < n; ++i) {
        arma::rowvec row = X.row(i);
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

// [[Rcpp::export]]
arma::vec calc_variogram_arma(const arma::mat& param_mat, int i, int burnin, int nr_iter ) {
    int n_cols = param_mat.n_cols;
    return arma::sum(arma::square(
            param_mat.cols(i, n_cols - 1) -
                param_mat.cols(0, n_cols - i - 1)
    ), 1);
}

// [[Rcpp::export]]
arma::mat variogram_chain(int T_max, const arma::mat& param_mat1, const arma::mat& param_mat2, int burnin, int nr_iter ) {
    int n_rows = param_mat1.n_rows;
    arma::mat result(n_rows, 2 * T_max);
    for (int i = 0; i < T_max; ++i) {
        result.col(i) = calc_variogram_arma(param_mat1, i + 1, burnin, nr_iter);
        result.col(i + T_max) = calc_variogram_arma(param_mat2, i + 1, burnin, nr_iter);
    }
    return result;
}

// [[Rcpp::export]]
arma::mat distance_matrix(const arma::vec& x) {
    int n = x.n_elem;
    arma::mat dist_mat(n, n, arma::fill::zeros);
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            double dist = std::abs(x(i) - x(j));
            dist_mat(i, j) = dist;
            dist_mat(j, i) = dist; // Since distance is symmetric
        }
    }
    return dist_mat;
}

// [[Rcpp::export]]
arma::mat create_A_cpp(const arma::vec& h) {
    int n = h.n_elem;
    arma::vec unique_h = arma::unique(h);
    int unique_n = unique_h.n_elem;
    arma::mat A = arma::zeros<arma::mat>(n, unique_n);
    A(0, 0) = 1;
    int i = 0;
    for (int ii = 1; ii < n; ++ii) {
        if (h(ii) == h(ii - 1)) {
            A(ii, i) = 1;
        } else {
            i += 1;
            A(ii, i) = 1;
        }
    }
    return A;
}

// [[Rcpp::export]]
double pri(const std::string& type, const arma::vec& args) {
    double result = 0.0;
    if (type == "c") {
        result = args(0) - std::exp(args(0)) * args(1);
    } else if (type == "sigma_eps2") {
        result = 0.5 * args(0) - std::exp(0.5 * args(0)) * args(1);
    } else if (type == "sigma_b") {
        result = args(0) - std::exp(args(0)) * args(1);
    } else if (type == "phi_b") {
        result = -0.5 * args(0) - args(1) * std::sqrt(0.5) * std::exp(-0.5 * args(0));
    } else if (type == "eta_1") {
        result = 0.5 * args(0) - std::exp(0.5 * args(0)) * args(1);
    } else if (type == "eta_minus1") {
        result = -0.5 * arma::dot(args, args);
    } else if (type == "sigma_eta") {
        result = args(0) - std::exp(args(0)) * args(1);
    }

    if (!std::isfinite(result)) {
        // Rcpp::Rcout << "Warning: Non-finite pri value. Type: " << type
        //             << ", Args: " << args.t() << std::endl;
        return -std::numeric_limits<double>::max();
    }

    return result;
}

double log_of_normal_pdf(double x, double mu, double sigma) {
    static const double log_2pi = std::log(2.0 * M_PI);
    return -0.5 * (log_2pi + 2.0 * std::log(sigma) + std::pow((x - mu) / sigma, 2));
}

// [[Rcpp::export]]
Rcpp::List chain_statistics_cpp(const arma::mat& chains) {
    int chains_length = chains.n_rows;
    int split_idx = chains_length / 2;

    // Adjust split_idx for odd lengths
    bool is_odd = chains_length % 2 != 0;
    if (is_odd) {
        split_idx += 1;
    }

    arma::mat chains1 = chains.rows(0, split_idx - 1);
    arma::mat chains2 = chains.rows(split_idx - (is_odd ? 1 : 0), chains_length - 1);

    arma::rowvec means1 = arma::mean(chains1, 0);
    arma::rowvec means2 = arma::mean(chains2, 0);
    arma::rowvec vars1 = arma::var(chains1, 0);
    arma::rowvec vars2 = arma::var(chains2, 0);

    double within_chain_var = (arma::mean(vars1) + arma::mean(vars2)) / 2;

    arma::rowvec mean_diff = means1 - means2;
    double between_chain_var = split_idx * arma::var(mean_diff);

    // Adjust n for var_hat calculation
    int n = is_odd ? split_idx - 1 : split_idx;
    double var_hat = ((n - 1) * within_chain_var + between_chain_var) / n;

    return Rcpp::List::create(
        Rcpp::Named("W") = within_chain_var,
        Rcpp::Named("var_hat") = var_hat
    );
}



