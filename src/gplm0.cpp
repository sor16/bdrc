// gplm0.cpp
// This file contains implementations of Generalized Power Law Model (GPLM0) functions.
// GPLM0 extends the PLM0 by incorporating spatial correlation in the beta parameters.
// Functions include density evaluation and prediction for both known and unknown c scenarios.

#include "cppFunctions.h"

using namespace Rcpp;



// [[Rcpp::export]]
Rcpp::List gplm0_me_density_evaluation_unknown_c_cpp(const arma::vec& theta,
                                                     const arma::vec& h,
                                                     const arma::vec& y,
                                                     const arma::mat& A,
                                                     const arma::mat& dist,
                                                     const arma::vec& epsilon,
                                                     const arma::vec& tau,
                                                     double h_min,
                                                     double nugget,
                                                     int n_unique,
                                                     const arma::vec& mu_x,
                                                     const arma::mat& Sig_ab,
                                                     const arma::mat& Z,
                                                     double lambda_c,
                                                     double lambda_se,
                                                     double lambda_sb,
                                                     double lambda_pb) {
    double zeta = theta(0);
    double log_sig_eps2 = theta(1);
    double log_sig_b = theta(2);
    double log_phi_b = theta(3);
    int n = h.n_elem;

    arma::vec l = arma::log(h - h_min + std::exp(zeta));
    arma::vec varr = epsilon * std::exp(log_sig_eps2);
    if (arma::any(varr > 100)) {
        return Rcpp::List::create(Rcpp::Named("p") = -1e9);
    }

    // Construct Sigma_eps using tau (transformed measurement errors) and add zero at the end for mean-zero GP
    arma::mat Sig_eps = arma::diagmat(arma::join_vert(arma::square(tau), arma::vec({0.0})));

    // Construct Sigma_u2 (for y_true) without adding zero at the end
    arma::mat Sig_u2 = arma::diagmat(varr);

    // Matern covariance
    arma::mat R_Beta = (1.0 + std::sqrt(5.0) * dist / std::exp(log_phi_b) +
        5.0 * arma::square(dist) / (3.0 * std::pow(std::exp(log_phi_b), 2))) %
        arma::exp(-std::sqrt(5.0) * dist / std::exp(log_phi_b));
    R_Beta.diag() += nugget;

    // Construct Sigma_x
    arma::mat Sig_x = arma::join_cols(
        arma::join_rows(Sig_ab, arma::zeros(2, n_unique), arma::zeros(2, n)),
        arma::join_rows(arma::zeros(n_unique, 2), std::exp(2 * log_sig_b) * R_Beta, arma::zeros(n_unique, n)),
        arma::join_rows(arma::zeros(n, 2 + n_unique), Sig_u2)
    );

    // Construct X matrix
    arma::mat X = arma::join_cols(
        arma::join_rows(arma::ones(n, 1), l, arma::diagmat(l) * A, arma::eye(n, n)),
        Z
    );

    // Compute L using Cholesky decomposition
    arma::mat M = X * Sig_x * X.t() + Sig_eps;
    M.diag() += nugget;
    arma::mat L = arma::chol(M, "lower");
    arma::vec w = arma::solve(L, y - X * mu_x, arma::solve_opts::fast);

    double p = -0.5 * arma::dot(w, w) - arma::sum(arma::log(L.diag())) +
        pri("c", arma::vec({zeta, lambda_c})) +
        pri("sigma_eps2", arma::vec({log_sig_eps2, lambda_se})) +
        pri("sigma_b", arma::vec({log_sig_b, lambda_sb})) +
        pri("phi_b", arma::vec({log_phi_b, lambda_pb}));

    // Compute posterior samples
    arma::mat W = arma::solve(arma::trimatl(L), X * Sig_x, arma::solve_opts::fast);
    arma::vec x_u = mu_x + arma::chol(Sig_x).t() * arma::randn(mu_x.n_elem);
    arma::vec sss = X * x_u - y + arma::join_vert(tau % arma::randn(n), arma::vec({0.0}));
    arma::vec x = x_u - W.t() * arma::solve(L, sss, arma::solve_opts::fast);

    // Extract subset of x (excluding delta values)
    arma::vec x_subset = x.subvec(0, 1 + n_unique);

    // Extract components from x
    arma::vec beta = x.subvec(0, 1);  // a_0 and b
    arma::vec u1 = x.subvec(2, 1 + n_unique);  // beta(h)
    arma::vec u2 = x.subvec(2 + n_unique, x.n_elem - 1);  // delta

    // Compute mu_post
    arma::vec mu_post = beta(0) + (beta(1) + A * u1) % l;

    // Compute y_true_post_pred (which includes only model error, not measurement error)
    arma::vec y_true_post_pred = mu_post + arma::randn(n) % arma::sqrt(varr);

    // Compute y_true (which includes model error)
    arma::vec y_true = mu_post + u2;

    // Replace NaN with -inf in mu_post, y_true and y_true_post_pred
    mu_post.elem(arma::find_nonfinite(mu_post)).fill(-arma::datum::inf);
    y_true.elem(arma::find_nonfinite(y_true)).fill(-arma::datum::inf);
    y_true_post_pred.elem(arma::find_nonfinite(y_true_post_pred)).fill(-arma::datum::inf);

    // Compute log_lik
    double log_lik = 0.0;
    for(size_t i = 0; i < n; ++i) {
        log_lik += log_of_normal_pdf(y(i), mu_post(i), std::sqrt(varr(i) + std::pow(tau(i), 2)));
    }

    return Rcpp::List::create(
        Rcpp::Named("p") = p,
        Rcpp::Named("x") = x_subset,
        Rcpp::Named("mu_post") = mu_post,
        Rcpp::Named("y_true") = y_true,
        Rcpp::Named("y_true_post_pred") = y_true_post_pred,
        Rcpp::Named("log_lik") = log_lik
    );
}

// [[Rcpp::export]]
Rcpp::List gplm0_density_evaluation_unknown_c_cpp(const arma::vec& theta,
                                                  const arma::vec& h,
                                                  const arma::vec& y,
                                                  const arma::mat& A,
                                                  const arma::mat& dist,
                                                  const arma::vec& epsilon,
                                                  double h_min,
                                                  double nugget,
                                                  int n_unique,
                                                  const arma::vec& mu_x,
                                                  const arma::mat& Sig_ab,
                                                  const arma::mat& Z,
                                                  double lambda_c,
                                                  double lambda_se,
                                                  double lambda_sb,
                                                  double lambda_pb) {
    double zeta = theta(0);
    double log_sig_eps2 = theta(1);
    double log_sig_b = theta(2);
    double log_phi_b = theta(3);
    int n = h.n_elem;

    arma::vec l = arma::log(h - h_min + std::exp(zeta));
    arma::vec varr = epsilon * std::exp(log_sig_eps2);
    if (arma::any(varr > 100)) {
        return Rcpp::List::create(Rcpp::Named("p") = -1e9);
    }

    arma::mat Sig_eps = arma::diagmat(arma::join_vert(varr, arma::vec({0.0})));

    arma::mat R_Beta = (1.0 + std::sqrt(5.0) * dist / std::exp(log_phi_b) +
        5.0 * arma::square(dist) / (3.0 * std::pow(std::exp(log_phi_b), 2))) %
        arma::exp(-std::sqrt(5.0) * dist / std::exp(log_phi_b));
    R_Beta.diag() += nugget;

    arma::mat m1 = arma::zeros(2, n_unique);
    arma::mat m2 = arma::zeros(n_unique, 2);

    arma::mat Sig_x = arma::join_cols(
        arma::join_rows(Sig_ab, m1),
        arma::join_rows(m2, std::exp(2 * log_sig_b) * R_Beta)
    );

    arma::mat X = arma::join_cols(
        arma::join_rows(arma::ones(n, 1), l, arma::diagmat(l) * A),
        Z
    );

    arma::mat M = X * Sig_x * X.t() + Sig_eps;
    M.diag() += nugget;
    arma::mat L = arma::chol(M, "lower");
    arma::vec w = arma::solve(L, y - X * mu_x, arma::solve_opts::fast);

    double p = -0.5 * arma::dot(w, w) - arma::sum(arma::log(L.diag())) +
        pri("c", arma::vec({zeta, lambda_c})) +
        pri("sigma_eps2", arma::vec({log_sig_eps2, lambda_se})) +
        pri("sigma_b", arma::vec({log_sig_b, lambda_sb})) +
        pri("phi_b", arma::vec({log_phi_b, lambda_pb}));

    arma::mat W = arma::solve(arma::trimatl(L), X * Sig_x, arma::solve_opts::fast);
    arma::vec x_u = mu_x + arma::chol(Sig_x).t() * arma::randn(n_unique + 2);
    arma::vec sss = X * x_u - y + arma::join_vert(arma::sqrt(varr) % arma::randn(n), arma::vec({0.0}));
    arma::mat x = x_u - W.t() * arma::solve(L, sss, arma::solve_opts::fast);

    arma::vec mu_post = arma::vec(X * x).subvec(0, n-1);
    arma::vec y_true_post_pred = mu_post + arma::randn(n) % arma::sqrt(varr);

    mu_post.elem(arma::find_nonfinite(mu_post)).fill(-arma::datum::inf);
    y_true_post_pred.elem(arma::find_nonfinite(y_true_post_pred)).fill(-arma::datum::inf);

    double log_lik = 0.0;
    for(size_t i = 0; i < n; ++i) {
        log_lik += log_of_normal_pdf(y(i), mu_post(i), std::sqrt(varr(i)));
    }

    return Rcpp::List::create(
        Rcpp::Named("p") = p,
        Rcpp::Named("x") = x,
        Rcpp::Named("mu_post") = mu_post,
        Rcpp::Named("y_true_post_pred") = y_true_post_pred,
        Rcpp::Named("log_lik") = log_lik
    );
}

// [[Rcpp::export]]
Rcpp::List gplm0_me_density_evaluation_known_c_cpp(const arma::vec& theta,
                                                   const arma::vec& h,
                                                   const arma::vec& y,
                                                   const arma::mat& A,
                                                   const arma::mat& dist,
                                                   const arma::vec& epsilon,
                                                   const arma::vec& tau,
                                                   double c,
                                                   double nugget,
                                                   int n_unique,
                                                   const arma::vec& mu_x,
                                                   const arma::mat& Sig_ab,
                                                   const arma::mat& Z,
                                                   double lambda_se,
                                                   double lambda_sb,
                                                   double lambda_pb) {
    double log_sig_eps2 = theta(0);
    double log_sig_b = theta(1);
    double log_phi_b = theta(2);
    int n = h.n_elem;

    arma::vec l = arma::log(h - c);
    arma::vec varr = epsilon * std::exp(log_sig_eps2);
    if (arma::any(varr > 100)) {
        return Rcpp::List::create(Rcpp::Named("p") = -1e9);
    }

    // Construct Sigma_eps using tau instead of varr
    arma::mat Sig_eps = arma::diagmat(arma::join_vert(arma::square(tau), arma::vec({0.0})));

    // Construct Sigma_u2 (for y_true) without adding zero at the end
    arma::mat Sig_u2 = arma::diagmat(varr);

    // Matern covariance
    arma::mat R_Beta = (1.0 + std::sqrt(5.0) * dist / std::exp(log_phi_b) +
        5.0 * arma::square(dist) / (3.0 * std::pow(std::exp(log_phi_b), 2))) %
        arma::exp(-std::sqrt(5.0) * dist / std::exp(log_phi_b));
    R_Beta.diag() += nugget;

    // Construct Sigma_x to include Sig_u2
    arma::mat Sig_x = arma::join_cols(
        arma::join_rows(Sig_ab, arma::zeros(2, n_unique), arma::zeros(2, n)),
        arma::join_rows(arma::zeros(n_unique, 2), std::exp(2 * log_sig_b) * R_Beta, arma::zeros(n_unique, n)),
        arma::join_rows(arma::zeros(n, 2 + n_unique), Sig_u2)
    );

    // Construct X matrix to include identity matrix for delta terms
    arma::mat X = arma::join_cols(
        arma::join_rows(arma::ones(n, 1), l, arma::diagmat(l) * A, arma::eye(n, n)),
        Z
    );

    // Compute L using Cholesky decomposition
    arma::mat M = X * Sig_x * X.t() + Sig_eps;
    M.diag() += nugget;
    arma::mat L = arma::chol(M, "lower");
    arma::vec w = arma::solve(L, y - X * mu_x, arma::solve_opts::fast);

    double p = -0.5 * arma::dot(w, w) - arma::sum(arma::log(L.diag())) +
        pri("sigma_eps2", arma::vec({log_sig_eps2, lambda_se})) +
        pri("sigma_b", arma::vec({log_sig_b, lambda_sb})) +
        pri("phi_b", arma::vec({log_phi_b, lambda_pb}));


    // Compute posterior samples
    arma::mat W = arma::solve(arma::trimatl(L), X * Sig_x, arma::solve_opts::fast);
    arma::vec x_u = mu_x + arma::chol(Sig_x).t() * arma::randn(mu_x.n_elem);
    arma::vec sss = X * x_u - y + arma::join_vert(tau % arma::randn(n), arma::vec({0.0}));
    arma::vec x = x_u - W.t() * arma::solve(L, sss, arma::solve_opts::fast);

    // Extract subset of x (excluding delta values)
    arma::vec x_subset = x.subvec(0, 1 + n_unique);

    // Extract components from x
    arma::vec beta = x.subvec(0, 1);  // a_0 and b
    arma::vec u1 = x.subvec(2, 1 + n_unique);  // beta(h)
    arma::vec u2 = x.subvec(2 + n_unique, x.n_elem - 1);  // delta

    // Compute mu_post instead of yp
    arma::vec mu_post = beta(0) + (beta(1) + A * u1) % l;

    // Compute y_true (which includes model error)
    arma::vec y_true = mu_post + u2;

    // Compute y_true_post_pred instead of ypo
    arma::vec y_true_post_pred = mu_post + arma::randn(n) % arma::sqrt(varr);

    // Replace NaN with -inf in mu_post, y_true and y_true_post_pred
    mu_post.elem(arma::find_nonfinite(mu_post)).fill(-arma::datum::inf);
    y_true.elem(arma::find_nonfinite(y_true)).fill(-arma::datum::inf);
    y_true_post_pred.elem(arma::find_nonfinite(y_true_post_pred)).fill(-arma::datum::inf);

    // Compute log_lik
    double log_lik = 0.0;
    for(size_t i = 0; i < n; ++i) {
        log_lik += log_of_normal_pdf(y(i), mu_post(i), std::sqrt(varr(i) + std::pow(tau(i), 2)));
    }

    return Rcpp::List::create(
        Rcpp::Named("p") = p,
        Rcpp::Named("x") = x_subset,
        Rcpp::Named("mu_post") = mu_post,
        Rcpp::Named("y_true") = y_true,
        Rcpp::Named("y_true_post_pred") = y_true_post_pred,
        Rcpp::Named("log_lik") = log_lik
    );
}

// [[Rcpp::export]]
Rcpp::List gplm0_density_evaluation_known_c_cpp(const arma::vec& theta,
                                                const arma::vec& h,
                                                const arma::vec& y,
                                                const arma::mat& A,
                                                const arma::mat& dist,
                                                const arma::vec& epsilon,
                                                double c,
                                                double nugget,
                                                int n_unique,
                                                const arma::vec& mu_x,
                                                const arma::mat& Sig_ab,
                                                const arma::mat& Z,
                                                double lambda_se,
                                                double lambda_sb,
                                                double lambda_pb) {
    double log_sig_eps2 = theta(0);
    double log_sig_b = theta(1);
    double log_phi_b = theta(2);
    int n = h.n_elem;

    arma::vec l = arma::log(h - c);
    arma::vec varr = epsilon * std::exp(log_sig_eps2);
    if (arma::any(varr > 100)) {
        return Rcpp::List::create(Rcpp::Named("p") = -1e9);
    }

    arma::mat Sig_eps = arma::diagmat(arma::join_vert(varr, arma::vec({0.0})));

    // Matern covariance
    arma::mat R_Beta = (1.0 + std::sqrt(5.0) * dist / std::exp(log_phi_b) +
        5.0 * arma::square(dist) / (3.0 * std::pow(std::exp(log_phi_b), 2))) %
        arma::exp(-std::sqrt(5.0) * dist / std::exp(log_phi_b));
    R_Beta.diag() += nugget;

    arma::mat m1 = arma::zeros(2, n_unique);
    arma::mat m2 = arma::zeros(n_unique, 2);
    arma::mat Sig_x = arma::join_cols(
        arma::join_rows(Sig_ab, m1),
        arma::join_rows(m2, std::exp(2 * log_sig_b) * R_Beta)
    );

    arma::mat X = arma::join_cols(
        arma::join_rows(arma::ones(n, 1), l, arma::diagmat(l) * A),
        Z
    );

    arma::mat M = X * Sig_x * X.t() + Sig_eps;
    M.diag() += nugget;
    arma::mat L = arma::chol(M, "lower");
    arma::vec w = arma::solve(L, y - X * mu_x, arma::solve_opts::fast);

    double p = -0.5 * arma::dot(w, w) - arma::sum(arma::log(L.diag())) +
        pri("sigma_eps2", arma::vec({log_sig_eps2, lambda_se})) +
        pri("sigma_b", arma::vec({log_sig_b, lambda_sb})) +
        pri("phi_b", arma::vec({log_phi_b, lambda_pb}));

    arma::mat W = arma::solve(arma::trimatl(L), X * Sig_x, arma::solve_opts::fast);
    arma::vec x_u = mu_x + arma::chol(Sig_x).t() * arma::randn(n_unique + 2);
    arma::vec sss = X * x_u - y + arma::join_vert(arma::sqrt(varr) % arma::randn(n), arma::vec({0.0}));
    arma::mat x = x_u - W.t() * arma::solve(L, sss, arma::solve_opts::fast);

    arma::vec mu_post = arma::vec(X * x).subvec(0, n-1);
    arma::vec y_true_post_pred = mu_post + arma::randn(n) % arma::sqrt(varr);

    mu_post.elem(arma::find_nonfinite(mu_post)).fill(-arma::datum::inf);
    y_true_post_pred.elem(arma::find_nonfinite(y_true_post_pred)).fill(-arma::datum::inf);

    double log_lik = 0.0;
    for(size_t i = 0; i < n; ++i) {
        log_lik += log_of_normal_pdf(y(i), mu_post(i), std::sqrt(varr(i)));
    }

    return Rcpp::List::create(
        Rcpp::Named("p") = p,
        Rcpp::Named("x") = x,
        Rcpp::Named("mu_post") = mu_post,
        Rcpp::Named("y_true_post_pred") = y_true_post_pred,
        Rcpp::Named("log_lik") = log_lik
    );
}

// [[Rcpp::export]]
Rcpp::List gplm0_predict_u_unknown_c_cpp(const arma::vec& theta,
                                         const arma::vec& x,
                                         const arma::vec& h_unique,
                                         const arma::vec& h_u,
                                         const arma::mat& dist_all,
                                         double h_min,
                                         double nugget,
                                         int n_unique,
                                         int n_u) {
    int n = n_unique;
    int m = n_u;

    arma::vec beta_u(m);
    arma::vec mu_post_u(m, arma::fill::value(-arma::datum::inf));
    arma::vec y_true_post_pred_u(m, arma::fill::value(-arma::datum::inf));
    arma::vec varr_u(m);

    double zeta = theta(0);
    double log_sig_eps2 = theta(1);
    double sig_b = std::exp(theta(2));
    double phi_b = std::exp(theta(3));

    // Repeated calculations
    double sqrt_5 = std::sqrt(5.0);

    // Get sample of data variance
    varr_u.fill(std::exp(log_sig_eps2));

    // Covariance of the joint prior for betas from data and beta unobserved
    arma::mat sigma_all = std::pow(sig_b, 2) *
        (1.0 + sqrt_5 * dist_all / phi_b + 5.0 * arma::square(dist_all) / (3.0 * std::pow(phi_b, 2))) %
        arma::exp(-sqrt_5 * dist_all / phi_b);

    sigma_all.diag() += nugget;
    arma::mat sigma_11 = sigma_all.submat(0, 0, n-1, n-1);
    arma::mat sigma_22 = sigma_all.submat(n, n, n+m-1, n+m-1);
    arma::mat sigma_12 = sigma_all.submat(0, n, n-1, n+m-1);
    arma::mat sigma_21 = sigma_all.submat(n, 0, n+m-1, n-1);

    // Parameters for the posterior of beta_u
    arma::vec x_subset = x.subvec(2, x.n_elem - 1);
    arma::vec mu_x_u = sigma_21 * arma::solve(sigma_11, x_subset);
    arma::mat Sigma_x_u = sigma_22 - sigma_21 * arma::solve(sigma_11, sigma_12);

    // A sample from posterior of beta_u drawn
    beta_u = mu_x_u + arma::chol(Sigma_x_u, "lower") * arma::randn(m);
    arma::uvec above_c = arma::find(-(std::exp(zeta) - h_min) < h_u);
    int m_above_c = above_c.n_elem;
    if (m_above_c > 0) {
        arma::vec l = arma::log(h_u.elem(above_c) - h_min + std::exp(zeta));
        arma::mat X;
        if (l.n_elem > 1) {
            X = arma::join_rows(arma::ones(m_above_c), l, arma::diagmat(l));
        } else {
            X = arma::mat({1, l(0), l(0)}).t();
        }
        arma::vec x_u(m_above_c + 2);
        x_u(0) = x(0);
        x_u(1) = x(1);
        for (int j = 0; j < m_above_c; ++j) {
            x_u(j+2) = beta_u(above_c(j));
        }
        arma::vec mu_post_u_temp = X * x_u;
        for (int j = 0; j < m_above_c; ++j) {
            mu_post_u(above_c(j)) = mu_post_u_temp(j);
            y_true_post_pred_u(above_c(j)) = mu_post_u_temp(j) + arma::randn() * std::sqrt(varr_u(above_c(j)));
        }
    }

    // Replace NaN with -inf in mu_post_u and y_true_post_pred_u
    mu_post_u.elem(arma::find_nonfinite(mu_post_u)).fill(-arma::datum::inf);
    y_true_post_pred_u.elem(arma::find_nonfinite(y_true_post_pred_u)).fill(-arma::datum::inf);

    return Rcpp::List::create(
        Rcpp::Named("x") = arma::trans(beta_u),
        Rcpp::Named("mu_post") = mu_post_u,
        Rcpp::Named("y_true_post_pred") = y_true_post_pred_u
    );
}

// [[Rcpp::export]]
Rcpp::List gplm0_predict_u_known_c_cpp(const arma::vec& theta,
                                       const arma::vec& x,
                                       const arma::vec& h_unique,
                                       const arma::vec& h_u,
                                       const arma::mat& dist_all,
                                       double c,
                                       double nugget,
                                       int n_unique,
                                       int n_u) {
    int n = n_unique;
    int m = n_u;
    double log_sig_eps2 = theta(0);
    double sig_b = std::exp(theta(1));
    double phi_b = std::exp(theta(2));

    // Repeated calculations
    double sqrt_5 = std::sqrt(5.0);

    // Get sample of data variance
    arma::vec varr_u = arma::vec(m, arma::fill::value(std::exp(log_sig_eps2)));

    // Covariance of the joint prior for betas from data and beta unobserved
    arma::mat sigma_all = std::pow(sig_b, 2) *
        (1.0 + sqrt_5 * dist_all / phi_b + 5.0 * arma::square(dist_all) / (3.0 * std::pow(phi_b, 2))) %
        arma::exp(-sqrt_5 * dist_all / phi_b);

    sigma_all.diag() += nugget;
    arma::mat sigma_11 = sigma_all.submat(0, 0, n-1, n-1);
    arma::mat sigma_22 = sigma_all.submat(n, n, n+m-1, n+m-1);
    arma::mat sigma_12 = sigma_all.submat(0, n, n-1, n+m-1);
    arma::mat sigma_21 = sigma_all.submat(n, 0, n+m-1, n-1);

    // Parameters for the posterior of beta_u
    arma::vec x_subset = x.subvec(2, x.n_elem - 1);
    arma::vec mu_x_u = sigma_21 * arma::solve(sigma_11, x_subset);
    arma::mat Sigma_x_u = sigma_22 - sigma_21 * arma::solve(sigma_11, sigma_12);

    // A sample from posterior of beta_u drawn
    arma::vec beta_u = mu_x_u + arma::chol(Sigma_x_u, "lower") * arma::randn(m);

    // Building blocks of the explanatory matrix X calculated
    arma::vec l = arma::log(h_u - c);
    arma::mat X;
    if (l.n_elem > 1) {
        X = arma::join_rows(arma::ones(m), l, arma::diagmat(l));
    } else {
        X = arma::mat({1, l(0), l(0)}).t();
    }
    arma::vec x_u(m + 2);
    x_u(0) = x(0);
    x_u(1) = x(1);
    x_u.subvec(2, m+1) = beta_u;

    // Sample from the posterior of discharge y
    arma::vec mu_post_u = X * x_u;

    // Make sure the log discharge at point of zero discharge is -Inf
    mu_post_u(0) = -arma::datum::inf;
    arma::vec y_true_post_pred_u = mu_post_u + arma::randn(m) % arma::sqrt(varr_u);

    // Replace NaN with -inf in mu_post_u and y_true_post_pred_u
    mu_post_u.elem(arma::find_nonfinite(mu_post_u)).fill(-arma::datum::inf);
    y_true_post_pred_u.elem(arma::find_nonfinite(y_true_post_pred_u)).fill(-arma::datum::inf);

    return Rcpp::List::create(
        Rcpp::Named("x") = arma::trans(beta_u),
        Rcpp::Named("mu_post") = mu_post_u,
        Rcpp::Named("y_true_post_pred") = y_true_post_pred_u
    );
}

