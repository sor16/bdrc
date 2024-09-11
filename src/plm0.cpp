// plm0.cpp
// This file contains implementations of Power Law Model (PLM0) functions.
// PLM0 is a simplified version of the rating curve model without the spline component.
// Functions include density evaluation and prediction for both known and unknown c scenarios.

#include "cppFunctions.h"

using namespace Rcpp;




// [[Rcpp::export]]
Rcpp::List plm0_me_density_evaluation_unknown_c_cpp(const arma::vec& theta,
                                                    const arma::vec& h,
                                                    const arma::vec& y,
                                                    const arma::vec& tau,
                                                    const arma::vec& epsilon,
                                                    const arma::mat& Sig_ab,
                                                    const arma::vec& mu_x,
                                                    double h_min,
                                                    double nugget,
                                                    double lambda_c,
                                                    double lambda_se) {
    double zeta = theta(0);
    double log_sig_eps2 = theta(1);
    int n = h.n_elem;

    arma::vec l = arma::log(h - h_min + std::exp(zeta));
    arma::vec varr = epsilon * std::exp(log_sig_eps2);
    if (arma::any(varr > 100)) {
        return Rcpp::List::create(Rcpp::Named("p") = -1e9);
    }

    arma::mat Sig_eps = arma::diagmat(arma::square(tau));

    arma::mat Sig_u2 = arma::diagmat(varr);

    arma::mat Sig_x = arma::join_cols(
        arma::join_rows(Sig_ab, arma::zeros(2, n)),
        arma::join_rows(arma::zeros(n, 2), Sig_u2)
    );

    arma::mat X = arma::join_rows(arma::ones(n), l, arma::eye(n, n));

    arma::mat M = X * Sig_x * X.t() + Sig_eps;
    M.diag() += nugget;
    arma::mat L = arma::chol(M, "lower");
    arma::vec w = arma::solve(L, y - X * mu_x, arma::solve_opts::fast);

    double p = -0.5 * arma::dot(w, w) - arma::sum(arma::log(L.diag())) +
        pri("c", arma::vec({zeta, lambda_c})) +
        pri("sigma_eps2", arma::vec({log_sig_eps2, lambda_se}));

    arma::mat W = arma::solve(arma::trimatl(L), X * Sig_x, arma::solve_opts::fast);
    arma::vec x_u = mu_x + arma::chol(Sig_x).t() * arma::randn(mu_x.n_elem);
    arma::vec sss = X * x_u - y + tau % arma::randn(n);
    arma::vec x = x_u - W.t() * arma::solve(L, sss, arma::solve_opts::fast);

    arma::vec beta = x.subvec(0, 1);
    arma::vec mu_post = X.cols(0, 1) * beta;

    arma::vec y_true = mu_post + x.subvec(2, x.n_elem - 1);
    arma::vec y_true_post_pred = mu_post + arma::randn(n) % arma::sqrt(varr);

    double log_lik = 0.0;
    for(size_t i = 0; i < n; ++i) {
        log_lik += log_of_normal_pdf(y(i), mu_post(i), std::sqrt(varr(i) + std::pow(tau(i), 2)));
    }

    return Rcpp::List::create(
        Rcpp::Named("p") = p,
        Rcpp::Named("x") = beta,
        Rcpp::Named("mu_post") = mu_post,
        Rcpp::Named("y_true") = y_true,
        Rcpp::Named("y_true_post_pred") = y_true_post_pred,
        Rcpp::Named("log_lik") = log_lik
    );
}

// [[Rcpp::export]]
Rcpp::List plm0_density_evaluation_unknown_c_cpp(const arma::vec& theta,
                                                 const arma::vec& h,
                                                 const arma::vec& y,
                                                 const arma::vec& epsilon,
                                                 const arma::mat& Sig_ab,
                                                 const arma::vec& mu_x,
                                                 double h_min,
                                                 double nugget,
                                                 double lambda_c,
                                                 double lambda_se) {
    double zeta = theta(0);
    double log_sig_eps2 = theta(1);
    int n = h.n_elem;

    arma::vec l = arma::log(h - h_min + std::exp(zeta));
    arma::vec varr = epsilon * std::exp(log_sig_eps2);
    if (arma::any(varr > 100)) {
        return Rcpp::List::create(Rcpp::Named("p") = -1e9);
    }

    arma::mat Sig_eps = arma::diagmat(varr);

    arma::mat Sig_x = Sig_ab;

    arma::mat X = arma::join_rows(arma::ones(n), l);

    arma::mat M = X * Sig_x * X.t() + Sig_eps;
    M.diag() += nugget;
    arma::mat L = arma::chol(M, "lower");
    arma::vec w = arma::solve(L, y - X * mu_x, arma::solve_opts::fast);

    double p = -0.5 * arma::dot(w, w) - arma::sum(arma::log(L.diag())) +
        pri("c", arma::vec({zeta, lambda_c})) +
        pri("sigma_eps2", arma::vec({log_sig_eps2, lambda_se}));

    arma::mat W = arma::solve(arma::trimatl(L), X * Sig_x, arma::solve_opts::fast);

    arma::vec x_u = mu_x + arma::chol(Sig_x).t() * arma::randn(mu_x.n_elem);
    arma::vec sss = X * x_u - y + arma::sqrt(varr) % arma::randn(n);
    arma::vec x = x_u - W.t() * arma::solve(L, sss, arma::solve_opts::fast);

    arma::vec mu_post = X * x;
    arma::vec y_true_post_pred = mu_post + arma::randn(n) % arma::sqrt(varr);

    // Compute log_lik
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
Rcpp::List plm0_me_density_evaluation_known_c_cpp(const arma::vec& theta,
                                                  const arma::vec& h,
                                                  const arma::vec& y,
                                                  const arma::vec& tau,
                                                  const arma::vec& epsilon,
                                                  const arma::mat& Sig_ab,
                                                  const arma::vec& mu_x,
                                                  double c,
                                                  double nugget,
                                                  double lambda_se) {
    double log_sig_eps2 = theta(0);
    int n = h.n_elem;

    arma::vec l = arma::log(h - c);
    arma::vec varr = epsilon * std::exp(log_sig_eps2);
    if (arma::any(varr > 100)) {
        return Rcpp::List::create(Rcpp::Named("p") = -1e9);
    }

    arma::mat Sig_eps = arma::diagmat(arma::square(tau));

    arma::mat Sig_u2 = arma::diagmat(varr);

    arma::mat Sig_x = arma::join_cols(
        arma::join_rows(Sig_ab, arma::zeros(2, n)),
        arma::join_rows(arma::zeros(n, 2), Sig_u2)
    );

    arma::mat X = arma::join_rows(arma::ones(n), l, arma::eye(n, n));

    arma::mat M = X * Sig_x * X.t() + Sig_eps;
    M.diag() += nugget;
    arma::mat L = arma::chol(M, "lower");
    arma::vec w = arma::solve(L, y - X * mu_x, arma::solve_opts::fast);

    double p = -0.5 * arma::dot(w, w) - arma::sum(arma::log(L.diag())) +
        pri("sigma_eps2", arma::vec({log_sig_eps2, lambda_se}));

    arma::mat W = arma::solve(arma::trimatl(L), X * Sig_x, arma::solve_opts::fast);
    arma::vec x_u = mu_x + arma::chol(Sig_x).t() * arma::randn(mu_x.n_elem);
    arma::vec sss = X * x_u - y + tau % arma::randn(n);
    arma::vec x = x_u - W.t() * arma::solve(L, sss, arma::solve_opts::fast);

    arma::vec beta = x.subvec(0, 1);
    arma::vec mu_post = X.cols(0, 1) * beta;

    arma::vec y_true = mu_post + x.subvec(2, x.n_elem - 1);
    arma::vec y_true_post_pred = mu_post + arma::randn(n) % arma::sqrt(varr);

    double log_lik = 0.0;
    for(size_t i = 0; i < n; ++i) {
        log_lik += log_of_normal_pdf(y(i), mu_post(i), std::sqrt(varr(i) + std::pow(tau(i), 2)));
    }

    return Rcpp::List::create(
        Rcpp::Named("p") = p,
        Rcpp::Named("x") = beta,
        Rcpp::Named("mu_post") = mu_post,
        Rcpp::Named("y_true") = y_true,
        Rcpp::Named("y_true_post_pred") = y_true_post_pred,
        Rcpp::Named("log_lik") = log_lik
    );
}

// [[Rcpp::export]]
Rcpp::List plm0_density_evaluation_known_c_cpp(const arma::vec& theta,
                                               const arma::vec& h,
                                               const arma::vec& y,
                                               const arma::vec& epsilon,
                                               const arma::mat& Sig_ab,
                                               const arma::vec& mu_x,
                                               double c,
                                               double nugget,
                                               double lambda_se) {
    double log_sig_eps2 = theta(0);
    int n = h.n_elem;

    arma::vec l = arma::log(h - c);
    arma::vec varr = epsilon * std::exp(log_sig_eps2);
    if (arma::any(varr > 100)) {
        return Rcpp::List::create(Rcpp::Named("p") = -1e9);
    }

    arma::mat Sig_eps = arma::diagmat(varr);

    arma::mat Sig_x = Sig_ab;

    arma::mat X = arma::join_rows(arma::ones(n), l);

    arma::mat M = X * Sig_x * X.t() + Sig_eps;
    M.diag() += nugget;
    arma::mat L = arma::chol(M, "lower");
    arma::vec w = arma::solve(L, y - X * mu_x, arma::solve_opts::fast);

    double p = -0.5 * arma::dot(w, w) - arma::sum(arma::log(L.diag())) +
        pri("sigma_eps2", arma::vec({log_sig_eps2, lambda_se}));

    arma::mat W = arma::solve(arma::trimatl(L), X * Sig_x, arma::solve_opts::fast);
    arma::vec x_u = mu_x + arma::chol(Sig_x).t() * arma::randn(mu_x.n_elem);
    arma::vec sss = X * x_u - y + arma::sqrt(varr) % arma::randn(n);
    arma::vec x = x_u - W.t() * arma::solve(L, sss, arma::solve_opts::fast);

    arma::vec mu_post = X * x;
    arma::vec y_true_post_pred = mu_post + arma::randn(n) % arma::sqrt(varr);

    // Compute log_lik
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
Rcpp::List plm0_predict_u_unknown_c_cpp(const arma::vec& theta,
                                        const arma::vec& x,
                                        const arma::vec& h_u,
                                        double h_min) {
    double zeta = theta(0);
    double log_sig_eps2 = theta(1);
    int m = h_u.n_elem;

    arma::uvec above_c = arma::find(h_min - std::exp(zeta) < h_u);
    int m_above_c = above_c.n_elem;
    arma::vec mu_post_u(m, arma::fill::value(-arma::datum::inf));
    arma::vec y_true_post_pred_u(m, arma::fill::value(-arma::datum::inf));
    if (m_above_c > 0) {
        arma::vec l = arma::log(h_u.elem(above_c) - h_min + std::exp(zeta));
        arma::mat X = arma::join_rows(arma::ones(m_above_c), l);
        arma::vec mu_post_u_temp = X * x;
        mu_post_u.elem(above_c) = mu_post_u_temp;
        y_true_post_pred_u.elem(above_c) = mu_post_u_temp + arma::randn(m_above_c) * std::sqrt(std::exp(log_sig_eps2));
    }
    return Rcpp::List::create(
        Rcpp::Named("mu_post") = mu_post_u,
        Rcpp::Named("y_true_post_pred") = y_true_post_pred_u
    );
}

// [[Rcpp::export]]
Rcpp::List plm0_predict_u_known_c_cpp(const arma::vec& theta,
                                      const arma::vec& x,
                                      const arma::vec& h_u,
                                      double c) {
    double log_sig_eps2 = theta(0);
    int m = h_u.n_elem;
    arma::vec l = arma::log(h_u - c);
    arma::mat X = arma::join_rows(arma::ones(m), l);
    arma::vec mu_post_u = X * x;
    arma::vec y_true_post_pred_u = mu_post_u + arma::randn(m) * std::sqrt(std::exp(log_sig_eps2));
    return Rcpp::List::create(
        Rcpp::Named("mu_post") = mu_post_u,
        Rcpp::Named("y_true_post_pred") = y_true_post_pred_u
    );
}


