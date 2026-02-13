// plm.cpp
// This file contains implementations of Power Law Model (PLM) functions.
// PLM incorporates a spline component into the rating curve model.
// Functions include density evaluation and prediction for both known and unknown c scenarios.

#include "cppFunctions.h"

using namespace Rcpp;





// [[Rcpp::export]]
Rcpp::List plm_density_evaluation_unknown_c_cpp(const arma::vec& theta,
                                                const arma::mat& P,
                                                const arma::vec& h,
                                                const arma::mat& B,
                                                const arma::vec& y,
                                                const arma::vec& epsilon,
                                                const arma::mat& Sig_x,
                                                const arma::vec& mu_x,
                                                double h_min,
                                                double nugget,
                                                double lambda_c,
                                                double lambda_eta_1,
                                                double lambda_seta) {
    double zeta = theta(0);
    double log_sig_eta = theta(1);
    double eta_1 = theta(2);
    arma::vec z = theta.subvec(3, 7);
    int n = h.n_elem;
    arma::vec lambda = P * arma::join_vert(arma::vec({eta_1}), std::exp(log_sig_eta) * z);
    arma::vec l = arma::log(h - h_min + std::exp(zeta));
    arma::vec varr = epsilon % arma::exp(B * lambda);
    if (arma::any(varr > 100)) {
        return Rcpp::List::create(Rcpp::Named("p") = -1e9);
    }
    arma::mat Sig_eps = arma::diagmat(varr);
    arma::mat X = arma::join_rows(arma::ones(n, 1), l);
    arma::mat M = X * Sig_x * X.t() + Sig_eps;
    M.diag() += nugget;
    arma::mat L = arma::chol(M, "lower");
    arma::vec w = arma::solve(L, y - X * mu_x, arma::solve_opts::fast);
    double p = -0.5 * arma::dot(w, w) - arma::sum(arma::log(L.diag())) +
        pri("c", arma::vec({zeta, lambda_c})) +
        pri("eta_1", arma::vec({eta_1, lambda_eta_1})) +
        pri("eta_minus1", z) +
        pri("sigma_eta", arma::vec({log_sig_eta, lambda_seta}));
    arma::mat W = arma::solve(arma::trimatl(L), X * Sig_x, arma::solve_opts::fast);
    arma::mat chol_Sig_x = arma::chol(Sig_x);
    arma::vec x_u = mu_x + chol_Sig_x.t() * arma::randn(mu_x.n_elem);
    arma::vec sss = X * x_u - y + arma::sqrt(varr) % arma::randn(n);
    arma::vec x = x_u - W.t() * arma::solve(L, sss, arma::solve_opts::fast);
    arma::vec yp = X * x;
    arma::vec ypo = yp + arma::randn(n) % arma::sqrt(varr);
    // Replace NaN with -inf in yp and ypo
    yp.elem(arma::find_nonfinite(yp)).fill(-arma::datum::inf);
    ypo.elem(arma::find_nonfinite(ypo)).fill(-arma::datum::inf);
    // Compute log_lik
    double log_lik = 0.0;
    for(int i = 0; i < n; ++i) {
        log_lik += log_of_normal_pdf(y(i), yp(i), std::sqrt(varr(i)));
    }
    return Rcpp::List::create(
        Rcpp::Named("p") = p,
        Rcpp::Named("x") = x,
        Rcpp::Named("y_post") = yp,
        Rcpp::Named("y_post_pred") = ypo,
        Rcpp::Named("sigma_eps") = varr,
        Rcpp::Named("log_lik") = log_lik
    );
}

// [[Rcpp::export]]
Rcpp::List plm_density_evaluation_known_c_cpp(const arma::vec& theta,
                                              const arma::mat& P,
                                              const arma::vec& h,
                                              const arma::mat& B,
                                              const arma::vec& y,
                                              const arma::vec& epsilon,
                                              const arma::mat& Sig_x,
                                              const arma::vec& mu_x,
                                              double c,
                                              double nugget,
                                              double lambda_eta_1,
                                              double lambda_seta) {
    double log_sig_eta = theta(0);
    double eta_1 = theta(1);
    arma::vec z = theta.subvec(2, 6);
    int n = h.n_elem;
    arma::vec lambda = P * arma::join_vert(arma::vec({eta_1}), std::exp(log_sig_eta) * z);
    arma::vec l = arma::log(h - c);
    arma::vec varr = epsilon % arma::exp(B * lambda);
    if (arma::any(varr > 100)) {
        return Rcpp::List::create(Rcpp::Named("p") = -1e9);
    }
    arma::mat Sig_eps = arma::diagmat(varr);
    arma::mat X = arma::join_rows(arma::ones(n, 1), l);
    arma::mat M = X * Sig_x * X.t() + Sig_eps;
    M.diag() += nugget;
    arma::mat L = arma::chol(M, "lower");
    arma::vec w = arma::solve(L, y - X * mu_x, arma::solve_opts::fast);
    double p = -0.5 * arma::dot(w, w) - arma::sum(arma::log(L.diag())) +
        pri("eta_1", arma::vec({eta_1, lambda_eta_1})) +
        pri("eta_minus1", z) +
        pri("sigma_eta", arma::vec({log_sig_eta, lambda_seta}));
    arma::mat W = arma::solve(arma::trimatl(L), X * Sig_x, arma::solve_opts::fast);
    arma::mat chol_Sig_x = arma::chol(Sig_x);
    arma::vec x_u = mu_x + chol_Sig_x.t() * arma::randn(mu_x.n_elem);
    arma::vec sss = X * x_u - y + arma::sqrt(varr) % arma::randn(n);
    arma::vec x = x_u - W.t() * arma::solve(L, sss, arma::solve_opts::fast);
    arma::vec yp = X * x;
    arma::vec ypo = yp + arma::randn(n) % arma::sqrt(varr);
    // Replace NaN with -inf in yp and ypo
    yp.elem(arma::find_nonfinite(yp)).fill(-arma::datum::inf);
    ypo.elem(arma::find_nonfinite(ypo)).fill(-arma::datum::inf);
    // Compute log_lik
    double log_lik = 0.0;
    for(int i = 0; i < n; ++i) {
        log_lik += log_of_normal_pdf(y(i), yp(i), std::sqrt(varr(i)));
    }
    return Rcpp::List::create(
        Rcpp::Named("p") = p,
        Rcpp::Named("x") = x,
        Rcpp::Named("y_post") = yp,
        Rcpp::Named("y_post_pred") = ypo,
        Rcpp::Named("sigma_eps") = varr,
        Rcpp::Named("log_lik") = log_lik
    );
}

// [[Rcpp::export]]
Rcpp::List plm_predict_u_unknown_c_cpp(const arma::vec& theta,
                                       const arma::vec& x,
                                       const arma::mat& P,
                                       const arma::mat& B_u,
                                       const arma::vec& h_u,
                                       double h_min,
                                       int n_u) {
    double zeta = theta(0);
    double log_sig_eta = theta(1);
    double eta_1 = theta(2);
    arma::vec z = theta.subvec(3, 7);
    int m = n_u;
    arma::vec lambda = P * arma::join_vert(arma::vec({eta_1}), std::exp(log_sig_eta) * z);
    arma::vec varr_u = arma::exp(B_u * lambda);
    arma::uvec above_c = arma::find(h_min - std::exp(zeta) < h_u);
    int m_above_c = above_c.n_elem;
    arma::vec yp_u(m, arma::fill::value(-arma::datum::inf));
    arma::vec ypo_u(m, arma::fill::value(-arma::datum::inf));
    if (m_above_c > 0) {
        arma::vec l = arma::log(h_u.elem(above_c) - h_min + std::exp(zeta));
        arma::mat X = arma::join_rows(arma::ones(m_above_c), l);

        arma::vec yp_u_temp = X * x.subvec(0, 1);
        yp_u.elem(above_c) = yp_u_temp;
        ypo_u.elem(above_c) = yp_u_temp + arma::randn(m_above_c) % arma::sqrt(varr_u.elem(above_c));
    }
    // Replace NaN with -inf in yp_u and ypo_u
    yp_u.elem(arma::find_nonfinite(yp_u)).fill(-arma::datum::inf);
    ypo_u.elem(arma::find_nonfinite(ypo_u)).fill(-arma::datum::inf);
    return Rcpp::List::create(
        Rcpp::Named("y_post") = yp_u,
        Rcpp::Named("y_post_pred") = ypo_u,
        Rcpp::Named("sigma_eps") = varr_u
    );
}

// [[Rcpp::export]]
Rcpp::List plm_predict_u_known_c_cpp(const arma::vec& theta,
                                     const arma::vec& x,
                                     const arma::mat& P,
                                     const arma::mat& B_u,
                                     const arma::vec& h_u,
                                     double c) {
    double log_sig_eta = theta(0);
    double eta_1 = theta(1);
    arma::vec z = theta.subvec(2, 6);
    int m = h_u.n_elem;
    arma::vec lambda = P * arma::join_vert(arma::vec({eta_1}), std::exp(log_sig_eta) * z);
    arma::vec varr_u = arma::exp(B_u * lambda);
    arma::vec l = arma::log(h_u - c);
    arma::mat X = arma::join_rows(arma::ones(m), l);
    arma::vec yp_u = X * x;
    arma::vec ypo_u = yp_u + arma::randn(m) % arma::sqrt(varr_u);
    // Replace NaN with -inf in yp_u and ypo_u
    yp_u.elem(arma::find_nonfinite(yp_u)).fill(-arma::datum::inf);
    ypo_u.elem(arma::find_nonfinite(ypo_u)).fill(-arma::datum::inf);
    return Rcpp::List::create(
        Rcpp::Named("y_post") = yp_u,
        Rcpp::Named("y_post_pred") = ypo_u,
        Rcpp::Named("sigma_eps") = varr_u
    );
}
