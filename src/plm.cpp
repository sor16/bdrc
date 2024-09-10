// plm.cpp
// This file contains implementations of Power Law Model (PLM) functions.
// PLM incorporates a spline component into the rating curve model.
// Functions include density evaluation and prediction for both known and unknown c scenarios.

#include "cppFunctions.h"

using namespace Rcpp;



// [[Rcpp::export]]
Rcpp::List plm_me_density_evaluation_unknown_c_cpp(const arma::vec& theta,
                                                   const arma::mat& P,
                                                   const arma::vec& h,
                                                   const arma::mat& B,
                                                   const arma::vec& y,
                                                   const arma::vec& tau,
                                                   const arma::vec& epsilon,
                                                   const arma::mat& Sig_ab,
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

    arma::mat Sig_eps = arma::diagmat(arma::square(tau));

    arma::mat Sig_u2 = arma::diagmat(varr);

    arma::mat Sig_x = arma::join_cols(
        arma::join_rows(Sig_ab, arma::zeros(2, n)),
        arma::join_rows(arma::zeros(n, 2), Sig_u2)
    );

    arma::mat X = arma::join_rows(arma::ones(n, 1), l, arma::eye(n, n));

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
    arma::vec x_u = mu_x + arma::chol(Sig_x).t() * arma::randn(mu_x.n_elem);
    arma::vec sss = X * x_u - y + tau % arma::randn(n);
    arma::vec x = x_u - W.t() * arma::solve(L, sss, arma::solve_opts::fast);

    arma::vec beta = x.subvec(0, 1);
    arma::vec u2 = x.subvec(2, x.n_elem - 1);
    arma::vec mu_post = beta(0) + beta(1) * l;
    arma::vec y_true = mu_post + u2;
    arma::vec y_true_post_pred = mu_post + arma::randn(n) % arma::sqrt(varr);

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
        Rcpp::Named("x") = beta,
        Rcpp::Named("mu_post") = mu_post,
        Rcpp::Named("y_true") = y_true,
        Rcpp::Named("y_true_post_pred") = y_true_post_pred,
        Rcpp::Named("sigma_eps") = varr,
        Rcpp::Named("log_lik") = log_lik
    );
}

// [[Rcpp::export]]
Rcpp::List plm_density_evaluation_unknown_c_cpp(const arma::vec& theta,
                                                const arma::mat& P,
                                                const arma::vec& h,
                                                const arma::mat& B,
                                                const arma::vec& y,
                                                const arma::vec& epsilon,
                                                const arma::mat& Sig_ab,
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

    arma::mat M = X * Sig_ab * X.t() + Sig_eps;
    M.diag() += nugget;
    arma::mat L = arma::chol(M, "lower");
    arma::vec w = arma::solve(L, y - X * mu_x, arma::solve_opts::fast);

    double p = -0.5 * arma::dot(w, w) - arma::sum(arma::log(L.diag())) +
        pri("c", arma::vec({zeta, lambda_c})) +
        pri("eta_1", arma::vec({eta_1, lambda_eta_1})) +
        pri("eta_minus1", z) +
        pri("sigma_eta", arma::vec({log_sig_eta, lambda_seta}));

    arma::mat W = arma::solve(arma::trimatl(L), X * Sig_ab, arma::solve_opts::fast);
    arma::mat chol_Sig_x = arma::chol(Sig_ab);
    arma::vec x_u = mu_x + chol_Sig_x.t() * arma::randn(mu_x.n_elem);
    arma::vec sss = X * x_u - y + arma::sqrt(varr) % arma::randn(n);
    arma::vec x = x_u - W.t() * arma::solve(L, sss, arma::solve_opts::fast);

    arma::vec mu_post = X * x;
    arma::vec y_true_post_pred = mu_post + arma::randn(n) % arma::sqrt(varr);

    // Replace NaN with -inf in mu_post and y_true_post_pred
    mu_post.elem(arma::find_nonfinite(mu_post)).fill(-arma::datum::inf);
    y_true_post_pred.elem(arma::find_nonfinite(y_true_post_pred)).fill(-arma::datum::inf);

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
        Rcpp::Named("sigma_eps") = varr,
        Rcpp::Named("log_lik") = log_lik
    );
}

// [[Rcpp::export]]
Rcpp::List plm_me_density_evaluation_known_c_cpp(const arma::vec& theta,
                                                 const arma::mat& P,
                                                 const arma::vec& h,
                                                 const arma::mat& B,
                                                 const arma::vec& y,
                                                 const arma::vec& tau,
                                                 const arma::vec& epsilon,
                                                 const arma::mat& Sig_ab,
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

    arma::mat Sig_eps = arma::diagmat(arma::square(tau));

    arma::mat Sig_u2 = arma::diagmat(varr);

    arma::mat Sig_x = arma::join_cols(
        arma::join_rows(Sig_ab, arma::zeros(2, n)),
        arma::join_rows(arma::zeros(n, 2), Sig_u2)
    );

    arma::mat X = arma::join_rows(arma::ones(n, 1), l, arma::eye(n, n));

    arma::mat M = X * Sig_x * X.t() + Sig_eps;
    M.diag() += nugget;
    arma::mat L = arma::chol(M, "lower");
    arma::vec w = arma::solve(L, y - X * mu_x, arma::solve_opts::fast);

    double p = -0.5 * arma::dot(w, w) - arma::sum(arma::log(L.diag())) +
        pri("eta_1", arma::vec({eta_1, lambda_eta_1})) +
        pri("eta_minus1", z) +
        pri("sigma_eta", arma::vec({log_sig_eta, lambda_seta}));

    arma::mat W = arma::solve(arma::trimatl(L), X * Sig_x, arma::solve_opts::fast);
    arma::vec x_u = mu_x + arma::chol(Sig_x).t() * arma::randn(mu_x.n_elem);
    arma::vec sss = X * x_u - y + tau % arma::randn(n);
    arma::vec x = x_u - W.t() * arma::solve(L, sss, arma::solve_opts::fast);

    arma::vec beta = x.subvec(0, 1);
    arma::vec u2 = x.subvec(2, x.n_elem - 1);
    arma::vec mu_post = beta(0) + beta(1) * l;
    arma::vec y_true = mu_post + u2;
    arma::vec y_true_post_pred = mu_post + arma::randn(n) % arma::sqrt(varr);

    // Replace NaN with -inf in yp and ypo
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
        Rcpp::Named("x") = beta,
        Rcpp::Named("mu_post") = mu_post,
        Rcpp::Named("y_true") = y_true,
        Rcpp::Named("y_true_post_pred") = y_true_post_pred,
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
                                              const arma::mat& Sig_ab,
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

    arma::mat M = X * Sig_ab * X.t() + Sig_eps;
    M.diag() += nugget;
    arma::mat L = arma::chol(M, "lower");
    arma::vec w = arma::solve(L, y - X * mu_x, arma::solve_opts::fast);

    double p = -0.5 * arma::dot(w, w) - arma::sum(arma::log(L.diag())) +
        pri("eta_1", arma::vec({eta_1, lambda_eta_1})) +
        pri("eta_minus1", z) +
        pri("sigma_eta", arma::vec({log_sig_eta, lambda_seta}));

    arma::mat W = arma::solve(arma::trimatl(L), X * Sig_ab, arma::solve_opts::fast);
    arma::mat chol_Sig_x = arma::chol(Sig_ab);
    arma::vec x_u = mu_x + chol_Sig_x.t() * arma::randn(mu_x.n_elem);
    arma::vec sss = X * x_u - y + arma::sqrt(varr) % arma::randn(n);
    arma::vec x = x_u - W.t() * arma::solve(L, sss, arma::solve_opts::fast);

    arma::vec mu_post = X * x;
    arma::vec y_true_post_pred = mu_post + arma::randn(n) % arma::sqrt(varr);

    // Replace NaN with -inf in mu_post and y_true_post_pred
    mu_post.elem(arma::find_nonfinite(mu_post)).fill(-arma::datum::inf);
    y_true_post_pred.elem(arma::find_nonfinite(y_true_post_pred)).fill(-arma::datum::inf);

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
    arma::vec mu_post_u(m, arma::fill::value(-arma::datum::inf));
    arma::vec y_true_post_pred_u(m, arma::fill::value(-arma::datum::inf));
    if (m_above_c > 0) {
        arma::vec l = arma::log(h_u.elem(above_c) - h_min + std::exp(zeta));
        arma::mat X = arma::join_rows(arma::ones(m_above_c), l);

        arma::vec mu_post_u_temp = X * x.subvec(0, 1);
        mu_post_u.elem(above_c) = mu_post_u_temp;
        y_true_post_pred_u.elem(above_c) = mu_post_u_temp + arma::randn(m_above_c) % arma::sqrt(varr_u.elem(above_c));
    }
    // Replace NaN with -inf in mu_post_u and y_true_post_pred_u
    mu_post_u.elem(arma::find_nonfinite(mu_post_u)).fill(-arma::datum::inf);
    y_true_post_pred_u.elem(arma::find_nonfinite(y_true_post_pred_u)).fill(-arma::datum::inf);
    return Rcpp::List::create(
        Rcpp::Named("mu_post") = mu_post_u,
        Rcpp::Named("y_true_post_pred") = y_true_post_pred_u,
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
    arma::vec mu_post_u = X * x;
    arma::vec y_true_post_pred_u = mu_post_u + arma::randn(m) % arma::sqrt(varr_u);
    // Replace NaN with -inf in mu_post_u and y_true_post_pred_u
    mu_post_u.elem(arma::find_nonfinite(mu_post_u)).fill(-arma::datum::inf);
    y_true_post_pred_u.elem(arma::find_nonfinite(y_true_post_pred_u)).fill(-arma::datum::inf);
    return Rcpp::List::create(
        Rcpp::Named("mu_post") = mu_post_u,
        Rcpp::Named("y_true_post_pred") = y_true_post_pred_u,
        Rcpp::Named("sigma_eps") = varr_u
    );
}
