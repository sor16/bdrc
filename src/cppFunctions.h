#ifndef CPP_FUNCTIONS_H
#define CPP_FUNCTIONS_H

#define _USE_MATH_DEFINES

#include <RcppArmadillo.h>
#include <cmath>

// Help function declarations
arma::mat matMult(const arma::mat& A, const arma::mat& B);
arma::mat choleskyDecomp(const arma::mat& X);
arma::vec solveArma(const arma::mat& A, const arma::vec& B);
arma::mat matInverse(const arma::mat& A);
arma::mat compute_L(const arma::mat& X, const arma::mat& Sig_x, const arma::mat& Sig_eps, double nugget);
arma::vec compute_w(const arma::mat& L, const arma::vec& y, const arma::mat& X, const arma::vec& mu_x);
Rcpp::DataFrame get_MCMC_summary_cpp(const arma::mat& X, const Rcpp::Nullable<Rcpp::NumericVector>& h = R_NilValue);
arma::vec calc_variogram_arma(const arma::mat& param_mat, int i, int burnin = 2000, int nr_iter = 20000);
arma::mat variogram_chain(int T_max, const arma::mat& param_mat1, const arma::mat& param_mat2, int burnin = 2000, int nr_iter = 20000);
arma::mat distance_matrix(const arma::vec& x);
arma::mat create_A_cpp(const arma::vec& h);
double pri(const std::string& type, const arma::vec& args);
double log_of_normal_pdf(double x, double mu, double sigma);
Rcpp::List chain_statistics_cpp(const arma::mat& chains);

// GPLM function declarations
Rcpp::List gplm_density_evaluation_unknown_c_cpp(const arma::vec& theta, const arma::mat& P, const arma::vec& h, const arma::mat& B, const arma::mat& dist, const arma::mat& A, const arma::vec& y, const arma::vec& epsilon, double h_min, double nugget, int n_unique, const arma::vec& mu_x, const arma::mat& Sig_ab, const arma::mat& Z, double lambda_c, double lambda_sb, double lambda_pb, double lambda_eta_1, double lambda_seta);
Rcpp::List gplm_density_evaluation_known_c_cpp(const arma::vec& theta, const arma::mat& P, const arma::vec& h, const arma::mat& B, const arma::mat& dist, const arma::mat& A, const arma::vec& y, const arma::vec& epsilon, double nugget, int n_unique, const arma::vec& mu_x, const arma::mat& Sig_ab, const arma::mat& Z, double lambda_sb, double lambda_pb, double lambda_eta_1, double lambda_seta, double c);
Rcpp::List gplm_predict_u_unknown_c_cpp(const arma::vec& theta, const arma::vec& x, const arma::mat& P, const arma::mat& B_u, const arma::vec& h_unique, const arma::vec& h_u, const arma::mat& dist_all, double h_min, double nugget, int n_unique, int n_u);
Rcpp::List gplm_predict_u_known_c_cpp(const arma::vec& theta, const arma::vec& x, const arma::mat& P, const arma::mat& B_u, const arma::vec& h_unique, const arma::vec& h_u, const arma::mat& dist_all, double c, double nugget, int n_unique, int n_u);

// GPLM0 function declarations
Rcpp::List gplm0_density_evaluation_unknown_c_cpp(const arma::vec& theta, const arma::vec& h, const arma::vec& y, const arma::mat& A, const arma::mat& dist, const arma::vec& epsilon, double h_min, double nugget, int n_unique, const arma::vec& mu_x, const arma::mat& Sig_ab, const arma::mat& Z, double lambda_c, double lambda_se, double lambda_sb, double lambda_pb);
Rcpp::List gplm0_density_evaluation_known_c_cpp(const arma::vec& theta, const arma::vec& h, const arma::vec& y, const arma::mat& A, const arma::mat& dist, const arma::vec& epsilon, double c, double nugget, int n_unique, const arma::vec& mu_x, const arma::mat& Sig_ab, const arma::mat& Z, double lambda_se, double lambda_sb, double lambda_pb);
Rcpp::List gplm0_predict_u_unknown_c_cpp(const arma::vec& theta, const arma::vec& x, const arma::vec& h_unique, const arma::vec& h_u, const arma::mat& dist_all, double h_min, double nugget, int n_unique, int n_u);
Rcpp::List gplm0_predict_u_known_c_cpp(const arma::vec& theta, const arma::vec& x, const arma::vec& h_unique, const arma::vec& h_u, const arma::mat& dist_all, double c, double nugget, int n_unique, int n_u);

// PLM function declarations
Rcpp::List plm_density_evaluation_unknown_c_cpp(const arma::vec& theta, const arma::mat& P, const arma::vec& h, const arma::mat& B, const arma::vec& y, const arma::vec& epsilon, const arma::mat& Sig_x, const arma::vec& mu_x, double h_min, double nugget, double lambda_c, double lambda_eta_1, double lambda_seta);
Rcpp::List plm_density_evaluation_known_c_cpp(const arma::vec& theta, const arma::mat& P, const arma::vec& h, const arma::mat& B, const arma::vec& y, const arma::vec& epsilon, const arma::mat& Sig_x, const arma::vec& mu_x, double c, double nugget, double lambda_eta_1, double lambda_seta);
Rcpp::List plm_predict_u_unknown_c_cpp(const arma::vec& theta, const arma::vec& x, const arma::mat& P, const arma::mat& B_u, const arma::vec& h_u, double h_min, int n_u);
Rcpp::List plm_predict_u_known_c_cpp(const arma::vec& theta, const arma::vec& x, const arma::mat& P, const arma::mat& B_u, const arma::vec& h_u, double c);

// PLM0 function declarations
Rcpp::List plm0_density_evaluation_unknown_c_cpp(const arma::vec& theta, const arma::vec& h, const arma::vec& y, const arma::vec& epsilon, const arma::mat& Sig_ab, const arma::vec& mu_x, double h_min, double nugget, double lambda_c, double lambda_se);
Rcpp::List plm0_density_evaluation_known_c_cpp(const arma::vec& theta, const arma::vec& h, const arma::vec& y, const arma::vec& epsilon, const arma::mat& Sig_ab, const arma::vec& mu_x, double c, double nugget, double lambda_se);
Rcpp::List plm0_predict_u_unknown_c_cpp(const arma::vec& theta, const arma::vec& x, const arma::vec& h_u, double h_min);
Rcpp::List plm0_predict_u_known_c_cpp(const arma::vec& theta, const arma::vec& x, const arma::vec& h_u, double c);

#endif // CPP_FUNCTIONS_H
