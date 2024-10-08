// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// gplm_density_evaluation_unknown_c_cpp
Rcpp::List gplm_density_evaluation_unknown_c_cpp(const arma::vec& theta, const arma::mat& P, const arma::vec& h, const arma::mat& B, const arma::mat& dist, const arma::mat& A, const arma::vec& y, const arma::vec& epsilon, double h_min, double nugget, int n_unique, const arma::vec& mu_x, const arma::mat& Sig_ab, const arma::mat& Z, double lambda_c, double lambda_sb, double lambda_pb, double lambda_eta_1, double lambda_seta);
RcppExport SEXP _bdrc_gplm_density_evaluation_unknown_c_cpp(SEXP thetaSEXP, SEXP PSEXP, SEXP hSEXP, SEXP BSEXP, SEXP distSEXP, SEXP ASEXP, SEXP ySEXP, SEXP epsilonSEXP, SEXP h_minSEXP, SEXP nuggetSEXP, SEXP n_uniqueSEXP, SEXP mu_xSEXP, SEXP Sig_abSEXP, SEXP ZSEXP, SEXP lambda_cSEXP, SEXP lambda_sbSEXP, SEXP lambda_pbSEXP, SEXP lambda_eta_1SEXP, SEXP lambda_setaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type P(PSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type h(hSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type B(BSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type dist(distSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< double >::type h_min(h_minSEXP);
    Rcpp::traits::input_parameter< double >::type nugget(nuggetSEXP);
    Rcpp::traits::input_parameter< int >::type n_unique(n_uniqueSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu_x(mu_xSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Sig_ab(Sig_abSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_c(lambda_cSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_sb(lambda_sbSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_pb(lambda_pbSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_eta_1(lambda_eta_1SEXP);
    Rcpp::traits::input_parameter< double >::type lambda_seta(lambda_setaSEXP);
    rcpp_result_gen = Rcpp::wrap(gplm_density_evaluation_unknown_c_cpp(theta, P, h, B, dist, A, y, epsilon, h_min, nugget, n_unique, mu_x, Sig_ab, Z, lambda_c, lambda_sb, lambda_pb, lambda_eta_1, lambda_seta));
    return rcpp_result_gen;
END_RCPP
}
// gplm_density_evaluation_known_c_cpp
Rcpp::List gplm_density_evaluation_known_c_cpp(const arma::vec& theta, const arma::mat& P, const arma::vec& h, const arma::mat& B, const arma::mat& dist, const arma::mat& A, const arma::vec& y, const arma::vec& epsilon, double nugget, int n_unique, const arma::vec& mu_x, const arma::mat& Sig_ab, const arma::mat& Z, double lambda_sb, double lambda_pb, double lambda_eta_1, double lambda_seta, double c);
RcppExport SEXP _bdrc_gplm_density_evaluation_known_c_cpp(SEXP thetaSEXP, SEXP PSEXP, SEXP hSEXP, SEXP BSEXP, SEXP distSEXP, SEXP ASEXP, SEXP ySEXP, SEXP epsilonSEXP, SEXP nuggetSEXP, SEXP n_uniqueSEXP, SEXP mu_xSEXP, SEXP Sig_abSEXP, SEXP ZSEXP, SEXP lambda_sbSEXP, SEXP lambda_pbSEXP, SEXP lambda_eta_1SEXP, SEXP lambda_setaSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type P(PSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type h(hSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type B(BSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type dist(distSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< double >::type nugget(nuggetSEXP);
    Rcpp::traits::input_parameter< int >::type n_unique(n_uniqueSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu_x(mu_xSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Sig_ab(Sig_abSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_sb(lambda_sbSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_pb(lambda_pbSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_eta_1(lambda_eta_1SEXP);
    Rcpp::traits::input_parameter< double >::type lambda_seta(lambda_setaSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(gplm_density_evaluation_known_c_cpp(theta, P, h, B, dist, A, y, epsilon, nugget, n_unique, mu_x, Sig_ab, Z, lambda_sb, lambda_pb, lambda_eta_1, lambda_seta, c));
    return rcpp_result_gen;
END_RCPP
}
// gplm_predict_u_unknown_c_cpp
Rcpp::List gplm_predict_u_unknown_c_cpp(const arma::vec& theta, const arma::vec& x, const arma::mat& P, const arma::mat& B_u, const arma::vec& h_unique, const arma::vec& h_u, const arma::mat& dist_all, double h_min, double nugget, int n_unique, int n_u);
RcppExport SEXP _bdrc_gplm_predict_u_unknown_c_cpp(SEXP thetaSEXP, SEXP xSEXP, SEXP PSEXP, SEXP B_uSEXP, SEXP h_uniqueSEXP, SEXP h_uSEXP, SEXP dist_allSEXP, SEXP h_minSEXP, SEXP nuggetSEXP, SEXP n_uniqueSEXP, SEXP n_uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type P(PSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type B_u(B_uSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type h_unique(h_uniqueSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type h_u(h_uSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type dist_all(dist_allSEXP);
    Rcpp::traits::input_parameter< double >::type h_min(h_minSEXP);
    Rcpp::traits::input_parameter< double >::type nugget(nuggetSEXP);
    Rcpp::traits::input_parameter< int >::type n_unique(n_uniqueSEXP);
    Rcpp::traits::input_parameter< int >::type n_u(n_uSEXP);
    rcpp_result_gen = Rcpp::wrap(gplm_predict_u_unknown_c_cpp(theta, x, P, B_u, h_unique, h_u, dist_all, h_min, nugget, n_unique, n_u));
    return rcpp_result_gen;
END_RCPP
}
// gplm_predict_u_known_c_cpp
Rcpp::List gplm_predict_u_known_c_cpp(const arma::vec& theta, const arma::vec& x, const arma::mat& P, const arma::mat& B_u, const arma::vec& h_unique, const arma::vec& h_u, const arma::mat& dist_all, double c, double nugget, int n_unique, int n_u);
RcppExport SEXP _bdrc_gplm_predict_u_known_c_cpp(SEXP thetaSEXP, SEXP xSEXP, SEXP PSEXP, SEXP B_uSEXP, SEXP h_uniqueSEXP, SEXP h_uSEXP, SEXP dist_allSEXP, SEXP cSEXP, SEXP nuggetSEXP, SEXP n_uniqueSEXP, SEXP n_uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type P(PSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type B_u(B_uSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type h_unique(h_uniqueSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type h_u(h_uSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type dist_all(dist_allSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    Rcpp::traits::input_parameter< double >::type nugget(nuggetSEXP);
    Rcpp::traits::input_parameter< int >::type n_unique(n_uniqueSEXP);
    Rcpp::traits::input_parameter< int >::type n_u(n_uSEXP);
    rcpp_result_gen = Rcpp::wrap(gplm_predict_u_known_c_cpp(theta, x, P, B_u, h_unique, h_u, dist_all, c, nugget, n_unique, n_u));
    return rcpp_result_gen;
END_RCPP
}
// gplm0_density_evaluation_unknown_c_cpp
Rcpp::List gplm0_density_evaluation_unknown_c_cpp(const arma::vec& theta, const arma::vec& h, const arma::vec& y, const arma::mat& A, const arma::mat& dist, const arma::vec& epsilon, double h_min, double nugget, int n_unique, const arma::vec& mu_x, const arma::mat& Sig_ab, const arma::mat& Z, double lambda_c, double lambda_se, double lambda_sb, double lambda_pb);
RcppExport SEXP _bdrc_gplm0_density_evaluation_unknown_c_cpp(SEXP thetaSEXP, SEXP hSEXP, SEXP ySEXP, SEXP ASEXP, SEXP distSEXP, SEXP epsilonSEXP, SEXP h_minSEXP, SEXP nuggetSEXP, SEXP n_uniqueSEXP, SEXP mu_xSEXP, SEXP Sig_abSEXP, SEXP ZSEXP, SEXP lambda_cSEXP, SEXP lambda_seSEXP, SEXP lambda_sbSEXP, SEXP lambda_pbSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type h(hSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type dist(distSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< double >::type h_min(h_minSEXP);
    Rcpp::traits::input_parameter< double >::type nugget(nuggetSEXP);
    Rcpp::traits::input_parameter< int >::type n_unique(n_uniqueSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu_x(mu_xSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Sig_ab(Sig_abSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_c(lambda_cSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_se(lambda_seSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_sb(lambda_sbSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_pb(lambda_pbSEXP);
    rcpp_result_gen = Rcpp::wrap(gplm0_density_evaluation_unknown_c_cpp(theta, h, y, A, dist, epsilon, h_min, nugget, n_unique, mu_x, Sig_ab, Z, lambda_c, lambda_se, lambda_sb, lambda_pb));
    return rcpp_result_gen;
END_RCPP
}
// gplm0_density_evaluation_known_c_cpp
Rcpp::List gplm0_density_evaluation_known_c_cpp(const arma::vec& theta, const arma::vec& h, const arma::vec& y, const arma::mat& A, const arma::mat& dist, const arma::vec& epsilon, double c, double nugget, int n_unique, const arma::vec& mu_x, const arma::mat& Sig_ab, const arma::mat& Z, double lambda_se, double lambda_sb, double lambda_pb);
RcppExport SEXP _bdrc_gplm0_density_evaluation_known_c_cpp(SEXP thetaSEXP, SEXP hSEXP, SEXP ySEXP, SEXP ASEXP, SEXP distSEXP, SEXP epsilonSEXP, SEXP cSEXP, SEXP nuggetSEXP, SEXP n_uniqueSEXP, SEXP mu_xSEXP, SEXP Sig_abSEXP, SEXP ZSEXP, SEXP lambda_seSEXP, SEXP lambda_sbSEXP, SEXP lambda_pbSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type h(hSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type dist(distSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    Rcpp::traits::input_parameter< double >::type nugget(nuggetSEXP);
    Rcpp::traits::input_parameter< int >::type n_unique(n_uniqueSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu_x(mu_xSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Sig_ab(Sig_abSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_se(lambda_seSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_sb(lambda_sbSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_pb(lambda_pbSEXP);
    rcpp_result_gen = Rcpp::wrap(gplm0_density_evaluation_known_c_cpp(theta, h, y, A, dist, epsilon, c, nugget, n_unique, mu_x, Sig_ab, Z, lambda_se, lambda_sb, lambda_pb));
    return rcpp_result_gen;
END_RCPP
}
// gplm0_predict_u_unknown_c_cpp
Rcpp::List gplm0_predict_u_unknown_c_cpp(const arma::vec& theta, const arma::vec& x, const arma::vec& h_unique, const arma::vec& h_u, const arma::mat& dist_all, double h_min, double nugget, int n_unique, int n_u);
RcppExport SEXP _bdrc_gplm0_predict_u_unknown_c_cpp(SEXP thetaSEXP, SEXP xSEXP, SEXP h_uniqueSEXP, SEXP h_uSEXP, SEXP dist_allSEXP, SEXP h_minSEXP, SEXP nuggetSEXP, SEXP n_uniqueSEXP, SEXP n_uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type h_unique(h_uniqueSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type h_u(h_uSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type dist_all(dist_allSEXP);
    Rcpp::traits::input_parameter< double >::type h_min(h_minSEXP);
    Rcpp::traits::input_parameter< double >::type nugget(nuggetSEXP);
    Rcpp::traits::input_parameter< int >::type n_unique(n_uniqueSEXP);
    Rcpp::traits::input_parameter< int >::type n_u(n_uSEXP);
    rcpp_result_gen = Rcpp::wrap(gplm0_predict_u_unknown_c_cpp(theta, x, h_unique, h_u, dist_all, h_min, nugget, n_unique, n_u));
    return rcpp_result_gen;
END_RCPP
}
// gplm0_predict_u_known_c_cpp
Rcpp::List gplm0_predict_u_known_c_cpp(const arma::vec& theta, const arma::vec& x, const arma::vec& h_unique, const arma::vec& h_u, const arma::mat& dist_all, double c, double nugget, int n_unique, int n_u);
RcppExport SEXP _bdrc_gplm0_predict_u_known_c_cpp(SEXP thetaSEXP, SEXP xSEXP, SEXP h_uniqueSEXP, SEXP h_uSEXP, SEXP dist_allSEXP, SEXP cSEXP, SEXP nuggetSEXP, SEXP n_uniqueSEXP, SEXP n_uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type h_unique(h_uniqueSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type h_u(h_uSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type dist_all(dist_allSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    Rcpp::traits::input_parameter< double >::type nugget(nuggetSEXP);
    Rcpp::traits::input_parameter< int >::type n_unique(n_uniqueSEXP);
    Rcpp::traits::input_parameter< int >::type n_u(n_uSEXP);
    rcpp_result_gen = Rcpp::wrap(gplm0_predict_u_known_c_cpp(theta, x, h_unique, h_u, dist_all, c, nugget, n_unique, n_u));
    return rcpp_result_gen;
END_RCPP
}
// matMult
arma::mat matMult(const arma::mat& A, const arma::mat& B);
RcppExport SEXP _bdrc_matMult(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(matMult(A, B));
    return rcpp_result_gen;
END_RCPP
}
// choleskyDecomp
arma::mat choleskyDecomp(const arma::mat& X);
RcppExport SEXP _bdrc_choleskyDecomp(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(choleskyDecomp(X));
    return rcpp_result_gen;
END_RCPP
}
// solveArma
arma::vec solveArma(const arma::mat& A, const arma::vec& B);
RcppExport SEXP _bdrc_solveArma(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(solveArma(A, B));
    return rcpp_result_gen;
END_RCPP
}
// matInverse
arma::mat matInverse(const arma::mat& A);
RcppExport SEXP _bdrc_matInverse(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(matInverse(A));
    return rcpp_result_gen;
END_RCPP
}
// compute_L
arma::mat compute_L(const arma::mat& X, const arma::mat& Sig_x, const arma::mat& Sig_eps, double nugget);
RcppExport SEXP _bdrc_compute_L(SEXP XSEXP, SEXP Sig_xSEXP, SEXP Sig_epsSEXP, SEXP nuggetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Sig_x(Sig_xSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Sig_eps(Sig_epsSEXP);
    Rcpp::traits::input_parameter< double >::type nugget(nuggetSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_L(X, Sig_x, Sig_eps, nugget));
    return rcpp_result_gen;
END_RCPP
}
// compute_w
arma::vec compute_w(const arma::mat& L, const arma::vec& y, const arma::mat& X, const arma::vec& mu_x);
RcppExport SEXP _bdrc_compute_w(SEXP LSEXP, SEXP ySEXP, SEXP XSEXP, SEXP mu_xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type L(LSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu_x(mu_xSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_w(L, y, X, mu_x));
    return rcpp_result_gen;
END_RCPP
}
// get_MCMC_summary_cpp
DataFrame get_MCMC_summary_cpp(const arma::mat& X, const Nullable<NumericVector>& h);
RcppExport SEXP _bdrc_get_MCMC_summary_cpp(SEXP XSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Nullable<NumericVector>& >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(get_MCMC_summary_cpp(X, h));
    return rcpp_result_gen;
END_RCPP
}
// calc_variogram_arma
arma::vec calc_variogram_arma(const arma::mat& param_mat, int i, int burnin, int nr_iter);
RcppExport SEXP _bdrc_calc_variogram_arma(SEXP param_matSEXP, SEXP iSEXP, SEXP burninSEXP, SEXP nr_iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type param_mat(param_matSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< int >::type nr_iter(nr_iterSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_variogram_arma(param_mat, i, burnin, nr_iter));
    return rcpp_result_gen;
END_RCPP
}
// variogram_chain
arma::mat variogram_chain(int T_max, const arma::mat& param_mat1, const arma::mat& param_mat2, int burnin, int nr_iter);
RcppExport SEXP _bdrc_variogram_chain(SEXP T_maxSEXP, SEXP param_mat1SEXP, SEXP param_mat2SEXP, SEXP burninSEXP, SEXP nr_iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type T_max(T_maxSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type param_mat1(param_mat1SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type param_mat2(param_mat2SEXP);
    Rcpp::traits::input_parameter< int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< int >::type nr_iter(nr_iterSEXP);
    rcpp_result_gen = Rcpp::wrap(variogram_chain(T_max, param_mat1, param_mat2, burnin, nr_iter));
    return rcpp_result_gen;
END_RCPP
}
// distance_matrix
arma::mat distance_matrix(const arma::vec& x);
RcppExport SEXP _bdrc_distance_matrix(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(distance_matrix(x));
    return rcpp_result_gen;
END_RCPP
}
// create_A_cpp
arma::mat create_A_cpp(const arma::vec& h);
RcppExport SEXP _bdrc_create_A_cpp(SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(create_A_cpp(h));
    return rcpp_result_gen;
END_RCPP
}
// pri
double pri(const std::string& type, const arma::vec& args);
RcppExport SEXP _bdrc_pri(SEXP typeSEXP, SEXP argsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type type(typeSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type args(argsSEXP);
    rcpp_result_gen = Rcpp::wrap(pri(type, args));
    return rcpp_result_gen;
END_RCPP
}
// chain_statistics_cpp
Rcpp::List chain_statistics_cpp(const arma::mat& chains);
RcppExport SEXP _bdrc_chain_statistics_cpp(SEXP chainsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type chains(chainsSEXP);
    rcpp_result_gen = Rcpp::wrap(chain_statistics_cpp(chains));
    return rcpp_result_gen;
END_RCPP
}
// plm_density_evaluation_unknown_c_cpp
Rcpp::List plm_density_evaluation_unknown_c_cpp(const arma::vec& theta, const arma::mat& P, const arma::vec& h, const arma::mat& B, const arma::vec& y, const arma::vec& epsilon, const arma::mat& Sig_x, const arma::vec& mu_x, double h_min, double nugget, double lambda_c, double lambda_eta_1, double lambda_seta);
RcppExport SEXP _bdrc_plm_density_evaluation_unknown_c_cpp(SEXP thetaSEXP, SEXP PSEXP, SEXP hSEXP, SEXP BSEXP, SEXP ySEXP, SEXP epsilonSEXP, SEXP Sig_xSEXP, SEXP mu_xSEXP, SEXP h_minSEXP, SEXP nuggetSEXP, SEXP lambda_cSEXP, SEXP lambda_eta_1SEXP, SEXP lambda_setaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type P(PSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type h(hSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type B(BSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Sig_x(Sig_xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu_x(mu_xSEXP);
    Rcpp::traits::input_parameter< double >::type h_min(h_minSEXP);
    Rcpp::traits::input_parameter< double >::type nugget(nuggetSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_c(lambda_cSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_eta_1(lambda_eta_1SEXP);
    Rcpp::traits::input_parameter< double >::type lambda_seta(lambda_setaSEXP);
    rcpp_result_gen = Rcpp::wrap(plm_density_evaluation_unknown_c_cpp(theta, P, h, B, y, epsilon, Sig_x, mu_x, h_min, nugget, lambda_c, lambda_eta_1, lambda_seta));
    return rcpp_result_gen;
END_RCPP
}
// plm_density_evaluation_known_c_cpp
Rcpp::List plm_density_evaluation_known_c_cpp(const arma::vec& theta, const arma::mat& P, const arma::vec& h, const arma::mat& B, const arma::vec& y, const arma::vec& epsilon, const arma::mat& Sig_x, const arma::vec& mu_x, double c, double nugget, double lambda_eta_1, double lambda_seta);
RcppExport SEXP _bdrc_plm_density_evaluation_known_c_cpp(SEXP thetaSEXP, SEXP PSEXP, SEXP hSEXP, SEXP BSEXP, SEXP ySEXP, SEXP epsilonSEXP, SEXP Sig_xSEXP, SEXP mu_xSEXP, SEXP cSEXP, SEXP nuggetSEXP, SEXP lambda_eta_1SEXP, SEXP lambda_setaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type P(PSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type h(hSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type B(BSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Sig_x(Sig_xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu_x(mu_xSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    Rcpp::traits::input_parameter< double >::type nugget(nuggetSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_eta_1(lambda_eta_1SEXP);
    Rcpp::traits::input_parameter< double >::type lambda_seta(lambda_setaSEXP);
    rcpp_result_gen = Rcpp::wrap(plm_density_evaluation_known_c_cpp(theta, P, h, B, y, epsilon, Sig_x, mu_x, c, nugget, lambda_eta_1, lambda_seta));
    return rcpp_result_gen;
END_RCPP
}
// plm_predict_u_unknown_c_cpp
Rcpp::List plm_predict_u_unknown_c_cpp(const arma::vec& theta, const arma::vec& x, const arma::mat& P, const arma::mat& B_u, const arma::vec& h_u, double h_min, int n_u);
RcppExport SEXP _bdrc_plm_predict_u_unknown_c_cpp(SEXP thetaSEXP, SEXP xSEXP, SEXP PSEXP, SEXP B_uSEXP, SEXP h_uSEXP, SEXP h_minSEXP, SEXP n_uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type P(PSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type B_u(B_uSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type h_u(h_uSEXP);
    Rcpp::traits::input_parameter< double >::type h_min(h_minSEXP);
    Rcpp::traits::input_parameter< int >::type n_u(n_uSEXP);
    rcpp_result_gen = Rcpp::wrap(plm_predict_u_unknown_c_cpp(theta, x, P, B_u, h_u, h_min, n_u));
    return rcpp_result_gen;
END_RCPP
}
// plm_predict_u_known_c_cpp
Rcpp::List plm_predict_u_known_c_cpp(const arma::vec& theta, const arma::vec& x, const arma::mat& P, const arma::mat& B_u, const arma::vec& h_u, double c);
RcppExport SEXP _bdrc_plm_predict_u_known_c_cpp(SEXP thetaSEXP, SEXP xSEXP, SEXP PSEXP, SEXP B_uSEXP, SEXP h_uSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type P(PSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type B_u(B_uSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type h_u(h_uSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(plm_predict_u_known_c_cpp(theta, x, P, B_u, h_u, c));
    return rcpp_result_gen;
END_RCPP
}
// plm0_density_evaluation_unknown_c_cpp
Rcpp::List plm0_density_evaluation_unknown_c_cpp(const arma::vec& theta, const arma::vec& h, const arma::vec& y, const arma::vec& epsilon, const arma::mat& Sig_ab, const arma::vec& mu_x, double h_min, double nugget, double lambda_c, double lambda_se);
RcppExport SEXP _bdrc_plm0_density_evaluation_unknown_c_cpp(SEXP thetaSEXP, SEXP hSEXP, SEXP ySEXP, SEXP epsilonSEXP, SEXP Sig_abSEXP, SEXP mu_xSEXP, SEXP h_minSEXP, SEXP nuggetSEXP, SEXP lambda_cSEXP, SEXP lambda_seSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type h(hSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Sig_ab(Sig_abSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu_x(mu_xSEXP);
    Rcpp::traits::input_parameter< double >::type h_min(h_minSEXP);
    Rcpp::traits::input_parameter< double >::type nugget(nuggetSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_c(lambda_cSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_se(lambda_seSEXP);
    rcpp_result_gen = Rcpp::wrap(plm0_density_evaluation_unknown_c_cpp(theta, h, y, epsilon, Sig_ab, mu_x, h_min, nugget, lambda_c, lambda_se));
    return rcpp_result_gen;
END_RCPP
}
// plm0_density_evaluation_known_c_cpp
Rcpp::List plm0_density_evaluation_known_c_cpp(const arma::vec& theta, const arma::vec& h, const arma::vec& y, const arma::vec& epsilon, const arma::mat& Sig_ab, const arma::vec& mu_x, double c, double nugget, double lambda_se);
RcppExport SEXP _bdrc_plm0_density_evaluation_known_c_cpp(SEXP thetaSEXP, SEXP hSEXP, SEXP ySEXP, SEXP epsilonSEXP, SEXP Sig_abSEXP, SEXP mu_xSEXP, SEXP cSEXP, SEXP nuggetSEXP, SEXP lambda_seSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type h(hSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Sig_ab(Sig_abSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu_x(mu_xSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    Rcpp::traits::input_parameter< double >::type nugget(nuggetSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_se(lambda_seSEXP);
    rcpp_result_gen = Rcpp::wrap(plm0_density_evaluation_known_c_cpp(theta, h, y, epsilon, Sig_ab, mu_x, c, nugget, lambda_se));
    return rcpp_result_gen;
END_RCPP
}
// plm0_predict_u_unknown_c_cpp
Rcpp::List plm0_predict_u_unknown_c_cpp(const arma::vec& theta, const arma::vec& x, const arma::vec& h_u, double h_min);
RcppExport SEXP _bdrc_plm0_predict_u_unknown_c_cpp(SEXP thetaSEXP, SEXP xSEXP, SEXP h_uSEXP, SEXP h_minSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type h_u(h_uSEXP);
    Rcpp::traits::input_parameter< double >::type h_min(h_minSEXP);
    rcpp_result_gen = Rcpp::wrap(plm0_predict_u_unknown_c_cpp(theta, x, h_u, h_min));
    return rcpp_result_gen;
END_RCPP
}
// plm0_predict_u_known_c_cpp
Rcpp::List plm0_predict_u_known_c_cpp(const arma::vec& theta, const arma::vec& x, const arma::vec& h_u, double c);
RcppExport SEXP _bdrc_plm0_predict_u_known_c_cpp(SEXP thetaSEXP, SEXP xSEXP, SEXP h_uSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type h_u(h_uSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(plm0_predict_u_known_c_cpp(theta, x, h_u, c));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_bdrc_gplm_density_evaluation_unknown_c_cpp", (DL_FUNC) &_bdrc_gplm_density_evaluation_unknown_c_cpp, 19},
    {"_bdrc_gplm_density_evaluation_known_c_cpp", (DL_FUNC) &_bdrc_gplm_density_evaluation_known_c_cpp, 18},
    {"_bdrc_gplm_predict_u_unknown_c_cpp", (DL_FUNC) &_bdrc_gplm_predict_u_unknown_c_cpp, 11},
    {"_bdrc_gplm_predict_u_known_c_cpp", (DL_FUNC) &_bdrc_gplm_predict_u_known_c_cpp, 11},
    {"_bdrc_gplm0_density_evaluation_unknown_c_cpp", (DL_FUNC) &_bdrc_gplm0_density_evaluation_unknown_c_cpp, 16},
    {"_bdrc_gplm0_density_evaluation_known_c_cpp", (DL_FUNC) &_bdrc_gplm0_density_evaluation_known_c_cpp, 15},
    {"_bdrc_gplm0_predict_u_unknown_c_cpp", (DL_FUNC) &_bdrc_gplm0_predict_u_unknown_c_cpp, 9},
    {"_bdrc_gplm0_predict_u_known_c_cpp", (DL_FUNC) &_bdrc_gplm0_predict_u_known_c_cpp, 9},
    {"_bdrc_matMult", (DL_FUNC) &_bdrc_matMult, 2},
    {"_bdrc_choleskyDecomp", (DL_FUNC) &_bdrc_choleskyDecomp, 1},
    {"_bdrc_solveArma", (DL_FUNC) &_bdrc_solveArma, 2},
    {"_bdrc_matInverse", (DL_FUNC) &_bdrc_matInverse, 1},
    {"_bdrc_compute_L", (DL_FUNC) &_bdrc_compute_L, 4},
    {"_bdrc_compute_w", (DL_FUNC) &_bdrc_compute_w, 4},
    {"_bdrc_get_MCMC_summary_cpp", (DL_FUNC) &_bdrc_get_MCMC_summary_cpp, 2},
    {"_bdrc_calc_variogram_arma", (DL_FUNC) &_bdrc_calc_variogram_arma, 4},
    {"_bdrc_variogram_chain", (DL_FUNC) &_bdrc_variogram_chain, 5},
    {"_bdrc_distance_matrix", (DL_FUNC) &_bdrc_distance_matrix, 1},
    {"_bdrc_create_A_cpp", (DL_FUNC) &_bdrc_create_A_cpp, 1},
    {"_bdrc_pri", (DL_FUNC) &_bdrc_pri, 2},
    {"_bdrc_chain_statistics_cpp", (DL_FUNC) &_bdrc_chain_statistics_cpp, 1},
    {"_bdrc_plm_density_evaluation_unknown_c_cpp", (DL_FUNC) &_bdrc_plm_density_evaluation_unknown_c_cpp, 13},
    {"_bdrc_plm_density_evaluation_known_c_cpp", (DL_FUNC) &_bdrc_plm_density_evaluation_known_c_cpp, 12},
    {"_bdrc_plm_predict_u_unknown_c_cpp", (DL_FUNC) &_bdrc_plm_predict_u_unknown_c_cpp, 7},
    {"_bdrc_plm_predict_u_known_c_cpp", (DL_FUNC) &_bdrc_plm_predict_u_known_c_cpp, 6},
    {"_bdrc_plm0_density_evaluation_unknown_c_cpp", (DL_FUNC) &_bdrc_plm0_density_evaluation_unknown_c_cpp, 10},
    {"_bdrc_plm0_density_evaluation_known_c_cpp", (DL_FUNC) &_bdrc_plm0_density_evaluation_known_c_cpp, 9},
    {"_bdrc_plm0_predict_u_unknown_c_cpp", (DL_FUNC) &_bdrc_plm0_predict_u_unknown_c_cpp, 4},
    {"_bdrc_plm0_predict_u_known_c_cpp", (DL_FUNC) &_bdrc_plm0_predict_u_known_c_cpp, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_bdrc(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
