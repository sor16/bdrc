% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gplm.R
\name{gplm}
\alias{gplm}
\title{Generalized power-law model with variance that varies with stage.}
\usage{
gplm(
  formula,
  data,
  c_param = NULL,
  h_max = NULL,
  parallel = TRUE,
  num_cores = NULL,
  forcepoint = rep(FALSE, nrow(data)),
  verbose = TRUE
)
}
\arguments{
\item{formula}{An object of class "formula", with discharge column name as response and stage column name as a covariate, i.e. of the form \code{y}~\code{x} where \code{y} is discharge in m\eqn{^3/}s and \code{x} is stage in m (it is very important that the data is in the correct units).}

\item{data}{A data.frame containing the variables specified in formula.}

\item{c_param}{The largest stage value for which there is zero discharge. If NULL, it is treated as unknown in the model and inferred from the data.}

\item{h_max}{The maximum stage to which the rating curve should extrapolate to. If NULL, the maximum stage value in the data is selected as an upper bound.}

\item{parallel}{A logical value indicating whether to run the MCMC in parallel or not. Defaults to TRUE.}

\item{num_cores}{An integer between 1 and 4 (number of MCMC chains) indicating how many cores to use. Only used if parallel=TRUE. If NULL, the number of cores available on the device is detected automatically.}

\item{forcepoint}{A logical vector of the same length as the number of rows in data. If an element at index \eqn{i} is TRUE it indicates that the rating curve should be forced through the \eqn{i}-th measurement. Use with care, as this will strongly influence the resulting rating curve.}

\item{verbose}{A logical value indicating whether to print progress and diagnostic information. If `TRUE`, the function will print messages as it runs. If `FALSE`, the function will run silently. Default is `TRUE`.}
}
\value{
gplm returns an object of class "gplm". An object of class "gplm" is a list containing the following components:
\describe{
  \item{\code{rating_curve}}{A data frame with 2.5\%, 50\% and 97.5\% percentiles of the posterior predictive distribution of the rating curve.}
  \item{\code{rating_curve_mean}}{A data frame with 2.5\%, 50\% and 97.5\% percentiles of the posterior distribution of the mean of the rating curve.}
  \item{\code{param_summary}}{A data frame with 2.5\%, 50\% and 97.5\% percentiles of the posterior distribution of latent- and hyperparameters. Additionally contains columns with r_hat and the effective number of samples for each parameter as defined in Gelman et al. (2013).}
  \item{\code{f_summary}}{A data frame with 2.5\%, 50\% and 97.5\% percentiles of the posterior distribution of \eqn{f(h)}.}
  \item{\code{beta_summary}}{A data frame with 2.5\%, 50\% and 97.5\% percentiles of the posterior distribution of \eqn{\beta(h)}.}
  \item{\code{sigma_eps_summary}}{A data frame with 2.5\%, 50\% and 97.5\% percentiles of the posterior distribution of \eqn{\sigma_\varepsilon(h)}.}
  \item{\code{posterior_log_likelihood_summary}}{A data frame with 2.5\%, 50\% and 97.5\% percentiles of the posterior log-likelihood values.}
  \item{\code{rating_curve_posterior}}{A matrix containing the full thinned posterior samples of the posterior predictive distribution of the rating curve excluding burn-in samples.}
  \item{\code{rating_curve_mean_posterior}}{A matrix containing the full thinned posterior samples of the posterior distribution of the mean of the rating curve excluding burn-in samples.}
  \item{\code{a_posterior}}{A numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{a} excluding burn-in samples.}
  \item{\code{b_posterior}}{A numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{b} excluding burn-in samples.}
  \item{\code{c_posterior}}{A numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{c} excluding burn-in samples.}
  \item{\code{sigma_beta_posterior}}{A numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{\sigma_\beta} excluding burn-in samples.}
  \item{\code{phi_beta_posterior}}{A numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{\phi_\beta} excluding burn-in samples.}
  \item{\code{sigma_eta_posterior}}{A numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{\sigma_\eta} excluding burn-in samples.}
  \item{\code{eta_1_posterior}}{A numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{\eta_1} excluding burn-in samples.}
  \item{\code{eta_2_posterior}}{A numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{\eta_2} excluding burn-in samples.}
  \item{\code{eta_3_posterior}}{A numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{\eta_3} excluding burn-in samples.}
  \item{\code{eta_4_posterior}}{A numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{\eta_4} excluding burn-in samples.}
  \item{\code{eta_5_posterior}}{A numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{\eta_5} excluding burn-in samples.}
  \item{\code{eta_6_posterior}}{A numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{\eta_6} excluding burn-in samples.}
  \item{\code{f_posterior}}{A numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{f(h)} excluding burn-in samples.}
  \item{\code{beta_posterior}}{A numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{\beta(h)} excluding burn-in samples.}
  \item{\code{sigma_eps_posterior}}{A numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{\sigma_\varepsilon(h)} excluding burn-in samples.}
  \item{\code{posterior_log_likelihood}}{A numeric vector containing the full thinned posterior log-likelihood values, excluding burn-in samples.}
  \item{\code{D_hat}}{A statistic defined as -2 times the log-likelihood evaluated at the median value of the parameters.}
  \item{\code{effective_num_param_DIC}}{The effective number of parameters, which is calculated as median(-2*posterior_log_likelihood) minus D_hat.}
  \item{\code{DIC}}{The Deviance Information Criterion for the model, calculated as D_hat plus 2*effective_num_parameters_DIC.}
  \item{\code{lppd}}{The log pointwise predictive density of the observed data under the model.}
  \item{\code{WAIC}}{The Widely Applicable Information Criterion for the model, defined as -2*( lppd - effective_num_param_WAIC ).}
  \item{\code{WAIC_i}}{The pointwise WAIC values, where WAIC := sum(WAIC_i).}
  \item{\code{effective_num_param_WAIC}}{The effective number of parameters, which is calculated by summing up the posterior variance of the log predictive density for each data point.}
  \item{\code{autocorrelation}}{A data frame with the autocorrelation of each parameter for different lags.}
  \item{\code{acceptance_rate}}{The proportion of accepted samples in the thinned MCMC chain (excluding burn-in).}
  \item{\code{formula}}{An object of type "formula" provided by the user.}
  \item{\code{data}}{The data provided by the user, ordered by stage.}
  \item{\code{run_info}}{The information about the input arguments and the specific parameters used in the MCMC chain.}
}
}
\description{
gplm is used to fit a discharge rating curve for paired measurements of stage and discharge using a generalized power-law model with variance that varies with stage as described in Hrafnkelsson et al. (2022).  See "Details" for a more elaborate description of the model.
}
\details{
The generalized power-law model is of the form
\deqn{Q=a(h-c)^{f(h)}}
where \eqn{Q} is discharge, \eqn{h} is stage, \eqn{a} and \eqn{c} are unknown constants and \eqn{f} is a function of \eqn{h}, referred to as the generalized power-law exponent.\cr\cr
The generalized power-law model is here inferred by using a Bayesian hierarchical model. The function \eqn{f} is modeled at the latent level as a fixed constant \eqn{b} plus a continuous stochastic process, \eqn{\beta(h)}, which is assumed to be twice differentiable. The model is on a logarithmic scale
\deqn{\log(Q_i) = \log(a) + (b + \beta(h_i)) \log(h_i - c) + \varepsilon_i,     i = 1,...,n}
where \eqn{\varepsilon_i} follows a normal distribution with mean zero and variance \eqn{\sigma_\varepsilon(h_i)^2} that varies with stage. The stochastic process \eqn{\beta(h)} is assumed a priori to be a Gaussian process governed by a Matern covariance function with smoothness parameter \eqn{\nu = 2.5}. The error variance, \eqn{\sigma_\varepsilon^2(h)}, of the log-discharge data  is modeled as an exponential of a B-spline curve, that is, a linear combination of six B-spline basis functions that are defined over the range of the stage observations. An efficient posterior simulation is achieved by sampling from the joint posterior density of the hyperparameters of the model, and then sampling from the density of the latent parameters conditional on the hyperparameters.\cr\cr
Bayesian inference is based on the posterior density and summary statistics such as the posterior mean and 95\% posterior intervals are based on the posterior density. Analytical formulas for these summary statistics are intractable in most cases and thus they are computed by generating samples from the posterior density using a Markov chain Monte Carlo simulation.
}
\examples{
\donttest{
data(norn)
set.seed(1)
gplm.fit <- gplm(formula=Q~W,data=norn,num_cores=2)
summary(gplm.fit)
}
}
\references{
Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A., and Rubin, D. B. (2013). Bayesian Data Analysis, Third Edition. Chapman & Hall/CRC Texts in Statistical Science. Taylor & Francis. doi: https://doi.org/10.1201/b16018

Hrafnkelsson, B., Sigurdarson, H., Rögnvaldsson, S., Jansson, A. Ö., Vias, R. D., and Gardarsson, S. M. (2022). Generalization of the power-law rating curve using hydrodynamic theory and Bayesian hierarchical modeling, Environmetrics, 33(2):e2711. doi: https://doi.org/10.1002/env.2711

Spiegelhalter, D., Best, N., Carlin, B., Van Der Linde, A. (2002). Bayesian measures of model complexity and fit. Journal of the Royal Statistical Society: Series B (Statistical Methodology) 64(4), 583–639. doi: https://doi.org/10.1111/1467-9868.00353

Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and widely applicable information criterion in singular learning theory. J. Mach. Learn. Res. 11, 3571–3594.
}
\seealso{
\code{\link{summary.gplm}} for summaries, \code{\link{predict.gplm}} for prediction and \code{\link{plot.gplm}} for plots. \code{\link{spread_draws}} and \code{\link{gather_draws}} are also useful to aid further visualization of the full posterior distributions.
}
