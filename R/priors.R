#' Prior parameter specification
#'
#' Specifies the prior parameters into \code{\link{model1BH}} and \code{\link{model2BH}}
#'@param country A string with the name of the country of which the prior parameters into the models should be specified for
#'@return Returns a list of prior parameters into \code{\link{model1BH}} and \code{\link{model2BH}}.
#'The priors are based from data from rivers of a given country.If you want to add your country to this function,
#' please contact the developers at sor16@@hi.is or aoj8@@hi.is.
#'@references Birgir Hrafnkelsson, Helgi Sigurdarson and Sigurdur M. Gardarson (2015) \emph{Bayesian Generalized Rating Curves}
priors <- function(country) {

    RC=list()
#Prior Parameters
if(country=="Iceland"){
    RC$mu_a=3.20;
    RC$mu_b=2.29;
    RC$sig_a=sqrt(1.21);
    RC$sig_b=sqrt(0.48);
    RC$p_ab=-0.61;
    RC$mu_c=1.9;
}

#Calculations
RC$Sig_x=rbind(c(RC$sig_a^2, RC$p_ab*RC$sig_a*RC$sig_b), c(RC$p_ab*RC$sig_a*RC$sig_b, RC$sig_b^2))
RC$mu_x=as.matrix(c(RC$mu_a, RC$mu_b))
RC$Sig_xinv=solve(RC$Sig_x);
RC$Lx=chol(RC$Sig_x);
RC$Sinvmu=RC$Sig_xinv%*%RC$mu_x;

return(RC)
}
