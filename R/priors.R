priors <- function(country) {

#Prior Parameters
if(country=="Iceland"){
    RC=list()
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

return(list("RC"=RC))
}
