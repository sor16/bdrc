Densevalm22_u <- function(param,RC){
    #collecting parameters from the MCMC sample
    th=param[1:9]
    x=param[10:length(param)]
    zeta=th[1]
    phi_b=exp(th[3])
    sig_b2=exp(th[2])
    lambda = th[4:9]
    #calculate spline variance from B_splines
    varr = c(exp(RC$Bsim %*% lambda))
    m=length(RC$W_u)
    n=RC$n
    #combine stages from data with unobserved stages
    W_all=c(RC$O,RC$W_u)
    #calculating distance matrix for W_all
    dist=abs(outer(W_all,W_all,FUN="-"))
    #defining the variance of the joint prior for betas from data and beta unobserved, that is p(beta,beta_u).
    #Matern covariance formula used for v=5/2
    sigma_all=sig_b2*(1 + sqrt(5)*dist/phi_b+(5*dist^2)/(3*phi_b^2))*exp(-sqrt(5)*dist/phi_b) + diag(length(W_all))*RC$nugget
    sigma_11=sigma_all[1:n,1:n]
    sigma_22=sigma_all[(n+1):(m+n),(n+1):(m+n)]
    sigma_12=sigma_all[1:n,(n+1):(n+m)]
    sigma_21=sigma_all[(n+1):(n+m),1:n]
    #parameters for the posterior of beta_u
    mu_u=sigma_21%*%solve(sigma_11,x[3:length(x)])
    Sigma_u=(sigma_22-sigma_21%*%solve(sigma_11,sigma_12))
    #a sample from posterior of beta_u drawn
    beta_u=as.numeric(mu_u) + rnorm(ncol(Sigma_u)) %*% chol(Sigma_u)
    #buidling blocks of the explanatory matrix X calculated
    c=min(RC$O)-exp(zeta)
    l=log(RC$W_u-c)
    X=cbind(rep(1,m),l,matrix(0,m,n),diag(l))
    #vector of parameters
    x=c(x,beta_u)
    #sample from the posterior of discharge y
    ypo = X%*%x + as.matrix(rnorm(m)) * sqrt(varr)
    return(c(ypo,x))
}
