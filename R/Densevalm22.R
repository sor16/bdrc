#'Evaluates the log density p at t for model 2
#'
#'Evaluates the log density log(p(t|y)) using p(t|y) propper to p(y|x,t)p(x|t)p(t)/p(x|y,t), which is independent of x
#'so x is set x=c(0 0).
#'@param th A vector with length 9 containing parameters
#'@param RC RC is a list containing prior parameters, matrices and the data.
#'@return Returns a list containing p,x,yp,ypo which are drawn out of the evaluated density and varr is a variance numeric vector.
#'@references Birgir Hrafnkelsson, Helgi Sigurdarson and Sigurdur M. Gardarson (2015) \emph{Bayesian Generalized Rating Curves}
Densevalm22 <- function(th,RC){
  phi_b=th[3]
  sig_b2=th[2]
  zeta=th[1]
  lambda=th[4:9]

  f=lambda[1:5]-lambda[6]
  l=c(log(RC$w_tild+exp(th[1])))

  varr=c(RC$epsilon*exp(RC$B%*%lambda))
  Sig_eps=diag(c(varr,0))
  R_Beta=(1+sqrt(5)*RC$dist/exp(phi_b)+5*RC$dist^2/(3*exp(phi_b)^2))*exp(-sqrt(5)*RC$dist/exp(phi_b))+diag(RC$n)*RC$nugget
  Sig_x=rbind(cbind(RC$Sig_ab,RC$m1),cbind(RC$m2,exp(sig_b2)*R_Beta))

  X=rbind(cbind(1,l,diag(l)%*%RC$A),RC$Z)#1337 micro
  L=t(chol(X%*%Sig_x%*%t(X)+Sig_eps))#2521 micro
  w=solve(L,RC$y-X%*%RC$mu_x)#754 micro
  p=-0.5%*%t(w)%*%w-sum(log(diag(L)))-
                 (RC$v+5-1)/2*log(RC$v*RC$s+f%*%RC$P%*%f)+
                 sig_b2-exp(sig_b2)/RC$mu_sb+zeta-exp(zeta)/RC$mu_c-0.5/RC$tau_pb2*(phi_b-RC$mu_pb)^2 #63 micro

  W=solve(L,X%*%Sig_x)#5618 micro
  x_u=RC$mu_x+t(chol(Sig_x))%*%rnorm(RC$n+2)
  sss=(X%*%x_u)-RC$y+rbind(sqrt(varr)*as.matrix(rnorm(RC$N)),0)#94
  x=as.matrix(x_u-t(W)%*%solve(L,sss)) #818 micro
  yp=X %*% x #24 micro
  yp2=yp[1:RC$N,] #7 micro
  ypo=yp2+as.matrix(rnorm(RC$N))*sqrt(varr)#60 micro

  #D=-2*sum(log(dlnorm(exp(RC$y[1:RC$N,]),yp,sqrt(varr))))#45.04

  return(list("p"=p,"x"=x,"yp"=yp,"ypo"=ypo,"varr"=varr))


}
