Densp2 <- function(th,RC){
  phi_b=th[3]
  sig_b2=th[2]
  zeta=th[1]
  lambda=th[4:9]
  
  f=lambda[1:5]-lambda[6]
  l=c(log(RC$w_tild+exp(th[1])))
  
  varr=c(exp(RC$B%*%lambda))
  Sig_eps=diag(c(varr,0))
  R_Beta=(1+sqrt(5)*RC$dist/exp(phi_b)+5*RC$dist^2/(3*exp(phi_b)^2))*exp(-sqrt(5)*RC$dist/exp(phi_b))+diag(RC$n)*RC$nugget
  Sig_x=rbind(cbind(RC$Sig_ab,RC$m1),cbind(RC$m2,exp(sig_b2)*R_Beta))
  
  X=rbind(cbind(1,l,diag(l)%*%RC$A),RC$Z) 
  L=t(chol(X%*%Sig_x%*%t(X)+Sig_eps))
  w=solve(L,RC$y-X%*%RC$mu_x)
  p=as.numeric(-0.5%*%t(w)%*%w-sum(log(diag(L)))-
                 (RC$v+5-1)/2*log(RC$v*RC$s+f%*%RC$P%*%f)+
                 sig_b2-exp(sig_b2)/RC$mu_sb+zeta-exp(zeta)/RC$mu_c-0.5/RC$tau_pb2*(phi_b-RC$mu_pb)^2) 
  
  return(p)
}