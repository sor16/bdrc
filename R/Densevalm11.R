#'Evaluates the log density p at t for model 1
#'
#'Evaluates the log density log(p(t|y)) using p(t|y) propper to p(y|x,t)p(x|t)p(t)/p(x|y,t), which is independent of x
#'so x is set x=c(0 0).
#'@param th A matrix with dim 2x1 containing parameters
#'@param RC RC is a list containing prior parameters, matrices and the data.
#'@return Returns a list containing pmin, minimum p value, and p,x,yp,ypo which are drawn out of the evaluated density.
Densevalm11 <- function(th,RC){

  #hugsanlega onnur breytunofn
  l=log(RC$w_tild+exp(th[1]))
  X=cbind(rep(1,length(l)),l)

  L=t(chol(RC$Sig_xinv + (t(X) %*% X)/exp(th[2])))

  q=solve(L,(RC$Sinvmu+t(X)%*%RC$y/exp(th[2])))
  #Solvi end

  #Solvi begin 27.mai
  p=0.5*sum(q^2)+log(L[1,1])+log(L[2,2])-
  0.5*sum(RC$y^2)/exp(th[2])-RC$n*th[2]/2 +
  th[1]-exp(th[1])*RC$mu_c-th[2]

  x=solve(t(L),(q+as.matrix(rnorm(2))))

  yp=X %*% x

  ypo=yp+as.matrix(rnorm(RC$n,sd=sqrt(exp(th[2]))))

  D=-2*sum(log(dlnorm(exp(RC$y),X%*%x,sqrt(exp(th[2])))))

  return(list("pmin"=-p,"p"=p,"x"=x,"yp"=yp,"ypo"=ypo,"D"=D))
}
