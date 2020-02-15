#' Generalized Bayesian Power Law Model
#'
#' Infers a rating curve for paired measurements of stage and discharge using a generalized power law model described in Hrafnkelsson et al.
#'@param formula formula with name of discharge column in data as response and name of stage column in data as the single covariate.
#'@param data data.frame containing the columns in formula
#'@param W_limits vector of length 2 setting the lower and upper bound of stage values at which a rating curve should be predicted. If NULL, the known value of c or the mle of c will be used as lower bound (depending on the value of the input parameter c) and maximum stage value in data as upper bound.
#'@param country Name of the country the prior parameters should be defined for, default value is "Iceland".
#'@param Wmin Positive numeric value for the lowest stage the user wants to calculate a rating curve. If input is an empty string (default) Wmin will
#'automatically be set to c_hat.
#'@param Wmax Positive numeric value for the highest stage the user wants to calculate a rating curve. If input is an empty string (default) Wmax will
#'automatically be set to the maximum stage of the data.
#'@return List containing information on the calculated rating curve,
#'the data frames observedData, betaData, completePrediction, observedPrediction, TableOfData, FitTable, LowerTable, UpperTable, plotTable.
#'@references Birgir Hrafnkelsson, Helgi Sigurdarson and Sigurdur M. Gardarson (2015) \emph{Bayesian Generalized Rating Curves}
#'@seealso \code{\link{clean}}
gbplm <- function(formula,data,c=NULL,W_limits=c(0,0),country="Iceland"){
    suppressPackageStartupMessages(require(doParallel))
    #TODO: add error message if length(formula)!=3 or if it contains more than one covariate. Also make sure that names in formula exist in data
    model_dat <- data[,all.vars(formula)]
    RC$formula = formula
    RC=priors(country)
    RC$nugget=10^-8
    RC$mu_sb=0.5
    RC$mu_pb=0.5
    RC$tau_pb2=0.25^2
    RC$s=3
    RC$v=5
    RC$y=rbind(as.matrix(log(model_dat[,1])),0)
    RC$w=as.matrix(model_dat[,2])
    RC$w_tild=RC$w-min(RC$w)
    Adist1 <- Adist(RC$w)
    RC$A=Adist1$A
    RC$dist=Adist1$dist
    RC$n=Adist1$n
    RC$N=Adist1$N
    RC$O=Adist1$O

    RC$Sig_ab= rbind(c(RC$sig_a^2, RC$p_ab*RC$sig_a*RC$sig_b), c(RC$p_ab*RC$sig_a*RC$sig_b, RC$sig_b^2))
    RC$mu_x=as.matrix(c(RC$mu_a,RC$mu_b, rep(0,RC$n)))

    RC$P=diag(nrow=5,ncol=5,6)-matrix(nrow=5,ncol=5,1)
    RC$B=B_splines(t(RC$w_tild)/RC$w_tild[length(RC$w_tild)])
    RC$epsilon=rep(1,RC$N)
    forceIndex=which('forcepoint'== observedData$Quality)
    forcepoint=model_dat[forceIndex,]
    if(any('forcepoint'== observedData$Quality)){
        RC$epsilon[forceIndex]=1/RC$N
    }

    RC$Z=cbind(t(rep(0,2)),t(rep(1,RC$n)))
    RC$m1=matrix(0,nrow=2,ncol=RC$n)
    RC$m2=matrix(0,nrow=RC$n,ncol=2)
    if(!is.null(c)){
        RC$c=c
        density_fun=density_evaluation_known_c
    }else{
        density_fun=density_evaluation_unknown_c
    }
    #determine proposal density
    theta_init=rep(0,9)
    loss_fun = function(th) {-density_fun(th,RC)$p}
    optim_obj=optim(par=theta_init,loss_fun,method="L-BFGS-B",hessian=TRUE)
    t_m =optim_obj$par
    H=optim_obj$hessian
    LH=t(chol(H))/0.8

    #make Wmin and Wmax divisable by 10 up, both in order to make rctafla and so l_m is defined
    if(is.null(W_limits)){
        Wmax=ceiling(max(RC$w)*10)/10
        Wmin=ceiling(10*ifelse(is.null(RC$c),min(RC$w)-exp(t_m[1]),RC$c))/10
    }else{
        Wmin=W_limits[1]
        Wmax=W_limits[2]
    }
    WFill=W_unobserved(c(RC$O),min=Wmin,max=Wmax)
    RC$W_u=WFill$W_u
    RC$W_u_tild=WFill$W_u_tild
    Bsiminput=t(RC$W_u_tild)/RC$W_u_tild[length(RC$W_u_tild)]
    Bsiminput[is.na(Bsiminput)]=0
    RC$Bsim=B_splines(Bsiminput)

    #MCMC parameters added, number of iterations,burnin and thin
    Nit=20000
    burnin=2000
    thin=5
    cl <- makeCluster(4)
    registerDoParallel(cl)
    MCMC <- foreach(i=1:4,.combine=cbind,.export=c("density_fun","predict_u")) %dopar% {
        ypo_obs=matrix(0,nrow=RC$N,ncol=Nit)
        param=matrix(0,nrow=length(t_m)+RC$n+2,ncol=Nit)
        t_old=as.matrix(t_m)
        Dens<-density_fun(t_old,RC)
        p_old=Dens$p
        ypo_old=Dens$ypo
        x_old=Dens$x

        for(j in 1:Nit){
            t_new=t_old+solve(t(LH),rnorm(9,0,1))
            Densnew <- density_fun(t_new,RC)
            x_new=Densnew$x
            ypo_new=Densnew$ypo
            p_new=Densnew$p
            logR=p_new-p_old

            if (logR>log(runif(1))){
                t_old=t_new
                p_old=p_new
                ypo_old=ypo_new
                x_old=x_new


            }
            ypo_obs[,j]=ypo_old
            param[,j]=rbind(t_old,x_old)
        }
        seq=seq(burnin,Nit,thin)
        ypo_obs=ypo_obs[,seq]
        param=param[,seq]
        unobserved=apply(param,2,FUN=function(x) predict_u(x,RC))
        output=rbind(ypo_obs,unobserved)

        return(output)
    }
    #TODO: create S3 object to store results from MCMC chain
    stopCluster(cl)
    MCMC[is.na(MCMC)]=-1000
    betasamples=apply(MCMC[(RC$N+length(RC$W_u)+1):nrow(MCMC),],2,FUN=function(x){x[2]+x[3:length(x)]})
    yposamples=MCMC[1:(RC$N+length(RC$W_u)),]
    completePrediction=as.data.frame(t(apply(yposamples,1,quantile, probs = c(0.025,0.5, 0.975),na.rm=T)))
    names(completePrediction)=c("lower","fit","upper")
    completePrediction$W=c(RC$w,RC$W_u)
    completePrediction$l_m=c(l,log(RC$W_u-min(RC$O)+exp(t_m[1])))
    observedPrediction=completePrediction[1:RC$N,]
    completePrediction=completePrediction[with(completePrediction,order(W)),]
    betaData=as.data.frame(t(apply(betasamples,1,quantile, probs = c(0.025,0.5, 0.975),na.rm=T)))
    names(betaData)=c("lower","fit","upper")
    betaData$W=c(RC$O,RC$W_u)
    betaData=betaData[with(betaData,order(W)),]
    observedPrediction$Q=RC$y[1:RC$N,]
    observedPrediction$residuals=(exp(observedPrediction$Q)-exp(observedPrediction$fit))
    observedPrediction$residupper=exp(observedPrediction$upper)-exp(observedPrediction$fit)
    observedPrediction$residlower=exp(observedPrediction$lower)-exp(observedPrediction$fit)
    observedPrediction$standardResiduals=(observedPrediction$Q-observedPrediction$fit)/sqrt(varr_m)

    TableOfData=observedData
    TableOfData$Q=round(TableOfData$Q,1)
    TableOfData$Qfit=round(exp(observedPrediction$fit),3)
    TableOfData$lower=round(exp(observedPrediction$lower),3)
    TableOfData$upper=round(exp(observedPrediction$upper),3)
    TableOfData$diffQ=TableOfData$Q-TableOfData$Qfit
    TableOfData$Qpercentage=round(100*TableOfData$diffQ/TableOfData$Q,1)
    names(TableOfData)=c("Date","Time","Quality","W","Q", "Q fit","Lower", "Upper","Q diff","Q%")
    TableOfData=TableOfData[with(TableOfData,order(Date)),]

    xout=seq(Wmin,-0.01+Wmax,by=0.01)

    fitInterpolation=approx(completePrediction$W,completePrediction$fit,xout=xout)
    FitTable=t(as.data.frame(split(x=fitInterpolation$y, f=ceiling(seq_along(fitInterpolation$y)/10))))
    colnames(FitTable)=0:9
    FitTable=round(exp(FitTable),3)
    Stage=seq(min(fitInterpolation$x),max(fitInterpolation$x),by=0.1)*100
    FitTable=as.data.frame(cbind(Stage,FitTable))
    names(FitTable)[1]="Stage (cm)"

    lowerInterpolation=approx(completePrediction$W,completePrediction$lower,xout=xout)
    LowerTable=t(as.data.frame(split(x=lowerInterpolation$y, f=ceiling(seq_along(lowerInterpolation$y)/10))))
    colnames(LowerTable)=0:9
    LowerTable=round(exp(LowerTable),3)
    Stage=seq(min(lowerInterpolation$x),max(lowerInterpolation$x),by=0.1)*100
    LowerTable=as.data.frame(cbind(Stage,LowerTable))
    names(LowerTable)[1]="Stage (cm)"

    upperInterpolation=approx(completePrediction$W,completePrediction$upper,xout=xout)
    UpperTable=t(as.data.frame(split(x=upperInterpolation$y, f=ceiling(seq_along(upperInterpolation$y)/10))))
    colnames(UpperTable)=0:9
    UpperTable=round(exp(UpperTable),3)
    Stage=seq(min(upperInterpolation$x),max(upperInterpolation$x),by=0.1)*100
    UpperTable=as.data.frame(cbind(Stage,UpperTable))
    names(UpperTable)[1]="Stage (cm)"

    plotTable=as.data.frame(cbind(lowerInterpolation$y,fitInterpolation$y,upperInterpolation$y))
    plotTable=exp(plotTable)
    names(plotTable)=c("Lower","Fit","Upper")
    plotTable$W=xout
    
    
    #S3 object gbplm Test
    
    print.gbplm <- function(x,...){
      cat("\nCall:\n",
              paste(deparse(x$formula), sep = "\n", collapse = "\n"), "\n\n", sep = "")
    }
    
    
    
    return(list("observedData"=observedData,"betaData"=betaData,"completePrediction"=completePrediction,"observedPrediction"=observedPrediction,"TableOfData"=TableOfData,
                "FitTable"=FitTable,"LowerTable"=LowerTable,"UpperTable"=UpperTable,"plotTable"=plotTable))
}

#'Density evaluation for model2
#'
#'Evaluates the log density of the posterior distribution of the parameters of .
#'@param th A vector with length 9 containing parameters
#'@param RC A list containing prior parameters, matrices and the data.
#'@return Returns a list containing predictive values of the parameters drawn out of the evaluated density.
#'@references Birgir Hrafnkelsson, Helgi Sigurdarson and Sigurdur M. Gardarson (2015) \emph{Bayesian Generalized Rating Curves}
density_evaluation_known_c <- function(th,RC){
    phi_b=th[2]
    sig_b2=th[1]
    lambda=th[4:8]

    f=lambda[1:5]-lambda[6]
    l=log(RC$w-RC$c)

    varr=c(RC$epsilon*exp(RC$B%*%lambda))
    Sig_eps=diag(c(varr,0))
    #Matern covariance
    R_Beta=(1+sqrt(5)*RC$dist/exp(phi_b)+5*RC$dist^2/(3*exp(phi_b)^2))*exp(-sqrt(5)*RC$dist/exp(phi_b))+diag(RC$n)*RC$nugget
    Sig_x=rbind(cbind(RC$Sig_ab,RC$m1),cbind(RC$m2,exp(sig_b2)*R_Beta))

    X=rbind(cbind(1,l,diag(l)%*%RC$A),RC$Z)#1337 micro
    L=t(chol(X%*%Sig_x%*%t(X)+Sig_eps))#2521 micro
    w=solve(L,RC$y-X%*%RC$mu_x)#754 micro
    p=-0.5%*%t(w)%*%w-sum(log(diag(L)))-
        (RC$v+5-1)/2*log(RC$v*RC$s+f%*%RC$P%*%f)+
        sig_b2-exp(sig_b2)/RC$mu_sb-0.5/RC$tau_pb2*(phi_b-RC$mu_pb)^2 #63 micro

    W=solve(L,X%*%Sig_x)
    x_u=RC$mu_x+t(chol(Sig_x))%*%rnorm(RC$n+2)
    sss=(X%*%x_u)-RC$y+rbind(sqrt(varr)*as.matrix(rnorm(RC$N)),0)
    x=as.matrix(x_u-t(W)%*%solve(L,sss))
    yp=(X %*% x)[1:RC$N,]
    #posterior predictive draw
    ypo=yp+as.matrix(rnorm(RC$N))*sqrt(varr)

    return(list("p"=p,"x"=x,"yp"=yp,"ypo"=ypo,"varr"=varr))
}

#'Density evaluation for model2
#'
#'Evaluates the log density of the posterior distribution of the parameters of model2BH.
#'@param th A vector with length 9 containing parameters
#'@param RC A list containing prior parameters, matrices and the data.
#'@return Returns a list containing predictive values of the parameters drawn out of the evaluated density.
#'@references Birgir Hrafnkelsson, Helgi Sigurdarson and Sigurdur M. Gardarson (2015) \emph{Bayesian Generalized Rating Curves}
density_evaluation_unknown_c <- function(th,RC){
    phi_b=th[3]
    sig_b2=th[2]
    zeta=th[1]
    lambda=th[4:9]

    f=lambda[1:5]-lambda[6]
    l=c(log(RC$w_tild+exp(th[1])))

    varr=c(RC$epsilon*exp(RC$B%*%lambda))
    Sig_eps=diag(c(varr,0))
    #Matern covariance
    R_Beta=(1+sqrt(5)*RC$dist/exp(phi_b)+5*RC$dist^2/(3*exp(phi_b)^2))*exp(-sqrt(5)*RC$dist/exp(phi_b))+diag(RC$n)*RC$nugget
    Sig_x=rbind(cbind(RC$Sig_ab,RC$m1),cbind(RC$m2,exp(sig_b2)*R_Beta))

    X=rbind(cbind(1,l,diag(l)%*%RC$A),RC$Z)#1337 micro
    L=t(chol(X%*%Sig_x%*%t(X)+Sig_eps))#2521 micro
    w=solve(L,RC$y-X%*%RC$mu_x)#754 micro
    p=-0.5%*%t(w)%*%w-sum(log(diag(L)))-
        (RC$v+5-1)/2*log(RC$v*RC$s+f%*%RC$P%*%f)+
        sig_b2-exp(sig_b2)/RC$mu_sb+zeta-exp(zeta)/RC$mu_c-0.5/RC$tau_pb2*(phi_b-RC$mu_pb)^2 #63 micro

    W=solve(L,X%*%Sig_x)
    x_u=RC$mu_x+t(chol(Sig_x))%*%rnorm(RC$n+2)
    sss=(X%*%x_u)-RC$y+rbind(sqrt(varr)*as.matrix(rnorm(RC$N)),0)
    x=as.matrix(x_u-t(W)%*%solve(L,sss))
    yp=(X %*% x)[1:RC$N,]
    #posterior predictive draw
    ypo=yp+as.matrix(rnorm(RC$N))*sqrt(varr)

    #D=-2*sum(log(dlnorm(exp(RC$y[1:RC$N,]),yp,sqrt(varr))))#45.04

    return(list("p"=p,"x"=x,"yp"=yp,"ypo"=ypo,"varr"=varr))
}


