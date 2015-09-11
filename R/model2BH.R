#'Calculates the rating curve using model2
#'
#'This function takes in clean data from the \code{\link{clean}} function
#'and calculates the rating curve using model 2.
#'@param clean List that is the output of \code{\link{clean}} function.
#'@param country Name of the country the prior parameters should be defined for, default value is "Iceland".
#'@param Wmin Positive numeric value for the lowest stage the user wants to calculate a rating curve. If input is an empty string (default) Wmin will
#'automatically be set to c_hat.
#'@param Wmax Positive numeric value for the highest stage the user wants to calculate a rating curve. If input is an empty string (default) Wmax will
#'automatically be set to the maximum stage of the data.
#'@return List containing information on the calculated rating curve,
#'the data frames observedData, betaData, completePrediction, observedPrediction, TableOfData, FitTable, LowerTable, UpperTable, plotTable.
#'@references Birgir Hrafnkelsson, Helgi Sigurdarson and Sigurdur M. Gardarson (2015) \emph{Bayesian Generalized Rating Curves}
#'@seealso \code{\link{clean}}
model2BH <- function(clean,country="Iceland",Wmin="",Wmax=""){
    suppressPackageStartupMessages(require(doParallel))
    require(RCmodels)
    list2env(clean,envir=environment())

    #MCMC parameters added, number of iterations,burnin and thin
    Nit=20000
    burnin=2000
    thin=5

    RC=priors(country)
    RC$nugget=10^-8
    RC$mu_sb=0.5
    RC$mu_pb=0.5
    RC$tau_pb2=0.25^2
    RC$s=3
    RC$v=5

    forceindex=which('forcepoint'== observedData$Quality)
    forcepoint=wq[forceindex,]

    RC$y=rbind(as.matrix(log(wq[,2])),0)
    RC$w=as.matrix(wq[,1])
    RC$w_tild=RC$w-min(RC$w)

    Adist1 <- Adist(RC$w)
    RC$A=Adist1$A
    RC$dist=Adist1$dist
    RC$n=Adist1$n
    RC$N=Adist1$N
    RC$O=Adist1$O

    RC$P=diag(nrow=5,ncol=5,6)-matrix(nrow=5,ncol=5,1)
    RC$Sig_ab= rbind(c(RC$sig_a^2, RC$p_ab*RC$sig_a*RC$sig_b), c(RC$p_ab*RC$sig_a*RC$sig_b, RC$sig_b^2))
    RC$mu_x=as.matrix(c(RC$mu_a,RC$mu_b, rep(0,RC$n))) #Setja i RC

    RC$B=B_splines(t(RC$w_tild)/RC$w_tild[length(RC$w_tild)])
    RC$epsilon=rep(1,RC$N)
    if(any('forcepoint'== observedData$Quality)){
        RC$epsilon[forceindex]=1/RC$N
    }

    RC$Z=cbind(t(rep(0,2)),t(rep(1,RC$n)))

    RC$m1=matrix(0,nrow=2,ncol=RC$n)
    RC$m2=matrix(0,nrow=RC$n,ncol=2)
    theta.init=rep(0,9)

    Dens = function(th) {-Densevalm22(th,RC)$p}
    Densmin=optim(par=theta.init,Dens,method="L-BFGS-B",hessian=TRUE)
    t_m =Densmin$par
    H=Densmin$hessian
    phi_b=t_m[3]
    sig_b2=t_m[2]
    zeta=t_m[1]
    lambda=t_m[4:9]
    l=log(RC$w_tild+exp(t_m[1]))
    varr_m=exp(RC$B%*%lambda)
    Sig_eps=diag(as.numeric(rbind(varr_m,0)))
    R_Beta=(1+sqrt(5)*RC$dist/exp(phi_b)+5*RC$dist^2/(3*exp(phi_b)^2))*exp(-sqrt(5)*RC$dist/exp(phi_b))+diag(1,RC$n,RC$n)*RC$nugget
    Sig_x=rbind(cbind(RC$Sig_ab,matrix(0,nrow=2,ncol=RC$n)),cbind(matrix(0,nrow=RC$n,ncol=2),exp(sig_b2)*R_Beta))

    X=rbind(cbind(matrix(1,dim(l)),l,diag(as.numeric(l))%*%RC$A),RC$Z)


    L=t(chol(as.matrix(X%*%Sig_x%*%t(X)+Sig_eps)))

    w=solve(L,(-RC$y+X%*%RC$mu_x))
    mu=RC$mu_x-Sig_x%*%(t(X)%*%(solve(t(L),w)))
    LH=t(chol(H))/0.8

    cl <- makeCluster(4)
    registerDoParallel(cl)
    Wmax=as.numeric(Wmax)
    if(is.na(Wmax)){
        Wmax=ceiling(max(RC$w)*10)/10
    }
    Wmin=as.numeric(Wmin)
    if(is.na(Wmin)){
        Wmin=min(RC$w)-exp(t_m[1])
    }
    #make Wmin and Wmax divisable by 10 up, both in order to make rctafla and so l_m is defined
    Wmax=ceiling(Wmax*10)/10
    Wmin=ceiling(Wmin*10)/10
    WFill=W_unobserved(c(RC$O),min=Wmin,max=Wmax)
    RC$W_u=WFill$W_u
    RC$W_u_tild=WFill$W_u_tild
    Bsiminput=t(RC$W_u_tild)/RC$W_u_tild[length(RC$W_u_tild)]
    Bsiminput[is.na(Bsiminput)]=0
    RC$Bsim=B_splines(Bsiminput)

    MCMC <- foreach(i=1:4,.combine=cbind,.export=c("Densevalm22","predict_u")) %dopar% {
        ypo_obs=matrix(0,nrow=RC$N,ncol=Nit)
        param=matrix(0,nrow=9+RC$n+2,ncol=Nit)
        t_old=as.matrix(t_m)
        Dens<-Densevalm22(t_old,RC)
        p_old=Dens$p
        ypo_old=Dens$ypo
        x_old=Dens$x

        for(j in 1:Nit){
            t_new=t_old+solve(t(LH),rnorm(9,0,1))
            Densnew<-Densevalm22(t_new,RC)
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
    names(TableOfData)=c("Date","Time","Quality","W","Q", "Q fit","Lower", "Upper","Q diff")
    TableOfData=TableOfData[with(TableOfData,order(Date)),]

    xout=seq(Wmin,-0.01+Wmax,by=0.01)

    fitInterpolation=approx(completePrediction$W,completePrediction$fit,xout=xout)
    FitTable=t(as.data.frame(split(x=fitInterpolation$y, f=ceiling(seq_along(fitInterpolation$y)/10))))
    colnames(FitTable)=0:9
    FitTable=round(exp(FitTable),3)
    Stage=seq(min(fitInterpolation$x),max(fitInterpolation$x),by=0.1)*100
    FitTable=as.data.frame(cbind(Stage,FitTable))
    names(FitTable)[1]="Stage (cm)"

    lowerInterpolationation=approx(completePrediction$W,completePrediction$lower,xout=xout)
    LowerTable=t(as.data.frame(split(x=lowerInterpolationation$y, f=ceiling(seq_along(lowerInterpolation$y)/10))))
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

    return(list("observedData"=observedData,"betaData"=betaData,"completePrediction"=completePrediction,"observedPrediction"=observedPrediction,"TableOfData"=TableOfData,
                "FitTable"=FitTable,"LowerTable"=LowerTable,"UpperTable"=UpperTable,"plotTable"=plotTable))
}
