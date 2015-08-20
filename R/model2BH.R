#'Calculates the rating curve using model2
#'
#'This function takes in the cleaned data from the \code{\link{clean}} function
#'and calculates the rating curve using model 2.
#'@param clean Is the output from the \code{\link{clean}} function. That is a list containg a matrix and a data frame.
#'The matrix containing the stage(in meters) and flow values and the data frame containing date, time, quality, stage and flow values.
#'@param country The input is a string with the name of the country you want to use the prior parameters from, e.g. country='Iceland'.
#'@param Wmin The input is either the default empty string or an integer for the lowest stage level. If you set as.numeric("") you get Na and then Wmin will be
#'automatically set to c_hat.
#'@param Wmax The input is either the default empty string or an integer for the highest stage level. If you set as.numeric("") you get Na and then Wmax will be
#'automatically set to the maximum stage of the data.
#'@return If all the parameters are used as described the output will be a list containing
#'the data frames qvdata, betadata, ypodata, realdata , tafla, fitrctafla, lowerrctafla, upperrctafla, plottafla.
#'@references Birgir Hrafnkelsson, Helgi Sigurdarson and Sigurdur M. Gardarson (2015) \emph{Bayesian Generalized Rating Curves}
#'@seealso \code{\link{clean}}
model2BH <- function(clean,country="Iceland",Wmin="",Wmax=""){
    require(doParallel)
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

    forceindex=which('forcepoint'== qvdata$Quality)
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
    if(any('forcepoint'== qvdata$Quality)){
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
    MCMC[is.na(MCMC)]=0
    betasamples=apply(MCMC[(RC$N+length(RC$W_u)+1):nrow(MCMC),],2,FUN=function(x){x[2]+x[3:length(x)]})
    yposamples=MCMC[1:(RC$N+length(RC$W_u)),]
    ypodata=as.data.frame(t(apply(yposamples,1,quantile, probs = c(0.025,0.5, 0.975),na.rm=T)))
    names(ypodata)=c("lower","fit","upper")
    ypodata$W=c(RC$w,RC$W_u)
    ypodata$l_m=c(l,log(RC$W_u-min(RC$O)+exp(t_m[1])))
    realdata=ypodata[1:RC$N,]
    ypodata=ypodata[with(ypodata,order(W)),]
    betadata=as.data.frame(t(apply(betasamples,1,quantile, probs = c(0.025,0.5, 0.975),na.rm=T)))
    names(betadata)=c("lower","fit","upper")
    betadata$W=c(RC$O,RC$W_u)
    betadata=betadata[with(betadata,order(W)),]
    realdata$Q=RC$y[1:RC$N,]
    realdata$residraun=(exp(realdata$Q)-exp(realdata$fit))
    realdata$residupper=exp(realdata$upper)-exp(realdata$fit)
    realdata$residlower=exp(realdata$lower)-exp(realdata$fit)
    realdata$residlog=(realdata$Q-realdata$fit)/sqrt(varr_m)

    tafla=qvdata
    tafla$Q=round(tafla$Q,1)
    tafla$Qfit=round(exp(realdata$fit),3)
    tafla$lower=round(exp(realdata$lower),3)
    tafla$upper=round(exp(realdata$upper),3)
    tafla$diffQ=tafla$Q-tafla$Qfit
    names(tafla)=c("Date","Time","Quality","W","Q", "Q fit","Lower", "Upper","Q diff")
    tafla=tafla[with(tafla,order(Date)),]

    xout=seq(Wmin,-0.01+Wmax,by=0.01)

    fitinterpol=approx(ypodata$W,ypodata$fit,xout=xout)
    fitrctafla=t(as.data.frame(split(x=fitinterpol$y, f=ceiling(seq_along(fitinterpol$y)/10))))
    colnames(fitrctafla)=0:9
    fitrctafla=round(exp(fitrctafla),3)
    Stage=seq(min(fitinterpol$x),max(fitinterpol$x),by=0.1)*100
    fitrctafla=as.data.frame(cbind(Stage,fitrctafla))

    lowerinterpol=approx(ypodata$W,ypodata$lower,xout=xout)
    lowerrctafla=t(as.data.frame(split(x=lowerinterpol$y, f=ceiling(seq_along(lowerinterpol$y)/10))))
    colnames(lowerrctafla)=0:9
    lowerrctafla=round(exp(lowerrctafla),3)
    Stage=seq(min(lowerinterpol$x),max(lowerinterpol$x),by=0.1)*100
    lowerrctafla=as.data.frame(cbind(Stage,lowerrctafla))

    upperinterpol=approx(ypodata$W,ypodata$upper,xout=xout)
    upperrctafla=t(as.data.frame(split(x=upperinterpol$y, f=ceiling(seq_along(upperinterpol$y)/10))))
    colnames(upperrctafla)=0:9
    upperrctafla=round(exp(upperrctafla),3)
    Stage=seq(min(upperinterpol$x),max(upperinterpol$x),by=0.1)*100
    upperrctafla=as.data.frame(cbind(Stage,upperrctafla))

    plottafla=as.data.frame(cbind(lowerinterpol$y,fitinterpol$y,upperinterpol$y))
    plottafla=exp(plottafla)
    names(plottafla)=c("Lower","Fit","Upper")
    plottafla$W=xout

    return(list("qvdata"=qvdata,"betadata"=betadata,"ypodata"=ypodata,"realdata"=realdata,"tafla"=tafla,
                "fitrctafla"=fitrctafla,"lowerrctafla"=lowerrctafla,"upperrctafla"=upperrctafla,"plottafla"=plottafla))
}
