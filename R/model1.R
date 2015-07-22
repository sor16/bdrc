model1BH <- function(clean,country="Iceland",Wmax=NA){
    require(RCmodels)
    list2env(clean,environment())
    RC=priors(country)

    RC$y=as.matrix(log(wq[,2]))
    RC$w=wq[,1]
    RC$w_tild=RC$w-min(RC$w)
    RC$n=length(RC$y)


    Dens <- function(th){ Densevalm11(th,RC)$pmin}
    Densmin=optim(par=c(0,0),Dens,hessian=TRUE)
    t_m=as.matrix(Densmin$par)
    H=Densmin$hessian


    l_m=as.matrix(log(RC$w_tild+exp(t_m[1,])))

    X_m=cbind(matrix(1,nrow(l_m),ncol(l_m)),l_m)

    L=t(chol(RC$Sig_xinv+t(X_m)%*%X_m/exp(t_m[2,])))

    mu=solve(t(L),(solve(L,(RC$Sinvmu+t(X_m)%*%RC$y/exp(t_m[2,])))))

    v_temp=X_m%*%solve(RC$Sig_xinv+t(X_m)%*%X_m/exp(t_m[2,]))%*%t(X_m)

    varappr=mean(as.matrix(diag(v_temp)+exp(t_m[2,])))

    RC$fit=X_m%*%mu

    RC$confinterval= cbind(X_m%*%mu+qnorm(0.025,0,sqrt(varappr)),X_m%*%mu+qnorm(0.975,0,sqrt(varappr)))

    realdata=data.frame(W=RC$w,Q=RC$y)
    realdata$l_m=l_m
    realdata$fit=RC$fit
    realdata$upper=RC$confinterval[,2]
    realdata$lower=RC$confinterval[,1]
    c_hat=min(realdata$W)-exp(t_m[1,])
    Wmax=as.numeric(Wmax)
    if(is.na(Wmax)){
        Wmax=max(RC$w)
    }
    simdata=data.frame(W=seq(ceiling(c_hat*10)/10,ceiling(Wmax*10)/10,length.out=1000))
    simdata$l_m = log(simdata$W-c_hat)
    simdata$fit=mu[1,]+mu[2,]*simdata$l_m
    simdata$upper=simdata$fit+qnorm(0.975,0,sqrt(varappr))
    simdata$lower=simdata$fit+qnorm(0.025,0,sqrt(varappr))
    realdata$residraun=(exp(realdata$Q)-exp(realdata$fit))
    realdata$residupper=exp(realdata$upper)-exp(realdata$fit)
    realdata$residlower=exp(realdata$lower)-exp(realdata$fit)
    realdata$residlog=(realdata$Q-realdata$fit)/sqrt(exp(t_m[2,]))

    tafla=qvdata
    tafla$W=tafla$W
    tafla$Q=round(tafla$Q,1)
    tafla$Qfit=as.numeric(round(exp(realdata$fit),3))
    tafla$lower=round(exp(realdata$lower),3)
    tafla$upper=round(exp(realdata$upper),3)
    tafla$diffQ=tafla$Q-tafla$Qfit
    names(tafla)=c("Date","Time","Quality","W","Q", "Q fit","Lower", "Upper","Q diff")
    tafla=tafla[with(tafla,order(Date)),]

    xout=seq(ceiling(c_hat*10)/10,-0.01+ceiling(Wmax*10)/10,by=0.01)

    fitinterpol=approx(simdata$W,simdata$fit,xout=xout)
    fitrctafla=t(as.data.frame(split(x=fitinterpol$y, f=ceiling(seq_along(fitinterpol$y)/10))))
    colnames(fitrctafla)=0:9
    fitrctafla=round(exp(fitrctafla),3)
    Stage=seq(min(fitinterpol$x),max(fitinterpol$x),by=0.1)*100
    fitrctafla=as.data.frame(cbind(Stage,fitrctafla))

    lowerinterpol=approx(simdata$W,simdata$lower,xout=xout)
    lowerrctafla=t(as.data.frame(split(x=lowerinterpol$y, f=ceiling(seq_along(lowerinterpol$y)/10))))
    colnames(lowerrctafla)=0:9
    lowerrctafla=round(exp(lowerrctafla),3)
    Stage=seq(min(lowerinterpol$x),max(lowerinterpol$x),by=0.1)*100
    lowerrctafla=as.data.frame(cbind(Stage,lowerrctafla))

    upperinterpol=approx(simdata$W,simdata$upper,xout=xout)
    upperrctafla=t(as.data.frame(split(x=upperinterpol$y, f=ceiling(seq_along(upperinterpol$y)/10))))
    colnames(upperrctafla)=0:9
    upperrctafla=round(exp(upperrctafla),3)
    Stage=seq(min(upperinterpol$x),max(upperinterpol$x),by=0.1)*100
    upperrctafla=as.data.frame(cbind(Stage,upperrctafla))

    plottafla=as.data.frame(cbind(lowerinterpol$y,fitinterpol$y,upperinterpol$y))
    plottafla=exp(plottafla)
    names(plottafla)=c("Lower","Fit","Upper")
    plottafla$W=xout

    return(list("varappr"=varappr,"qvdata"=qvdata,"simdata"=simdata,"realdata"=realdata,
                "tafla"=tafla,"mu"=mu,"c_hat"=c_hat,"fitrctafla"=fitrctafla,"lowerrctafla"=lowerrctafla,"upperrctafla"=upperrctafla,"plottafla"=plottafla))
}
