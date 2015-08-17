#'Calculates the rating curve using model1
#'
#'This function takes in the cleaned data from the \code{\link{clean}} function
#'and calculates the rating curve using model 1.
#'@param clean Is the output from the \code{\link{clean}} function. That is a list containg a matrix and a data frame.
#'The matrix containing the stage(in meters) and flow values and the data frame containing date, time, quality, stage and flow values.
#'@param country The input is a string with the name of the country you want to use the prior parameters from, e.g. country='Iceland'.
#'@param Wmin The input is either the default empty string or an integer for the lowest stage level. If you set as.numeric("") you get Na and then Wmin will be
#'automatically set to c_hat.
#'@param Wmax The input is either the default empty string or an integer for the highest stage level. If you set as.numeric("") you get Na and then Wmax will be
#'automatically set to the maximum stage of the data.
#'@return If all the parameters are used as described the output will be a list containing the integers varappr and c_hat, the matrix mu and
#'the data frames qvdata, simdata, realdata , tafla, fitrctafla, lowerrctafla, upperrctafla, plottafla.
model1BH <- function(clean,country="Iceland",Wmin="",Wmax=""){
    require(RCmodels)
    list2env(clean,environment())
    RC=priors(country)

    RC$y=as.matrix(log(wq[,2]))
    RC$w=wq[,1]
    RC$w_tild=RC$w-min(RC$w)
    RC$n=length(RC$y)
    epsilon=0.00001

    forceindex=which('forcepoint'== qvdata$Quality)
    forcepoint=wq[forceindex,]

    Dens <- function(th){ Densevalm11(th,RC)$pmin}
    Densmin=optim(par=c(0,0),Dens,hessian=TRUE)
    t_m=as.matrix(Densmin$par)
    H=Densmin$hessian

    sigma=rep(exp(t_m[2,]),nrow(wq))

    if(any('forcepoint'== qvdata$Quality)){
    sigma[forceindex]=epsilon
    }

    l_m=as.matrix(log(RC$w_tild+exp(t_m[1,])))

    X_m=cbind(rep(1,nrow(l_m)),l_m)

    L=t(chol(RC$Sig_xinv+t(X_m)%*%(X_m/sigma)))

    mu=solve(t(L),(solve(L,(RC$Sinvmu+t(X_m)%*%(RC$y/sigma)))))

    v_temp=X_m%*%solve(RC$Sig_xinv+t(X_m)%*%(X_m/sigma))%*%t(X_m)

    varappr=as.matrix(diag(v_temp)+sigma)

    constvar=mean(as.matrix(diag(v_temp)+exp(t_m[2,])))

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
    Wmin=as.numeric(Wmin)
    if(is.na(Wmin)){
        Wmin=c_hat
    }
    Wmax=ceiling(Wmax*10)/10
    Wmin=ceiling(Wmin*10)/10

    W_grid=sort(c(seq(Wmin,Wmax,length.out=1000),forcepoint[1]))
    forcesimindex=which(forcepoint[1]==W_grid)

    simdata=data.frame(W=W_grid)
    simdata$l_m = log(simdata$W-c_hat)
    simdata$fit=mu[1,]+mu[2,]*simdata$l_m

    simvar=rep(constvar,nrow(simdata))
    simvar[forcesimindex]=epsilon

    simdata$upper=simdata$fit+qnorm(0.975,0,sqrt(simvar))
    simdata$lower=simdata$fit+qnorm(0.025,0,sqrt(simvar))
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

    #-0.01 inorder for rctafla to be of the right length
    xout=seq(Wmin,-0.01+Wmax,by=0.01)

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

    return(list("varappr"=constvar,"qvdata"=qvdata,"simdata"=simdata,"realdata"=realdata,
                "tafla"=tafla,"mu"=mu,"c_hat"=c_hat,"fitrctafla"=fitrctafla,"lowerrctafla"=lowerrctafla,"upperrctafla"=upperrctafla,"plottafla"=plottafla))
}
