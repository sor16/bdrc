#'Calculates the rating curve using model1
#'
#'This function takes in the cleaned data from the \code{\link{clean}} function
#'and calculates the rating curve using model 1.
#'@param clean List that is the output of \code{\link{clean}} function.
#'@param country A string of t the country the prior parameters should be defined, default value is Iceland.
#'@param Wmin positive numeric value for the lowest stage the user wants to calculate a rating curve. If input is an empty string (default) Wmin will
#'automatically be set to c_hat.
#'@param Wmax positive numeric value for the highest stage the user wants to calculate a rating curve. If input is an empty string (default) Wmax will
#'automatically be set to the maximum stage of the data.
#'@return List containing information on the calculated rating curve, parameters varappr and c_hat, the matrix mu and
#'the data frames observedData, completePrediction, observedPrediction , TableOfData, FitTable, LowerTable, UpperTable, plotTable.
#'@references Birgir Hrafnkelsson, Helgi Sigurdarson and Sigurdur M. Gardarson (2015) \emph{Bayesian Generalized Rating Curves}
#'@seealso \code{\link{clean}}
model1BH <- function(clean,country="Iceland",Wmin="",Wmax=""){
    require(RCmodels)
    list2env(clean,environment())
    RC=priors(country)

    RC$y=as.matrix(log(wq[,2]))
    RC$w=wq[,1]
    RC$w_tild=RC$w-min(RC$w)
    RC$n=length(RC$y)
    epsilon=0.00001

    forceIndex=which('forcepoint'== observedData$Quality)
    forcepoint=wq[forceIndex,]

    Dens <- function(th){ Densevalm11(th,RC)$pmin}
    Densmin=optim(par=c(0,0),Dens,hessian=TRUE)
    t_m=as.matrix(Densmin$par)
    H=Densmin$hessian
    sigma=rep(exp(t_m[2,]),nrow(wq))
    if(any('forcepoint'== observedData$Quality)){
    sigma[forceIndex]=epsilon
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

    observedPrediction=data.frame(W=RC$w,Q=RC$y)
    observedPrediction$l_m=l_m
    observedPrediction$fit=RC$fit
    observedPrediction$upper=RC$confinterval[,2]
    observedPrediction$lower=RC$confinterval[,1]
    c_hat=min(observedPrediction$W)-exp(t_m[1,])
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
    indexOfForcepoint=which(forcepoint[1]==W_grid)

    completePrediction=data.frame(W=W_grid)
    completePrediction$l_m = log(completePrediction$W-c_hat)
    completePrediction$fit=mu[1,]+mu[2,]*completePrediction$l_m

    completeVariance=rep(constvar,nrow(completePrediction))
    completeVariance[indexOfForcepoint]=epsilon

    completePrediction$upper=completePrediction$fit+qnorm(0.975,0,sqrt(completeVariance))
    completePrediction$lower=completePrediction$fit+qnorm(0.025,0,sqrt(completeVariance))
    observedPrediction$residuals=(exp(observedPrediction$Q)-exp(observedPrediction$fit))
    observedPrediction$residupper=exp(observedPrediction$upper)-exp(observedPrediction$fit)
    observedPrediction$residlower=exp(observedPrediction$lower)-exp(observedPrediction$fit)
    observedPrediction$standardResiduals=(observedPrediction$Q-observedPrediction$fit)/sqrt(exp(t_m[2,]))

    TableOfData=observedData
    TableOfData$W=TableOfData$W
    TableOfData$Q=round(TableOfData$Q,1)
    TableOfData$Qfit=as.numeric(round(exp(observedPrediction$fit),3))
    TableOfData$lower=round(exp(observedPrediction$lower),3)
    TableOfData$upper=round(exp(observedPrediction$upper),3)
    TableOfData$diffQ=TableOfData$Q-TableOfData$Qfit
    TableOfData$Qpercentage=round(100*TableOfData$diffQ/TableOfData$Q,1)
    names(TableOfData)=c("Date","Time","Quality","W","Q", "Q fit","Lower", "Upper","Q diff","Q%")
    TableOfData=TableOfData[with(TableOfData,order(Date)),]

    #-0.01 in order to fix the length of the interpolation tables
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

    return(list("varappr"=constvar,"observedData"=observedData,"completePrediction"=completePrediction,"observedPrediction"=observedPrediction,
                "TableOfData"=TableOfData,"mu"=mu,"c_hat"=c_hat,"FitTable"=FitTable,"LowerTable"=LowerTable,"UpperTable"=UpperTable,"plotTable"=plotTable))
}
