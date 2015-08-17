#'Plots the rating curve for model 2
#'
#'This function takes in the filename of the data that is used to calculate the ratingcurve It uses the \code{\link{clean}} function
#'and the \code{\link{model2BH}} function to calculate the rating curve and then plot it.
#'@param filename Input is a string with the name of the txt file containing the stage and flow data. See also input in \code{\link{clean}}.
#'@param eps TRUE or FALSE whether or not to save the plots as eps files or not.
#'@param png TRUE or FALSE whether or not to save the plots as png files or not.
#'@param realscale Logical constant, whether or not to plot the real scale image.
#'@param logscale Logical constant, whether or not to plot the log scale image.
#'@param logresiduals Logical constant, whether or not to plot the residuals for the log scale.
#'@param realresiduals Logical constant, whether or not to plot the residuals for the real scale.
#'@param beta Logical constant, whether or not to plot the beta image.
#'@return The output is a list which contains the ggplots the user has chosen to plot.
plotmodel2 <- function(filename,eps=F,png=F,realscale=T,logscale=T,logresiduals=T,realresiduals=T,beta=T){

    library(ggplot2)
    if (is.null(filename)){
        return(NULL)
    }

    require(RCmodels)
    list2env(clean,environment())

    data <- model2BH(clean(filename))
    betadata=data$betadata
    ypodata=data$ypodata
    simdata=data$simdata
    realdata=data$realdata

    plotlist=list()

    if(realscale==TRUE){
        rcreal=ggplot(ypodata)+theme_bw()+geom_point(data=realdata,aes(exp(Q),W))+geom_line(aes(exp(fit),W))+
            geom_line(aes(exp(lower),W),linetype="dashed")+geom_line(aes(exp(upper),W),linetype="dashed")+
            ggtitle(paste("Rating curve for",gsub("\\.[^.]*$","",filename)))+ylab("W [m]")+xlab(expression(paste("Q  [",m^3,'/s]',sep='')))+
            theme(plot.title = element_text(vjust=2))

        plotlist$rcreal=rcreal
    }
    if(logscale==TRUE){
        rclog=ggplot(realdata)+geom_line(mapping=aes(fit,l_m))+theme_bw()+geom_point(mapping=aes(Q,l_m))+geom_line(mapping=aes(lower,l_m),linetype="dashed")+
            geom_line(mapping=aes(upper,l_m),linetype="dashed")+ggtitle(paste("Rating curve for",gsub("\\.[^.]*$","",filename),"(log scale)"))+
            ylab(expression(log(W-hat(c))))+xlab("log(Q)")+theme(plot.title = element_text(vjust=2))

        plotlist$rclog=rclog
    }
    if(logresiduals==TRUE){
        rcrealresid=ggplot(realdata)+geom_point(aes(W,residraun),color="red")+theme_bw()+geom_abline(intercept = 0, slope = 0)+
            geom_line(aes(W,residupper),linetype="dashed")+geom_line(aes(W,residlower),linetype="dashed")+ylab(expression(paste("Q - ",hat(Q) ,"  [",m^3,'/s]',sep='')))+
            ggtitle("Residual plot")+xlab("W  [cm]")+theme(plot.title = element_text(vjust=2))

        plotlist$rcrealresid=rcrealresid
    }
    if(realresiduals==TRUE){
        rclogresid=ggplot(realdata)+geom_point(aes(l_m,residlog),color="red")+theme_bw()+geom_abline(intercept = 0, slope = 0)+
            geom_abline(intercept = 2, slope = 0,linetype="dashed")+geom_abline(intercept = -2, slope = 0,linetype="dashed")+
            ylab(expression(epsilon[i]))+ggtitle("Residual plot (log scale)")+xlab(expression(log(W-hat(c))))+
            theme(plot.title = element_text(vjust=2))


        plotlist$rclogresid=rclogresid
    }
    if(beta==TRUE){
        smoothbeta=ggplot(data=betadata)+geom_line(aes(W,fit))+
                   geom_line(aes(W,lower),linetype="dashed")+geom_line(aes(W,upper),linetype="dashed")+
                   ylab(expression(b+beta(W)))+ggtitle("b parameter as a function of stage W")+xlab("W [m]")+
                   theme(plot.title = element_text(vjust=2))+theme_bw()
        plotlist$smoothbeta=smoothbeta
    }

    if(png==TRUE){
        if(length(plotlist)!=0){
            for(i in 1:length(plotlist)){
                png(file = paste(gsub("\\.[^.]*$","",filename),'_',i,'_','model2','.png', sep=""))
                print(plotlist[i])
                dev.off()
            }
        }
    }
    if(eps==TRUE){
        if(length(plotlist)!=0){
            for(i in 1:length(plotlist)){
                postscript(file = paste(gsub("\\.[^.]*$","",filename),'_',i,'_','model2','.eps', sep=""), width = 10, height = 7.5)
                print(plotlist[i])
                dev.off()
            }
        }
    }

    return(plotlist)
}
