#'Plots the rating curve for model 1
#'
#'This function takes in the filename of the data that is used to calculate the ratingcurve It uses the \code{\link{clean}} function
#'and the \code{\link{model1BH}} function to calculate the rating curve and then plot it.
#'@param filename Input is a string with the name of the txt file containing the stage and flow data. See also input in \code{\link{clean}}.
#'@param eps TRUE or FALSE whether or not to save the plots as eps files or not.
#'@param png TRUE or FALSE whether or not to save the plots as png files or not.
#'@param realscale Logical constant, whether or not to plot the real scale image.
#'@param logscale Logical constant, whether or not to plot the log scale image.
#'@param logresiduals Logical constant, whether or not to plot the residuals for the log scale.
#'@param realresiduals Logical constant, whether or not to plot the residuals for the real scale.
#'@return The output is a list which contains the ggplots the user has chosen to plot.
#'@references Birgir Hrafnkelsson, Helgi Sigurdarson and Sigurdur M. Gardarson (2015) \emph{Bayesian Generalized Rating Curves}
#'@seealso \code{\link{clean}} , \code{\link{model1BH}}
plotmodel1 <- function(filename,eps=F,png=F,realscale=T,logscale=T,logresiduals=T,realresiduals=T){
    require(RCmodels)
    require(ggplot2)
    if (is.null(filename)){
        return(NULL)
    }

    data <- model1BH(clean(filename))
    simdata=data$simdata
    realdata=data$realdata

    plotlist=list()

    if(realscale==TRUE){
        rcreal=ggplot(simdata)+theme_bw()+geom_point(data=realdata,aes(exp(Q),W))+geom_path(aes(exp(fit),W))+
               geom_path(aes(exp(lower),W),linetype="dashed")+geom_path(aes(exp(upper),W),linetype="dashed")+
               ggtitle(paste("Rating curve for",gsub("\\.[^.]*$","",filename)))+ylab("W  [m]")+xlab(expression(paste("Q  [",m^3,'/s]',sep='')))+
               theme(plot.title = element_text(vjust=2))


        plotlist$rcreal=rcreal
    }
    if(logscale==TRUE){
        rclog=ggplot(realdata)+geom_path(mapping=aes(fit,l_m))+theme_bw()+geom_point(mapping=aes(Q,l_m))+geom_path(mapping=aes(lower,l_m),linetype="dashed")+
              geom_path(mapping=aes(upper,l_m),linetype="dashed")+ggtitle(paste("Rating curve for",gsub("\\.[^.]*$","",filename),"(log scale)"))+
              ylab(expression(log(W-hat(c))))+xlab("log(Q)")+theme(plot.title = element_text(vjust=2))

        plotlist$rclog=rclog
    }
    if(logresiduals==TRUE){
        max=max(abs(realdata$residlog))
        if(max>4){
            ylim=c(-(max+0.2),max+0.2)
        }else{
            ylim=c(-4,4)
        }
        rcrealresid=ggplot(realdata)+geom_point(aes(W,residraun),color="red")+theme_bw()+geom_abline(intercept = 0, slope = 0)+
            geom_path(aes(W,residupper),linetype="dashed")+geom_path(aes(W,residlower),linetype="dashed")+ylab(expression(paste("Q - ",hat(Q) ,"  [",m^3,'/s]',sep='')))+
            ggtitle("Residual plot")+xlab("W  [cm]")+theme(plot.title = element_text(vjust=2))+ylim(ylim)

        plotlist$rcrealresid=rcrealresid
    }
    if(realresiduals==TRUE){
        rclogresid=ggplot(realdata)+geom_point(aes(l_m,residlog),color="red")+theme_bw()+geom_abline(intercept = 0, slope = 0)+
            geom_abline(intercept = 2, slope = 0,linetype="dashed")+geom_abline(intercept = -2, slope = 0,linetype="dashed")+
            ylab(expression(epsilon[i]))+ggtitle("Residual plot (log scale)")+xlab(expression(log(W-hat(c))))+
            theme(plot.title = element_text(vjust=2))


        plotlist$rclogresid=rclogresid
    }
    if(png==TRUE){
        if(length(plotlist)!=0){
            for(i in 1:length(plotlist)){
            png(file = paste(gsub("\\.[^.]*$","",filename),'_',i,'_','model1','.png', sep=""))
            print(plotlist[i])
            dev.off()
            }
        }
    }
    if(eps==TRUE){
        if(length(plotlist)!=0){
            for(i in 1:length(plotlist)){
            postscript(file = paste(gsub("\\.[^.]*$","",filename),'_',i,'_','model1','.eps', sep=""), width = 10, height = 7.5)
            print(plotlist[i])
            dev.off()
            }
        }
    }

    return(plotlist)
}
