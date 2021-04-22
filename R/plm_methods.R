print_fun <- function(x){
    cat("\nCall:\n",
        paste(deparse(x$formula), sep = "\n", collapse = "\n"), "\n\n", sep = "")
}

summary_fun <- function(x){
    param_summary <- x$param_summary[,c('lower','median','upper')]
    names(param_summary) <- paste0(names(param_summary),c('-2.5%','-50%','-97.5%'))
    cat("\nFormula: \n",
        paste(deparse(x$formula), sep = "\n", collapse = "\n"))
    cat("\nLatent parameters:\n")
    print(param_summary[1:2,],row.names = T,digits=3,right=F)
    cat("\nHyperparameters:\n")
    print(param_summary[3:nrow(param_summary),],row.names = T,digits=3,right=F)
    cat("\nDIC:",x$DIC)
}

#' Custom bdrc theme
#'
#' @return returns a theme object for the package
#' @export
#' @importFrom ggplot2 %+replace% theme_classic theme element_text element_blank
theme_bdrc <- function(){
    theme_classic() %+replace%
        theme( #text = element_text(family="Times", face="plain"),
               strip.background = element_blank(),
               strip.text.x = element_text(size = 16),
               axis.title.x = element_text(size=16),
               axis.title.y = element_text(size=16,angle=90),
               axis.text.x = element_text(size=12),
               axis.text.y = element_text(size=12),
               legend.text=element_text(size=12),
               legend.title=element_text(size=16),
               plot.title=element_text(size=18))
}

#' Plot bdrc model objects
#'
#' Visualize results from model ojbects in bdrc, plm0, plm, gplm0,gplm
#' @param x an object of class "plm0","plm","gplm0" or "gplm".
#' @param type a character denoting what type of plot should be drawn. Defaults to "rating_curve". Possible types are
#'                    \itemize{
#'                       \item{"rating_curve"}{ to plot the rating curve.}
#'                       \item{"rating_curve_mean"}{ to plot the posterior mean of the rating curve.}
#'                       \item{"f"}{ to plot the power-law exponent.}
#'                       \item{"beta"}{ to plot the random effect in the power-law exponent.}
#'                       \item{"sigma_eps"}{ to plot the standard deviation on the data level.}
#'                       \item{"residuals"}{ to plot the log residuals.}
#'                       \item{"trace"}{ to plot trace plots of parameters given in param.}
#'                       \item{"histogram"}{ to plot histograms of parameters given in param.}
#'                    }
#' @param param a character vector with the parameters to plot. Defaults to NULL and is only used if type is "trace" or "histogram". Allowed values are the parameters given in the model summary of x as well as "hyperparameters" or "latent_parameters" for specific groups of parameters.
#' @param transformed a logical value indicating whether the quantity should be plotted on a transformed scale used during the Bayesian inference. Defaults to FALSE.
#' @param title a character denoting the title of the plot. Defaults to NULL, i.e. no title.
#' @return returns an object of class ggplot2 or Grob object.
#' @importFrom ggplot2 ggplot aes geom_point geom_path geom_histogram geom_abline facet_wrap scale_color_manual scale_x_continuous scale_y_continuous label_parsed ggtitle xlab ylab
#' @importFrom rlang .data
#' @importFrom stats median
plot_fun <- function(x,type='rating_curve',param=NULL,transformed=F,title=NULL){
    legal_types <- c('rating_curve','rating_curve_mean','f','beta','sigma_eps','residuals','trace','histogram')
    types_with_param <- c('trace','histogram')
    if(!(type %in% legal_types)){
        stop(cat(paste('Type argument not recognized. Possible types are:\n -',paste(legal_types,collapse='\n - '))))
    }else if(type %in% types_with_param & is.null(param)){
        stop('If type histogram or trace, param must be non-null, should be a character vector of the names parameters to visualize.')
    }else if(type=='trace'){
        plot_dat <- gather_draws(x,param,transformed=transformed)
        if('h' %in% names(plot_dat)){
            stop('Plots of type "trace" can only be of stage-independent parameters')
        }
        params <- unique(plot_dat$name)
        if(length(params)>1){
            param_levels <- get_parameter_levels(params)
            plot_dat$name_expr <- factor(plot_dat$name,levels=param_levels,labels=sapply(param_levels,get_param_expression))
            plot_dat$chain <- factor(as.character(plot_dat$chain),levels=1:max(plot_dat$chain))
            p <- ggplot(plot_dat,aes(x=.data$iter,y=.data$value,col=.data$chain)) +
                geom_path() +
                facet_wrap(~name_expr,scales='free',labeller = label_parsed) +
                scale_color_manual(values=c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF"),
                                   name='Chain number') +
                xlab('Iteration') +
                ylab('') +
                theme_bdrc()
        }else{
            param_expr <- get_param_expression(params)
            plot_dat$chain_name <- paste0('Chain nr ',plot_dat$chain)
            p <- ggplot(plot_dat,aes(x=.data$iter,y=.data$value)) +
                geom_path(col="#0072B5FF") +
                facet_wrap(~chain_name,scales='free') +
                xlab('Iteration') +
                ylab(parse(text=param_expr)) +
                theme_bdrc()
        }
    }else if(type=='histogram'){
        plot_dat <- gather_draws(x,param,transformed=transformed)
        if('h' %in% names(plot_dat)){
            stop('Plots of type "histogram" can only be of stage-independent parameters')
        }
        params <- unique(plot_dat$name)
        param_levels <- get_parameter_levels(params)
        plot_dat$name_expr <- factor(plot_dat$name,levels=param_levels,labels=sapply(param_levels,get_param_expression))
        plot_dat$chain <- factor(as.character(plot_dat$chain),levels=1:max(plot_dat$chain))
        p <- ggplot(plot_dat,aes(x=.data$value)) +
            geom_histogram(bins=50,fill="#0072B5FF") +
            facet_wrap(~name_expr,scales='free',labeller = label_parsed) +
            xlab('') +
            ylab('') +
            theme_bdrc()
    }else if(type=='rating_curve' | type=='rating_curve_mean'){
        if(transformed){
            #to generate label - latex2exp::TeX('$\\log(\\textit{h-\\hat{c}})$','character')
            x_lab <- "paste('','',log,,,,'(','',italic(paste('h-',hat(paste('c')))),')','','')"
            #to generate label - latex2exp::TeX('$\\log(\\textit{Q})$','character')
            y_lab <- "paste('','',log,,,,'(','',italic(paste('Q')),')','','')"
            c_hat <- if(is.null(x$run_info$c_param)) median(x$c_posterior) else x$run_info$c_param
            plot_dat <- merge(x[[type]],x$data,by.x='h',by.y=all.vars(x$formula)[2])
            plot_dat[,'log(h-c_hat)'] <- log(plot_dat$h-c_hat)

            # plot_dat$log_Q <- log(plot_dat$Q)
            plot_dat$log_Q <- log(plot_dat[, all.vars(x$formula)[1]])

            plot_dat$log_lower <- log(plot_dat$lower)
            plot_dat$log_median <- log(plot_dat$median)
            plot_dat$log_upper <- log(plot_dat$upper)
            p <- ggplot(data=plot_dat) +
                geom_point(aes(x=.data$`log(h-c_hat)`,y=.data$log_Q)) +
                geom_path(aes(x=.data$`log(h-c_hat)`,y=.data$log_median)) +
                geom_path(aes(x=.data$`log(h-c_hat)`,y=.data$log_lower),linetype='dashed') +
                geom_path(aes(x=.data$`log(h-c_hat)`,y=.data$log_upper),linetype='dashed') +
                scale_x_continuous(expand=c(0.01,0)) +
                xlab(parse(text=x_lab)) +
                ylab(parse(text=y_lab)) +
                theme_bdrc()
        }else{
            #to generate label - latex2exp::TeX('$\\textit{Q}\\lbrack\\textit{m^3/s}\\rbrack$','character')
            x_lab <- "paste('','',italic(paste('Q')),paste('['),italic(paste('m',phantom() ^ {paste('3')},'/s')),paste(']'),'')"
            #to generate label - latex2exp::TeX('$\\textit{h}\\lbrack\\textit{m}\\rbrack$','character')
            y_lab <- "paste('','',italic(paste('h')),paste('['),italic(paste('m')),paste(']'),'')"
            p <- ggplot(data=x[[type]]) +
                geom_point(data=x$data,aes(.data[[all.vars(x$formula)[1]]],.data[[all.vars(x$formula)[2]]])) +
                geom_path(aes(x=.data$median,y=.data$h)) +
                geom_path(aes(x=.data$lower,y=.data$h),linetype='dashed') +
                geom_path(aes(x=.data$upper,y=.data$h),linetype='dashed') +
                scale_x_continuous(expand=c(0,0)) +
                scale_y_continuous(expand=c(0.01,0)) +
                xlab(parse(text=x_lab)) +
                ylab(parse(text=y_lab)) +
                theme_bdrc()
        }
    }else if(type=='sigma_eps'){
        #to generate label - latex2exp::TeX('$\\textit{h}\\lbrack\\textit{m}\\rbrack$','character')
        x_lab <- "paste('','',italic(paste('h')),paste('['),italic(paste('m')),paste(']'),'')"
        h_in_data <- x$data[,all.vars(x$formula)[2],drop=T]
        if('sigma_eps_summary' %in% names(x)){
            #to generate label - latex2exp::TeX('$\\sigma_{\\epsilon}(\\textit{h})$','character')
            y_lab <- "paste('','',sigma,,,,phantom() [ {paste('',epsilon,,,)} ],'(','',italic(paste('h')),')','','')"
            plot_dat <- x$sigma_eps_summary[x$sigma_eps_summary$h>=min(h_in_data) & x$sigma_eps_summary$h<=max(h_in_data),]
        }else{
            #to generate label - latex2exp::TeX('$\\sigma_{\\epsilon}$','character')
            y_lab <- "paste('','',sigma,,,,phantom() [ {paste('',epsilon,,,)} ],'')"
            plot_dat <- data.frame(h=x$data[,all.vars(x$formula)[2],drop=T],
                                   lower=x$param_summary['sigma_eps','lower'],
                                   median=x$param_summary['sigma_eps','median'],
                                   upper=x$param_summary['sigma_eps','upper'])
        }
        p <- ggplot(data=plot_dat) +
            geom_path(aes(x=.data$h,y=.data$median)) +
            geom_path(aes(x=.data$h,y=.data$lower),linetype='dashed') +
            geom_path(aes(x=.data$h,y=.data$upper),linetype='dashed') +
            xlab(parse(text=x_lab)) +
            ylab(parse(text=y_lab)) +
            theme_bdrc()
    }else if(type=='beta'){
        if(!('beta_summary' %in% names(x))){
            stop('Plots of type "beta" are only for models with stage dependent power law exponent, s.a. "gplm0" and "gplm"')
        }
        #to generate label - latex2exp::TeX('$\\textit{h}\\lbrack\\textit{m}\\rbrack$','character')
        x_lab <- "paste('','',italic(paste('h')),paste('['),italic(paste('m')),paste(']'),'')"
        #to generate label - latex2exp::TeX('$\\beta(\\textit{h})$','character')
        y_lab <- "paste('','',beta,,,,'(','',italic(paste('h')),')','','')"
        h_in_data <- x$data[,all.vars(x$formula)[2],drop=T]
        p <- ggplot(data=x$beta_summary[x$beta_summary$h>=min(h_in_data) & x$beta_summary$h<=max(h_in_data),]) +
            geom_path(aes(.data$h,.data$median)) +
            geom_path(aes(.data$h,.data$lower),linetype='dashed') +
            geom_path(aes(.data$h,.data$upper),linetype='dashed') +
            xlab(parse(text=x_lab)) +
            ylab(parse(text=y_lab)) +
            theme_bdrc()
    }else if(type=='f'){
        #to generate label - latex2exp::TeX('$\\textit{h}\\lbrack\\textit{m}\\rbrack$','character')
        x_lab <- "paste('','',italic(paste('h')),paste('['),italic(paste('m')),paste(']'),'')"
        h_in_data <- x$data[,all.vars(x$formula)[2],drop=T]
        if('f_summary' %in% names(x)){
            #to generate label - latex2exp::TeX('$\\textit{b+\\beta(h)}$','character')
            y_lab <- "paste('','',italic(paste('b+',beta,,,,'(','h',')','')),'')"
            plot_dat <- x$f_summary[x$f_summary$h>=min(h_in_data) & x$f_summary$h<=max(h_in_data),]
        }else{
            #to generate label - latex2exp::TeX('$\\textit{b}$','character')
            y_lab <- "paste('','',italic(paste('b')),'')"
            plot_dat <- data.frame(h=x$data[,all.vars(x$formula)[2],drop=T],
                                   lower=x$param_summary['b','lower'],
                                   median=x$param_summary['b','median'],
                                   upper=x$param_summary['b','upper'])
        }
        p <- ggplot(data=plot_dat) +
            geom_path(aes(.data$h,.data$median)) +
            geom_path(aes(.data$h,.data$lower),linetype='dashed') +
            geom_path(aes(.data$h,.data$upper),linetype='dashed') +
            xlab(parse(text=x_lab)) +
            ylab(parse(text=y_lab)) +
            theme_bdrc()
    }else if(type=='residuals'){
        resid_dat <- merge(x$rating_curve[,c('h','median')],x$data,by.x='h',by.y=all.vars(x$formula)[2],)
        if('sigma_eps_summary' %in% names(x)){
            resid_dat <- merge(resid_dat,x$sigma_eps_summary[,c('h','median')],by = 'h')
            names(resid_dat) <- c('h','median','Q','sigma_eps')
        }else{
            resid_dat$sigma_eps <- x$param_summary['sigma_eps','median']
        }
        c_hat <- if(is.null(x$run_info$c_param)) median(x$c_posterior) else x$run_info$c_param
        resid_dat[,'log(h-c_hat)'] <- log(resid_dat$h-c_hat)
        resid_dat$r_median <- log(resid_dat$Q)-log(resid_dat$median)
        resid_dat$r_lower <- -1.96*resid_dat$sigma_eps
        resid_dat$r_upper <- 1.96*resid_dat$sigma_eps
        #to generate label - latex2exp::TeX("$log(\\textit{Q})-log(\\textit{\\hat{Q}})$",'character')
        y_lab <- "paste('','log','(','',italic(paste('Q')),')','','-log','(','',italic(paste('',hat(paste('Q')))),')','','')"
        #to generate label - latex2exp::TeX("$log(\\textit{h - \\hat{c}})$",'character')
        x_lab <- "paste('','log','(','',italic(paste('h',phantom() - phantom(),'',hat(paste('c')))),')','','')"
        p <- ggplot(data=resid_dat) +
            geom_point(aes(.data$`log(h-c_hat)`,.data$r_median),size=2) +
            geom_path(aes(.data$`log(h-c_hat)`,.data$r_lower),linetype='dashed') +
            geom_path(aes(.data$`log(h-c_hat)`,.data$r_upper),linetype='dashed') +
            geom_abline(intercept=0,slope=0,size=1.1) +
            xlab(parse(text=x_lab)) +
            ylab(parse(text=y_lab)) +
            theme_bdrc()
    }
    if(!is.null(title)){
        p <- p + ggtitle(title)
    }
    return(p)
}


#' @importFrom gridExtra arrangeGrob
plot_collage <- function(x,transformed=F){
    types <- c('rating_curve','residuals','f','sigma_eps')
    plot_list <- lapply(types,function(ty){
        plot_fun(x,type=ty,transformed=transformed)
    })
    p <- do.call(arrangeGrob,c(plot_list,ncol=round(sqrt(length(types)))))
    return(p)
}

predict_fun <- function(object,newdata=NULL){
    if(is.null(newdata)){
        merged_data <- merge(object$rating_curve,object$data,by.x='h',by.y=all.vars(object$formula)[2])
        pred_dat <- merged_data[,c('h','lower','median','upper')]
    }else{
        if(class(newdata) !='numeric'){
            stop('newdata must be a vector of type "numeric" or NULL')
        }
        if(any(is.na(newdata))){
            stop('newdata must include NA')
        }
        if(any(newdata<min(object$rating_curve$h) | newdata>max(object$rating_curve$h))){
            stop('newdata must contain values within the range of stage values used to fit the rating curve. See "h_max" option to extrapolate the rating curve to higher stages')
        }
        lower_pred <- stats::approx(object$rating_curve$h,object$rating_curve$lower,xout=newdata)$y
        median_pred <- stats::approx(object$rating_curve$h,object$rating_curve$median,xout=newdata)$y
        upper_pred <- stats::approx(object$rating_curve$h,object$rating_curve$upper,xout=newdata)$y
        pred_dat <- data.frame(h=newdata,lower=lower_pred,median=median_pred,upper=upper_pred)
    }
    return(pred_dat)
}

#' Print plm0 object
#'
#' Print the results of a plm0 object
#' @param x an object of class "plm0"
#' @param ... not used in this function
#' @seealso \code{\link{plm0}} for fitting the plm0 model, \code{\link{summary.plm0}} for summaries, \code{\link{predict.plm0}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{plot.plm0}} to help visualize the full posterior distributions.
#' @examples
#' \dontrun{
#' data(V316_river)
#' f <- Q~W
#' plm0.fit <- plm0(f,V316_river)
#' print(plm0.fit)
#' }
#' @export
print.plm0 <- function(x,...){
    print_fun(x)
}

#' Summarizing plm0 fit
#'
#' Summarize the results of a plm0 object
#' @param object an object of class "plm0"
#' @param ... Not used for this function
#' @seealso \code{\link{plm0}} for fitting the plm0 model, \code{\link{predict.plm0}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{plot.plm0}} to help visualize the full posterior distributions.
#' @examples
#' \dontrun{
#' data(V316_river)
#' f <- Q~W
#' plm0.fit <- plm0(f,V316_river)
#' summary(plm0.fit)
#' }
#' @export
summary.plm0 <- function(object,...){
    summary_fun(object)
}

#' Autoplot plm0 fit
#'
#' Uses ggplot2 to plot plm0 object
#' @param x an object of class "plm0".
#' @param type a character denoting what type of plot should be drawn. Defaults to "rating_curve". Possible types are
#'                    \itemize{
#'                       \item{"rating_curve"}{ to plot the rating curve.}
#'                       \item{"rating_curve_mean"}{ to plot the posterior mean of the rating curve.}
#'                       \item{"f"}{ to plot the power-law exponent.}
#'                       \item{"beta"}{ to plot the random effect in the power-law exponent.}
#'                       \item{"sigma_eps"}{ to plot the standard deviation on the data level.}
#'                       \item{"residuals"}{ to plot the log residuals.}
#'                       \item{"trace"}{ to plot trace plots of parameters given in param.}
#'                       \item{"histogram"}{ to plot histograms of parameters given in param.}
#'                    }
#' @param param a character vector with the parameters to plot. Defaults to NULL and is only used if type is "trace" or "histogram". Allowed values are the parameters given in the model summary of x as well as "hyperparameters" or "latent_parameters" for specific groups of parameters.
#' @param transformed a logical value indicating whether the quantity should be plotted on a transformed scale used during the Bayesian inference. Defaults to FALSE.
#' @param title a character denoting the title of the plot. Defaults to NULL, i.e. no title.
#' @param ... further arguments passed to other methods (currently unused).
#' @return returns an object of class ggplot2.
#' @seealso \code{\link{plm0}} for fitting the plm0 model,\code{\link{summary.plm0}} for summaries of model parameters, \code{\link{predict.plm0}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{plot.plm0}} to help visualize the full posterior distributions.
#' @examples
#' \dontrun{
#' data(V316_river)
#' f <- Q~W
#' plm0.fit <- plm0(f,V316_river)
#' autoplot(plm0.fit)
#' }
#' @export
autoplot.plm0 <- function(x,type='rating_curve',param=NULL,transformed=F,title=NULL,...){
    plot_fun(x,type=type,param=param,transformed=transformed,title=title)
}

#' Plot plm0 fit
#'
#' Print the results of a  object
#' @param x an object of class "plm0".
#' @param type a character denoting what type of plot should be drawn. Defaults to "rating_curve". Possible types are
#'                    \itemize{
#'                       \item{"rating_curve"}{ to plot the rating curve.}
#'                       \item{"rating_curve_mean"}{ to plot the posterior mean of the rating curve.}
#'                       \item{"f"}{ to plot the power-law exponent.}
#'                       \item{"beta"}{ to plot the random effect in the power-law exponent.}
#'                       \item{"sigma_eps"}{ to plot the standard deviation on the data level.}
#'                       \item{"residuals"}{ to plot the log residuals.}
#'                       \item{"trace"}{ to plot trace plots of parameters given in param.}
#'                       \item{"histogram"}{ to plot histograms of parameters given in param.}
#'                    }
#' @param param a character vector with the parameters to plot. Defaults to NULL and is only used if type is "trace" or "histogram". Allowed values are the parameters given in the model summary of x as well as "hyperparameters" or "latent_parameters" for specific groups of parameters.
#' @param transformed a logical value indicating whether the quantity should be plotted on a transformed scale used during the Bayesian inference. Defaults to FALSE.
#' @param title a character denoting the title of the plot. Defaults to NULL, i.e. no title.
#' @param ... further arguments passed to other methods (currently unused).
#' @seealso \code{\link{plm0}} for fitting the plm0 model,\code{\link{summary.plm0}} for summaries, \code{\link{predict.plm0}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{plot.plm0}} to help visualize the full posterior distributions.
#' @examples
#' \dontrun{
#' data(V316_river)
#' f <- Q~W
#' plm0.fit <- plm0(f,V316_river)
#' plot(plm0.fit)
#' }
#' @export
#' @importFrom grid grid.draw
#' @importFrom ggplot2 autoplot
plot.plm0 <- function(x,type='rating_curve',param=NULL,transformed=F,title=NULL,...){
    if(is.null(type) || type!='collage'){
        p <- autoplot(x,type=type,param=param,transformed=transformed,title=title,...)
        print(p)
    }else{
        p <- plot_collage(x,transformed=transformed)
        grid.draw(p)
    }
}

#' Predict method for plm0 fit
#'
#' Print the results of a  object
#' @param object an object of class "plm0"
#' @param newdata a numeric vector of stage values for which to predict. If omitted, the stage values in the data are used.
#' @param ... not used in this function
#' @return numeric vector of discharge values for the stage values given in newdata
#' @seealso \code{\link{plm0}} for fitting the plm0 model,\code{\link{summary.plm0}} for summaries, \code{\link{predict.plm0}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{plot.plm0}} to help visualize the full posterior distributions.
#' @examples
#' \dontrun{
#' data(V316_river)
#' f <- Q~W
#' plm0.fit <- plm0(f,V316_river,h_max=2)
#' #predict rating curve on a equally 1 cm spaced grid from 1 to 2 meters
#' predict(plm0.fit,newdata=seq(1,2,by=0.01))
#' }
#' @export
predict.plm0 <- function(object,newdata=NULL,...){
    predict_fun(object,newdata)
}

#' Print plm object
#'
#' Print the results of a plm object
#' @param x an object of class "plm"
#' @param ... not used in this function
#' @return gplm0 returns an object of class "plm"\cr\cr
#' @seealso \code{\link{summary.plm}} for summaries, \code{\link{predict.plm}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{plot.plm}} to help visualize the full posterior distributions.
#' @examples
#' \dontrun{
#' data(V316_river)
#' f <- Q~W
#' plm.fit <- plm(f,V316_river)
#' print(plm.fit)
#' }
#' @export
#'
print.plm <- function(x,...){
    print_fun(x)
}

#' Summarizing plm fit
#'
#' Summarize the results of a plm object
#' @param object an object of class "plm"
#' @param ... not used in this function
#' @return gplm0 returns an object of class "plm"\cr\cr
#' @seealso \code{\link{plm}} for fitting the plm model,\code{\link{summary.plm}} for summaries, \code{\link{predict.plm}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{plot.plm}} to help visualize the full posterior distributions.
#' @examples
#' \dontrun{
#' data(V316_river)
#' f <- Q~W
#' plm.fit <- plm(f,V316_river)
#' summary(plm.fit)
#' }
#' @export
summary.plm <- function(object,...){
    summary_fun(object)
}

#' Autoplot plm fit
#'
#' Uses ggplot2 to plot plm object
#' @param x an object of class "plm"
#' @param type a character denoting what type of plot should be drawn. Defaults to "rating_curve". Possible types are
#'                    \itemize{
#'                       \item{"rating_curve"}{ to plot the rating curve.}
#'                       \item{"rating_curve_mean"}{ to plot the posterior mean of the rating curve.}
#'                       \item{"f"}{ to plot the power-law exponent.}
#'                       \item{"beta"}{ to plot the random effect in the power-law exponent.}
#'                       \item{"sigma_eps"}{ to plot the standard deviation on the data level.}
#'                       \item{"residuals"}{ to plot the log residuals.}
#'                       \item{"trace"}{ to plot trace plots of parameters given in param.}
#'                       \item{"histogram"}{ to plot histograms of parameters given in param.}
#'                    }
#' @param param a character vector with the parameters to plot. Defaults to NULL and is only used if type is "trace" or "histogram". Allowed values are the parameters given in the model summary of x as well as "hyperparameters" or "latent_parameters" for specific groups of parameters.
#' @param transformed a logical value indicating whether the quantity should be plotted on a transformed scale used during the Bayesian inference. Defaults to FALSE.
#' @param title a character denoting the title of the plot. Defaults to NULL, i.e. no title.
#' @param ... further arguments passed to other methods (currently unused).
#' @return returns an object of class ggplot2
#' @seealso \code{\link{plm}} for fitting the plm model,\code{\link{summary.plm}} for summaries of model parameters, \code{\link{predict.plm}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{plot.plm}} to help visualize the full posterior distributions.
#' @examples
#' \dontrun{
#' data(V316_river)
#' f <- Q~W
#' plm.fit <- plm(f,V316_river)
#' autoplot(plm.fit)
#' }
#' @export
autoplot.plm <- function(x,type='rating_curve',param=NULL,transformed=F,title=NULL,...){
    plot_fun(x,type=type,param=param,transformed=transformed,title=title)
}

#' Plot plm fit
#'
#' Print the results of a  object
#' @param x an object of class "plm"
#' @param type a character denoting what type of plot should be drawn. Defaults to "rating_curve". Possible types are
#'                    \itemize{
#'                       \item{"rating_curve"}{ to plot the rating curve.}
#'                       \item{"rating_curve_mean"}{ to plot the posterior mean of the rating curve.}
#'                       \item{"f"}{ to plot the power-law exponent.}
#'                       \item{"beta"}{ to plot the random effect in the power-law exponent.}
#'                       \item{"sigma_eps"}{ to plot the standard deviation on the data level.}
#'                       \item{"residuals"}{ to plot the log residuals.}
#'                       \item{"trace"}{ to plot trace plots of parameters given in param.}
#'                       \item{"histogram"}{ to plot histograms of parameters given in param.}
#'                    }
#' @param param a character vector with the parameters to plot. Defaults to NULL and is only used if type is "trace" or "histogram". Allowed values are the parameters given in the model summary of x as well as "hyperparameters" or "latent_parameters" for specific groups of parameters.
#' @param transformed a logical value indicating whether the quantity should be plotted on a transformed scale used during the Bayesian inference. Defaults to FALSE.
#' @param title a character denoting the title of the plot. Defaults to NULL, i.e. no title.
#' @param ... further arguments passed to other methods (currently unused).
#' @seealso \code{\link{plm}} for fitting the plm model,\code{\link{summary.plm}} for summaries, \code{\link{predict.plm}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{plot.plm}} to help visualize the full posterior distributions.
#' @examples
#' \dontrun{
#' data(V316_river)
#' f <- Q~W
#' plm.fit <- plm(f,V316_river)
#' plot(plm.fit)
#' }
#' @export
#' @importFrom grid grid.draw
#' @importFrom ggplot2 autoplot
plot.plm <- function(x,type='rating_curve',param=NULL,transformed=F,title=NULL,...){
    if(is.null(type) || type!='collage'){
        p <- autoplot(x,type=type,param=param,transformed=transformed,title=title,...)
        print(p)
    }else{
        p <- plot_collage(x,transformed=transformed)
        grid.draw(p)
    }
}

#' Predict method for plm fit
#'
#' Print the results of a  object
#' @param object an object of class "plm"
#' @param newdata a numeric vector of stage values for which to predict. If omitted, the stage values in the data are used.
#' @param ... not used in this function
#' @return numeric vector of discharge values for the stage values given in newdata
#' @seealso \code{\link{plm}} for fitting the plm model,\code{\link{summary.plm}} for summaries, \code{\link{predict.plm}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{plot.plm}} to help visualize the full posterior distributions.
#' @examples
#' \dontrun{
#' data(V316_river)
#' f <- Q~W
#' plm.fit <- plm(f,V316_river,h_max=2)
#' #predict rating curve on a equally 1 cm spaced grid from 1 to 2 meters
#' predict(plm.fit,newdata=seq(1,2,by=0.01))
#' }
#' @export
predict.plm <- function(object,newdata=NULL,...){
    predict_fun(object,newdata)
}

#' Print gplm0 object
#'
#' Print the results of a gplm0 object
#' @param x an object of class "gplm0"
#' @param ... not used in this function
#' @seealso \code{\link{summary.gplm0}} for summaries, \code{\link{predict.gplm0}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{plot.gplm0}} to help visualize the full posterior distributions.
#' @examples
#' \dontrun{
#' data(V316_river)
#' f <- Q~W
#' gplm0.fit <- gplm0(f,V316_river)
#' print(gplm0.fit)
#' }
#' @export
print.gplm0 <- function(x,...){
    print_fun(x)
}

#' Summarizing gplm0 fit
#'
#' Summarize the results of a gplm0 object
#' @param object an object of class "gplm0"
#' @param ... not used in this function
#' @return gplm0 returns an object of class "plm"\cr\cr
#' @seealso \code{\link{gplm0}} for fitting the gplm0 model,\code{\link{summary.gplm0}} for summaries, \code{\link{predict.gplm0}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{plot.gplm0}} to help visualize the full posterior distributions.
#' @examples
#' \dontrun{
#' data(V316_river)
#' f <- Q~W
#' gplm0.fit <- gplm0(f,V316_river)
#' summary(gplm0.fit)
#' }
#' @export
summary.gplm0 <- function(object,...){
    summary_fun(object)
}

#' Autoplot gplm0 fit
#'
#' Uses ggplot2 to plot gplm0 object
#' @param x an object of class "gplm0"
#' @param type a character denoting what type of plot should be drawn. Defaults to "rating_curve". Possible types are
#'                    \itemize{
#'                       \item{"rating_curve"}{ to plot the rating curve.}
#'                       \item{"rating_curve_mean"}{ to plot the posterior mean of the rating curve.}
#'                       \item{"f"}{ to plot the power-law exponent.}
#'                       \item{"beta"}{ to plot the random effect in the power-law exponent.}
#'                       \item{"sigma_eps"}{ to plot the standard deviation on the data level.}
#'                       \item{"residuals"}{ to plot the log residuals.}
#'                       \item{"trace"}{ to plot trace plots of parameters given in param.}
#'                       \item{"histogram"}{ to plot histograms of parameters given in param.}
#'                    }
#' @param param a character vector with the parameters to plot. Defaults to NULL and is only used if type is "trace" or "histogram". Allowed values are the parameters given in the model summary of x as well as "hyperparameters" or "latent_parameters" for specific groups of parameters.
#' @param transformed a logical value indicating whether the quantity should be plotted on a transformed scale used during the Bayesian inference. Defaults to FALSE.
#' @param title a character denoting the title of the plot. Defaults to NULL, i.e. no title.
#' @param ... further arguments passed to other methods (currently unused).
#' @return returns an object of class ggplot2
#' @seealso \code{\link{gplm0}} for fitting the gplm0 model,\code{\link{summary.gplm0}} for summaries of model parameters, \code{\link{predict.gplm0}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{plot.gplm0}} to help visualize the full posterior distributions.
#' @examples
#' \dontrun{
#' data(V316_river)
#' f <- Q~W
#' gplm0.fit <- gplm0(f,V316_river)
#' autoplot(gplm0.fit)
#' }
#' @export
autoplot.gplm0 <- function(x,type='rating_curve',param=NULL,transformed=F,title=NULL,...){
    plot_fun(x,type=type,param=param,transformed=transformed,title=title)
}

#' Plot gplm0 fit
#'
#' Print the results of a  object
#' @param x an object of class "gplm0"
#' @param type a character denoting what type of plot should be drawn. Defaults to "rating_curve". Possible types are
#'                    \itemize{
#'                       \item{"rating_curve"}{ to plot the rating curve.}
#'                       \item{"rating_curve_mean"}{ to plot the posterior mean of the rating curve.}
#'                       \item{"f"}{ to plot the power-law exponent.}
#'                       \item{"beta"}{ to plot the random effect in the power-law exponent.}
#'                       \item{"sigma_eps"}{ to plot the standard deviation on the data level.}
#'                       \item{"residuals"}{ to plot the log residuals.}
#'                       \item{"trace"}{ to plot trace plots of parameters given in param.}
#'                       \item{"histogram"}{ to plot histograms of parameters given in param.}
#'                    }
#' @param param a character vector with the parameters to plot. Defaults to NULL and is only used if type is "trace" or "histogram". Allowed values are the parameters given in the model summary of x as well as "hyperparameters" or "latent_parameters" for specific groups of parameters.
#' @param transformed a logical value indicating whether the quantity should be plotted on a transformed scale used during the Bayesian inference. Defaults to FALSE.
#' @param title a character denoting the title of the plot. Defaults to NULL, i.e. no title.
#' @param ... further arguments passed to other methods (currently unused).
#' @seealso \code{\link{gplm0}} for fitting the gplm0 model,\code{\link{summary.gplm0}} for summaries, \code{\link{predict.gplm0}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{plot.gplm0}} to help visualize the full posterior distributions.
#' @examples
#' \dontrun{
#' data(V316_river)
#' f <- Q~W
#' gplm0.fit <- gplm0(f,V316_river)
#' plot(gplm0.fit)
#' }
#' @export
#' @importFrom grid grid.draw
#' @importFrom ggplot2 autoplot
plot.gplm0 <- function(x,type='rating_curve',param=NULL,transformed=F,title=NULL,...){
    if(is.null(type) || type!='collage'){
        p <- autoplot(x,type=type,param=param,transformed=transformed,title=title,...)
        print(p)
    }else{
        p <- plot_collage(x,transformed=transformed)
        grid.draw(p)
    }
}

#' Predict method for gplm0 fit
#'
#' Print the results of a  object
#' @param object an object of class "gplm0"
#' @param newdata a numeric vector of stage values for which to predict. If omitted, the stage values in the data are used.
#' @param ... not used in this function
#' @return numeric vector of discharge values for the stage values given in newdata
#' @seealso \code{\link{gplm0}} for fitting the gplm0 model,\code{\link{summary.gplm0}} for summaries, \code{\link{predict.gplm0}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{plot.gplm0}} to help visualize the full posterior distributions.
#' @examples
#' \dontrun{
#' data(V316_river)
#' f <- Q~W
#' gplm0.fit <- gplm0(f,V316_river,h_max=2)
#' #predict rating curve on a equally 1 cm spaced grid from 1 to 2 meters
#' predict(gplm0.fit,newdata=seq(1,2,by=0.01))
#' }
#' @export
predict.gplm0 <- function(object,newdata=NULL,...){
    predict_fun(object,newdata)
}

#' Print gplm object
#'
#' Print the results of a gplm object
#' @param x an object of class "gplm"
#' @param ... not used in this function
#' @seealso \code{\link{summary.gplm}} for summaries, \code{\link{predict.gplm}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{plot.gplm}} to help visualize the full posterior distributions.
#' @examples
#' \dontrun{
#' data(V316_river)
#' f <- Q~W
#' gplm.fit <- gplm(f,V316_river)
#' print(gplm.fit)
#' }
#' @export
print.gplm <- function(x,...){
    print_fun(x)
}

#' Summarizing plm fit
#'
#' Summarize the results of a gplm object
#' @param object an object of class "gplm"
#' @param ... not used in this function
#' @return gplm0 returns an object of class "plm"\cr\cr
#' @seealso \code{\link{gplm}} for fitting the gplm model,\code{\link{summary.gplm}} for summaries, \code{\link{predict.gplm}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{plot.gplm}} to help visualize the full posterior distributions.
#' @examples
#' \dontrun{
#' data(V316_river)
#' f <- Q~W
#' gplm.fit <- gplm(f,V316_river)
#' summary(gplm.fit)
#' }
#' @export
summary.gplm <- function(object,...){
    summary_fun(object)
}

#' Autoplot gplm fit
#'
#' Uses ggplot2 to plot gplm object
#'
#' @param x an object of class "gplm"
#' @param type a character denoting what type of plot should be drawn. Defaults to "rating_curve". Possible types are
#'                    \itemize{
#'                       \item{"rating_curve"}{ to plot the rating curve.}
#'                       \item{"rating_curve_mean"}{ to plot the posterior mean of the rating curve.}
#'                       \item{"f"}{ to plot the power-law exponent.}
#'                       \item{"beta"}{ to plot the random effect in the power-law exponent.}
#'                       \item{"sigma_eps"}{ to plot the standard deviation on the data level.}
#'                       \item{"residuals"}{ to plot the log residuals.}
#'                       \item{"trace"}{ to plot trace plots of parameters given in param.}
#'                       \item{"histogram"}{ to plot histograms of parameters given in param.}
#'                    }
#' @param param a character vector with the parameters to plot. Defaults to NULL and is only used if type is "trace" or "histogram". Allowed values are the parameters given in the model summary of x as well as "hyperparameters" or "latent_parameters" for specific groups of parameters.
#' @param transformed a logical value indicating whether the quantity should be plotted on a transformed scale used during the Bayesian inference. Defaults to FALSE.
#' @param title a character denoting the title of the plot. Defaults to NULL, i.e. no title.
#' @param ... further arguments passed to other methods (currently unused).
#' @return returns an object of class ggplot2
#' @seealso \code{\link{gplm}} for fitting the gplm model,\code{\link{summary.gplm}} for summaries of model parameters, \code{\link{predict.gplm}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{plot.gplm}} to help visualize the full posterior distributions.
#' @examples
#' \dontrun{
#' data(V316_river)
#' f <- Q~W
#' gplm.fit <- gplm(f,V316_river)
#' autoplot(gplm.fit)
#' }
#' @export
autoplot.gplm <- function(x,type='rating_curve',param=NULL,transformed=F,title=NULL,...){
    plot_fun(x,type=type,param=param,transformed=transformed,title=title)
}

#' Plot gplm fit
#'
#' Print the results of a  object
#'
#' @param x an object of class "gplm"
#' @param type a character denoting what type of plot should be drawn. Defaults to "rating_curve". Possible types are
#'                    \itemize{
#'                       \item{"rating_curve"}{ to plot the rating curve.}
#'                       \item{"rating_curve_mean"}{ to plot the posterior mean of the rating curve.}
#'                       \item{"f"}{ to plot the power-law exponent.}
#'                       \item{"beta"}{ to plot the random effect in the power-law exponent.}
#'                       \item{"sigma_eps"}{ to plot the standard deviation on the data level.}
#'                       \item{"residuals"}{ to plot the log residuals.}
#'                       \item{"trace"}{ to plot trace plots of parameters given in param.}
#'                       \item{"histogram"}{ to plot histograms of parameters given in param.}
#'                    }
#' @param param a character vector with the parameters to plot. Defaults to NULL and is only used if type is "trace" or "histogram". Allowed values are the parameters given in the model summary of x as well as "hyperparameters" or "latent_parameters" for specific groups of parameters.
#' @param transformed a logical value indicating whether the quantity should be plotted on a transformed scale used during the Bayesian inference. Defaults to FALSE.
#' @param title a character denoting the title of the plot. Defaults to NULL, i.e. no title.
#' @param ... further arguments passed to other methods (currently unused).
#' @seealso \code{\link{gplm}} for fitting the gplm model,\code{\link{summary.gplm}} for summaries, \code{\link{predict.gplm}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{plot.gplm}} to help visualize the full posterior distributions.
#' @examples
#' \dontrun{
#' data(V316_river)
#' f <- Q~W
#' gplm.fit <- gplm(f,V316_river)
#' plot(gplm.fit)
#' }
#' @export
#' @importFrom grid grid.draw
#' @importFrom ggplot2 autoplot
plot.gplm <- function(x,type='rating_curve',param=NULL,transformed=F,title=NULL,...){
    if(is.null(type) || type!='collage'){
        p <- autoplot(x,type=type,param=param,transformed=transformed,title=title,...)
        print(p)
    }else{
        p <- plot_collage(x,transformed=transformed)
        grid.draw(p)
    }
}

#' Predict method for gplm fit
#'
#' Print the results of a  object
#' @param object an object of class "gplm"
#' @param newdata a numeric vector of stage values for which to predict. If omitted, the stage values in the data are used.
#' @param ... not used in this function
#' @return numeric vector of discharge values for the stage values given in newdata
#' @seealso \code{\link{gplm}} for fitting the gplm model,\code{\link{summary.gplm}} for summaries, \code{\link{predict.gplm}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{plot.gplm}} to help visualize the full posterior distributions.
#' @examples
#' \dontrun{
#' data(V316_river)
#' f <- Q~W
#' gplm0.fit <- gplm0(f,V316_river,h_max=2)
#' #predict rating curve on a equally 1 cm spaced grid from 1 to 2 meters
#' predict(gplm0.fit,newdata=seq(1,2,by=0.01))
#' }
#' @export
predict.gplm <- function(object,newdata=NULL,...){
    predict_fun(object,newdata)
}
