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


plot_fun <- function(x,type=NULL,...,transformed=F,title=NULL){
    args <- c(...)
    if(type=='trace'){
        plot_dat <- gather_draws(x,args,transformed=transformed)
        if('h' %in% names(plot_dat)){
            stop('Plots of type "trace" can only be of stage-independent parameters')
        }
        params <- unique(plot_dat$name)
        if(length(params)>1){
            param_levels <- get_parameter_levels(params)
            plot_dat$name_expr <- factor(plot_dat$name,levels=param_levels,labels=sapply(param_levels,get_param_expression))
            plot_dat$chain <- factor(as.character(plot_dat$chain),levels=1:max(plot_dat$chain))
            p <- ggplot(plot_dat,aes(iter,value,col=chain)) +
                geom_line() +
                facet_wrap(~name_expr,scales='free',labeller = label_parsed) +
                scale_color_manual(values=c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF"),
                                   name='Chain number') +
                xlab('Iteration') +
                ylab('') +
                theme_classic() +
                theme(strip.background = element_blank(),
                      strip.text.x = element_text(size = 16),
                      axis.title.x = element_text(size=16),
                      axis.text.x = element_text(size=12),
                      axis.text.y = element_text(size=12),
                      legend.text=element_text(size=12),
                      legend.title=element_text(size=16),
                )
        }else{
            param_expr <- get_param_expression(params)
            plot_dat$chain_name <- paste0('Chain nr ',plot_dat$chain)
            p <- ggplot(plot_dat,aes(iter,value)) +
                geom_line(col="#0072B5FF") +
                facet_wrap(~chain_name,scales='free') +
                xlab('Iteration') +
                ylab(parse(text=param_expr)) +
                theme_classic() +
                theme(strip.background = element_blank(),
                      strip.text.x = element_text(size = 16),
                      axis.title.x = element_text(size=16),
                      axis.title.y = element_text(size=16),
                      axis.text.x = element_text(size=12),
                      axis.text.y = element_text(size=12),
                      legend.text=element_text(size=16),
                      legend.title=element_text(size=16),
                )
        }
    }else if(type=='histogram'){
        plot_dat <- gather_draws(x,args,transformed=transformed)
        if('h' %in% names(plot_dat)){
            stop('Plots of type "histogram" can only be of stage-independent parameters')
        }
        params <- unique(plot_dat$name)
        param_levels <- get_parameter_levels(params)
        plot_dat$name_expr <- factor(plot_dat$name,levels=param_levels,labels=sapply(param_levels,get_param_expression))
        plot_dat$chain <- factor(as.character(plot_dat$chain),levels=1:max(plot_dat$chain))
        p <- ggplot(plot_dat,aes(value)) +
            geom_histogram(bins=100,fill="#0072B5FF") +
            facet_wrap(~name_expr,scales='free',labeller = label_parsed) +
            xlab('') +
            ylab('') +
            theme_classic() +
            theme(strip.background = element_blank(),
                  strip.text.x = element_text(size = 16),
                  axis.text.x = element_text(size=12),
                  axis.text.y = element_text(size=12))
    }else if(type=='rating_curve' | type=='rating_curve_mean'){
        if(transformed){
            x_lab <- latex2exp::TeX('$\\log(\\textit{h-\\hat{c}})$','character')
            y_lab <- latex2exp::TeX('$\\log(\\textit{Q})$','character')
            c_hat <- if(is.null(x$run_info_c_param)) median(x$c_posterior) else x$run_info$c_param
            plot_dat <- merge(x[[type]],x$data,by.x='h',by.y=all.vars(x$formula)[2])
            plot_dat[,'log(h-c_hat)'] <- log(plot_dat$h-c_hat)
            p <- ggplot(data=plot_dat) +
                geom_point(aes(`log(h-c_hat)`,log(Q))) +
                geom_line(aes(`log(h-c_hat)`,log(median))) +
                geom_line(aes(`log(h-c_hat)`,log(lower)),linetype='dashed') +
                geom_line(aes(`log(h-c_hat)`,log(upper)),linetype='dashed') +
                xlab(parse(text=x_lab)) +
                ylab(parse(text=y_lab)) +
                theme_classic() +
                theme(axis.text.x = element_text(size = 12),
                      axis.text.y = element_text(size = 12),
                      axis.title.x = element_text(size = 16),
                      axis.title.y = element_text(size = 16))
        }else{
            x_lab <- latex2exp::TeX('$\\textit{Q}\\lbrack\\textit{m^3/s}\\rbrack$','character')
            y_lab <- latex2exp::TeX('$\\textit{h}\\lbrack\\textit{m}\\rbrack$','character')
            p <- ggplot(data=x[[type]]) +
                geom_point(data=x$data,aes_string(all.vars(x$formula)[1],all.vars(x$formula)[2])) +
                geom_line(aes(median,h)) +
                geom_line(aes(lower,h),linetype='dashed') +
                geom_line(aes(upper,h),linetype='dashed') +
                scale_x_continuous(expand=c(0,0)) +
                scale_y_continuous(expand=c(0,0)) +
                xlab(parse(text=x_lab)) +
                ylab(parse(text=y_lab)) +
                theme_classic() +
                theme(axis.text.x = element_text(size = 12),
                      axis.text.y = element_text(size = 12),
                      axis.title.x = element_text(size = 16),
                      axis.title.y = element_text(size = 16))
        }
    }else if(type=='sigma_eps'){
        x_lab <- latex2exp::TeX('$\\textit{h}\\lbrack\\textit{m}\\rbrack$','character')
        if('sigma_eps_summary' %in% names(x)){
            y_lab <- latex2exp::TeX('$\\sigma_{\\epsilon}(\\textit{h})$','character')
            plot_dat <- x$sigma_eps_summary
        }else{
            y_lab <- latex2exp::TeX('$\\sigma_{\\epsilon}$','character')
            plot_dat <- data.frame(h=x$rating_curve$h,
                                   lower=x$param_summary['sigma_eps','lower'],
                                   median=x$param_summary['sigma_eps','median'],
                                   upper=x$param_summary['sigma_eps','upper'])
        }
        #y_lab <- 'sigma[epsilon](h)'
        p <- ggplot(data=plot_dat) +
            geom_line(aes(h,median)) +
            geom_line(aes(h,lower),linetype='dashed') +
            geom_line(aes(h,upper),linetype='dashed') +
            xlab(parse(text=x_lab)) +
            ylab(parse(text=y_lab)) +
            theme_classic() +
            theme(axis.text.x = element_text(size = 12),
                  axis.text.y = element_text(size = 12),
                  axis.title.x = element_text(size = 16),
                  axis.title.y = element_text(size = 16))
    }else if(type=='beta'){
        if(!('beta_summary' %in% names(x))){
            stop('Plots of type "beta" are only for models with stage dependent power law exponent, s.a. "bgplm0" and "bgplm"')
        }
        x_lab <- latex2exp::TeX('$\\textit{h}\\lbrack\\textit{m}\\rbrack$','character')
        y_lab <- latex2exp::TeX('$\\beta(\\textit{h})$')
        p <- ggplot(data=x$beta_summary) +
            geom_line(aes(h,median)) +
            geom_line(aes(h,lower),linetype='dashed') +
            geom_line(aes(h,upper),linetype='dashed') +
            xlab(parse(text=x_lab)) +
            ylab(parse(text=y_lab)) +
            theme_classic() +
            theme(axis.text.x = element_text(size = 12),
                  axis.text.y = element_text(size = 12),
                  axis.title.x = element_text(size = 16),
                  axis.title.y = element_text(size = 16))
    }else if(type=='f'){
        x_lab <- latex2exp::TeX('$\\textit{h}\\lbrack\\textit{m}\\rbrack$','character')
        if('f_summary' %in% names(x)){
            y_lab <- latex2exp::TeX('$\\textit{b+\\beta(h)}$')
            plot_dat <- x$f_summary
        }else{
            y_lab <- latex2exp::TeX('$\\textit{b}$')
            plot_dat <- data.frame(h=x$rating_curve$h,
                                   lower=x$param_summary['b','lower'],
                                   median=x$param_summary['b','median'],
                                   upper=x$param_summary['b','upper'])
        }
        p <- ggplot(data=plot_dat) +
            geom_line(aes(h,median)) +
            geom_line(aes(h,lower),linetype='dashed') +
            geom_line(aes(h,upper),linetype='dashed') +
            xlab(parse(text=x_lab)) +
            ylab(parse(text=y_lab)) +
            theme_classic() +
            theme(axis.text.x = element_text(size = 12),
                  axis.text.y = element_text(size = 12),
                  axis.title.x = element_text(size = 16),
                  axis.title.y = element_text(size = 16))
    }else if(type=='residuals'){
        resid_dat <- merge(x$rating_curve[,c('h','median')],x$data,by.x='h',by.y=all.vars(x$formula)[2],)
        if('sigma_eps_summary' %in% names(x)){
            resid_dat <- merge(resid_dat,x$sigma_eps_summary[,c('h','median')],by = 'h')
            names(resid_dat) <- c('h','median','Q','sigma_eps')
        }else{
            resid_dat$sigma_eps <- x$param_summary['sigma_eps','median']
        }
        c_hat <- if(is.null(x$run_info_c_param)) median(x$c_posterior) else x$run_info$c_param
        resid_dat[,'log(h-c_hat)'] <- log(resid_dat$h-c_hat)
        resid_dat$r_median <- log(resid_dat$Q)-log(resid_dat$median)
        resid_dat$r_lower <- -1.96*resid_dat$sigma_eps
        resid_dat$r_upper <- 1.96*resid_dat$sigma_eps
        y_lab <- TeX("$log(\\textit{Q})-log(\\textit{\\hat{Q}})$")
        x_lab <- TeX("$log(\\textit{h - \\hat{c}})$")
        p <- ggplot(data=resid_dat) +
            geom_point(aes(`log(h-c_hat)`,r_median),size=2) +
            geom_line(aes(`log(h-c_hat)`,r_lower),linetype='dashed') +
            geom_line(aes(`log(h-c_hat)`,r_upper),linetype='dashed') +
            geom_abline(intercept=0,slope=0,size=1.1) +
            xlab(parse(text=x_lab)) +
            ylab(parse(text=y_lab)) +
            theme_classic() +
            theme(axis.text.x = element_text(size = 12),
                  axis.text.y = element_text(size = 12),
                  axis.title.x = element_text(size = 16),
                  axis.title.y = element_text(size = 16))
    }else if(type=='collage'){
        args <- unique(c(...))
        if(length(args)==0){
            args <- c('rating_curve','residuals','f','sigma_eps')
        }
        legal_args <- c('rating_curve','rating_curve_mean','sigma_eps','f','beta','residuals')
        if(all(args %in% legal_args)){
            plot_list <- lapply(args,function(y){
                plot(x,type=y,transformed=transformed)
            })
            p <- do.call(gridExtra::grid.arrange,c(plot_list,ncol=round(sqrt(length(args)))))
        }else{
            stop(cat(paste('For type collage, arguments must be in the following list:\n -',paste(legal_args,collapse='\n -'))))
        }
    }
    if(type!='collage' & !is.null(title)){
        p <- p + ggtitle(title) + theme(plot.title=element_text(size=18))
    }
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

#' Print bplm0 object
#'
#' Print the results of a bplm0 object
#' @param x an object of class "bplm0"
#' @seealso \code{\link{bplm0}} for fitting the bplm0 model, \code{\link{summary.bplm0}} for summaries, \code{\link{predict.bplm0}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{bplm0.plot}} to help visualize the full posterior distributions.
#' @examples
#' data(V316_river)
#' f <- Q~W
#' bplm.fit <- bplm(f,V316_river)
#' print(bplm0.fit)
#' @export
print.bplm0 <- function(x,...){
    print_fun(x)
}

#' Summarizing bplm0 fit
#'
#' Print the results of a  object
#' @param x an object of class "bplm0"
#' @seealso \code{\link{bplm0}} for fitting the bplm0 model, \code{\link{predict.bplm0}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{bplm0.plot}} to help visualize the full posterior distributions.
#' @examples
#' data(V316_river)
#' f <- Q~W
#' bplm0.fit <- bplm0(f,V316_river)
#' summary(bplm0.fit)
#' @export
summary.bplm0 <- function(object,...){
    summary_fun(object)
}

#' Plot bplm0 fit
#'
#' Print the results of a  object
#' @param x an object of class "bplm0"
#' @param type a character to denote the plot type. One of 'trace','histogram','rating_curve','rating_curve_mean','residuals'
#' @return bplm0 returns an object of class ggplot2
#' @seealso \code{\link{bplm0}} for fitting the bplm0 model,\code{\link{summary.bplm0}} for summaries, \code{\link{predict.bplm0}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{bplm0.plot}} to help visualize the full posterior distributions.
#' @examples
#' data(V316_river)
#' f <- Q~W
#' bplm0.fit <- bplm0(f,V316_river)
#' plot(bplm0.fit)
#' @export
#' @import ggplot2
#' @importFrom latex2exp TeX
#' @importFrom gridExtra grid.arrange
plot.bplm0 <- function(x,type='rating_curve',...,transformed=F,title=NULL){
    plot_fun(x,type,...,transformed=transformed,title=title)
}

#' Predict method for bplm0 fit
#'
#' Print the results of a  object
#' @param object an object of class "bplm0"
#' @param newdata a numeric vector of stage values for which to predict. If omitted, the stage values in the data are used.
#' @return numeric vector of discharge values for the stage values given in newdata
#' @seealso \code{\link{bplm0}} for fitting the bplm0 model,\code{\link{summary.bplm0}} for summaries, \code{\link{predict.bplm0}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{bplm0.plot}} to help visualize the full posterior distributions.
#' @examples
#' data(V316_river)
#' f <- Q~W
#' bplm0.fit <- bplm0(f,V316_river)
#' plot(bplm0.fit)
#' @export
predict.bplm0 <- function(object,newdata=NULL){
    predict_fun(object,newdata)
}

#' Print bplm object
#'
#' Print the results of a bplm object
#' @param x an object of class "bplm"
#' @return bgplm0 returns an object of class "bplm"\cr\cr
#' @seealso \code{\link{summary.bplm}} for summaries, \code{\link{predict.bplm}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{bplm.plot}} to help visualize the full posterior distributions.
#' @examples
#' data(V316_river)
#' f <- Q~W
#' bplm.fit <- bplm(f,V316_river)
#' print(bplm.fit)
#' @export
#'
print.bplm <- function(x,...){
    print_fun(x)
}

#' Summarizing bplm fit
#'
#' Print the results of a  object
#' @param x an object of class "bplm0"
#' @return bgplm0 returns an object of class "bplm"\cr\cr
#' @seealso \code{\link{bplm0}} for fitting the bplm0 model,\code{\link{summary.bplm0}} for summaries, \code{\link{predict.bplm0}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{bplm0.plot}} to help visualize the full posterior distributions.
#' @examples
#' data(V316_river)
#' f <- Q~W
#' bplm0.fit <- bplm0(f,V316_river)
#' summary(bplm0.fit)
#' @export
summary.bplm <- function(object,...){
    summary_fun(object)
}

#' Plot bplm0 fit
#'
#' Print the results of a  object
#' @param x an object of class "bplm0"
#' @param type a character to denote the plot type. One of 'trace','histogram','rating_curve','rating_curve_mean','residuals'
#' @return bplm0 returns an object of class "bplm0"\cr\cr
#' @seealso \code{\link{bplm0}} for fitting the bplm0 model,\code{\link{summary.bplm0}} for summaries, \code{\link{predict.bplm0}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{bplm0.plot}} to help visualize the full posterior distributions.
#' @examples
#' data(V316_river)
#' f <- Q~W
#' bplm0.fit <- bplm0(f,V316_river)
#' plot(bplm0.fit)
#' @export
#' @import ggplot2
#' @importFrom latex2exp TeX
#' @importFrom gridExtra grid.arrange
plot.bplm <- function(x,type='rating_curve',...,transformed=F,title=NULL){
    plot_fun(x,type,...,transformed=transformed,title=title)
}

#' Predict method for bplm0 fit
#'
#' Print the results of a  object
#' @param object an object of class "bplm0"
#' @param newdata a numeric vector of stage values for which to predict. If omitted, the stage values in the data are used.
#' @return numeric vector of discharge values for the stage values given in newdata
#' @seealso \code{\link{bplm0}} for fitting the bplm0 model,\code{\link{summary.bplm0}} for summaries, \code{\link{predict.bplm0}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{bplm0.plot}} to help visualize the full posterior distributions.
#' @examples
#' data(V316_river)
#' f <- Q~W
#' bplm0.fit <- bplm0(f,V316_river)
#' plot(bplm0.fit)
#' @export
predict.bplm <- function(object,newdata=NULL){
    predict_fun(object,newdata)
}

#' Print bgplm0 object
#'
#' Print the results of a bgplm0 object
#' @param x an object of class "bgplm0"
#' @seealso \code{\link{summary.bgplm0}} for summaries, \code{\link{predict.bgplm0}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{bgplm0.plot}} to help visualize the full posterior distributions.
#' @examples
#' data(V316_river)
#' f <- Q~W
#' bgplm0.fit <- bgplm0(f,V316_river)
#' print(bgplm0.fit)
#' @export
print.bgplm0 <- function(x,...){
    print_fun(x)
}

#' Summarizing bplm fit
#'
#' Print the results of a  object
#' @param x an object of class "bplm0"
#' @return bgplm0 returns an object of class "bplm"\cr\cr
#' @seealso \code{\link{bplm0}} for fitting the bplm0 model,\code{\link{summary.bplm0}} for summaries, \code{\link{predict.bplm0}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{bplm0.plot}} to help visualize the full posterior distributions.
#' @examples
#' data(V316_river)
#' f <- Q~W
#' bplm0.fit <- bplm0(f,V316_river)
#' summary(bplm0.fit)
#' @export
summary.bgplm0 <- function(object,...){
    summary_fun(object)
}

#' Plot bplm0 fit
#'
#' Print the results of a  object
#' @param x an object of class "bplm0"
#' @param type a character to denote the plot type. One of 'trace','histogram','rating_curve','rating_curve_mean','residuals'
#' @return bplm0 returns an object of class "bplm0"\cr\cr
#' @seealso \code{\link{bplm0}} for fitting the bplm0 model,\code{\link{summary.bplm0}} for summaries, \code{\link{predict.bplm0}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{bplm0.plot}} to help visualize the full posterior distributions.
#' @examples
#' data(V316_river)
#' f <- Q~W
#' bplm0.fit <- bplm0(f,V316_river)
#' plot(bplm0.fit)
#' @export
#' @import ggplot2
#' @importFrom latex2exp TeX
#' @importFrom gridExtra grid.arrange
plot.bgplm0 <- function(x,type='rating_curve',...,transformed=F,title=NULL){
    plot_fun(x,type,...,transformed=transformed,title=title)
}

#' Predict method for bplm0 fit
#'
#' Print the results of a  object
#' @param object an object of class "bplm0"
#' @param newdata a numeric vector of stage values for which to predict. If omitted, the stage values in the data are used.
#' @return numeric vector of discharge values for the stage values given in newdata
#' @seealso \code{\link{bplm0}} for fitting the bplm0 model,\code{\link{summary.bplm0}} for summaries, \code{\link{predict.bplm0}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{bplm0.plot}} to help visualize the full posterior distributions.
#' @examples
#' data(V316_river)
#' f <- Q~W
#' bplm0.fit <- bplm0(f,V316_river)
#' plot(bplm0.fit)
#' @export
predict.bgplm0 <- function(object,newdata=NULL){
    predict_fun(object,newdata)
}

#' Print bgplm object
#'
#' Print the results of a bgplm object
#' @param x an object of class "bgplm"
#' @seealso \code{\link{summary.bgplm}} for summaries, \code{\link{predict.bgplm}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{bgplm.plot}} to help visualize the full posterior distributions.
#' @examples
#' data(V316_river)
#' f <- Q~W
#' bgplm.fit <- bgplm(f,V316_river)
#' print(bgplm.fit)
#' @export
print.bgplm <- function(x,...){
    print_fun(x)
}

#' Summarizing bplm fit
#'
#' Print the results of a  object
#' @param x an object of class "bplm0"
#' @return bgplm0 returns an object of class "bplm"\cr\cr
#' @seealso \code{\link{bplm0}} for fitting the bplm0 model,\code{\link{summary.bplm0}} for summaries, \code{\link{predict.bplm0}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{bplm0.plot}} to help visualize the full posterior distributions.
#' @examples
#' data(V316_river)
#' f <- Q~W
#' bplm0.fit <- bplm0(f,V316_river)
#' summary(bplm0.fit)
#' @export
summary.bgplm <- function(object,...){
    summary_fun(object)
}

#' Plot bplm0 fit
#'
#' Print the results of a  object
#' @param x an object of class "bplm0"
#' @param type a character to denote the plot type. One of 'trace','histogram','rating_curve','rating_curve_mean','residuals'
#' @return bplm0 returns an object of class "bplm0"\cr\cr
#' @seealso \code{\link{bplm0}} for fitting the bplm0 model,\code{\link{summary.bplm0}} for summaries, \code{\link{predict.bplm0}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{bplm0.plot}} to help visualize the full posterior distributions.
#' @examples
#' data(V316_river)
#' f <- Q~W
#' bplm0.fit <- bplm0(f,V316_river)
#' plot(bplm0.fit)
#' @export
#' @import ggplot2
#' @importFrom latex2exp TeX
#' @importFrom gridExtra grid.arrange
plot.bgplm <- function(x,type='rating_curve',...,transformed=F,title=NULL){
    plot_fun(x,type,...,transformed=transformed,title=title)
}

#' Predict method for bplm0 fit
#'
#' Print the results of a  object
#' @param object an object of class "bplm0"
#' @param newdata a numeric vector of stage values for which to predict. If omitted, the stage values in the data are used.
#' @return numeric vector of discharge values for the stage values given in newdata
#' @seealso \code{\link{bplm0}} for fitting the bplm0 model,\code{\link{summary.bplm0}} for summaries, \code{\link{predict.bplm0}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{bplm0.plot}} to help visualize the full posterior distributions.
#' @examples
#' data(V316_river)
#' f <- Q~W
#' bplm0.fit <- bplm0(f,V316_river)
#' plot(bplm0.fit)
#' @export
predict.bgplm <- function(object,newdata=NULL){
    predict_fun(object,newdata)
}

# gather_draws(bgplm.fit, 'rating_curve') %>%
#     group_by(h) %>%
#     summarise(lower_50 = quantile(value, 0.25),
#               upper_50 = quantile(value, 0.75),
#               lower_60 = quantile(value, 0.2),
#               upper_60 = quantile(value, 0.8),
#               lower_70 = quantile(value, 0.15),
#               upper_70 = quantile(value, 0.85),
#               lower_80 = quantile(value, 0.1),
#               upper_80 = quantile(value, 0.9),
#               lower_90 = quantile(value, 0.05),
#               upper_90 = quantile(value, 0.95),
#               lower_95 = quantile(value, 0.025),
#               upper_95 = quantile(value, 0.975)) %>%
#     pivot_longer(c(-h), names_to = c("which", "prob"), names_sep = "_") %>%
#     pivot_wider(names_from = which, values_from = value) %>%
#     mutate(prob = as.numeric(prob)) %>%
#     ggplot() +
#     geom_ribbon(aes(h, ymin = lower, ymax = upper,fill = factor(-prob)), alpha = 0.7) +
#     geom_point(data=bgplm.fit$data,aes(W,Q))+
#     scale_fill_brewer() +
#     scale_x_continuous(expand = c(0,0))+
#     scale_y_continuous(expand = c(0,0)) +
#     theme_classic() +
#     theme(legend.position = "none")
