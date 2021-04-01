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
#' Visualize results from model ojbects in bdrc, bplm0, bplm, bgplm0,bgplm
#' @param x an object of class "bplm0","bplm","bgplm0" or "bgplm"
#' @param ... Arguments to be passed to other methods. The following arguments are supported:
#' \itemize{
#'  \item{\code{type}}{ a character denoting what type of plot should be drawn. Possible types are
#'                    \itemize{
#'                       \item{"rating_curve"}{ to plot the rating curve on original scale.}
#'                       \item{"rating_curve_mean"}{ to plot the rating curve on a log scale.}
#'                       \item{"f"}{ to plot the power-law exponent}
#'                       \item{"beta"}{ to plot the random effect in the power-law exponent}
#'                       \item{"sigma_eps"}{ to plot the standard deviation on the data level}
#'                       \item{"residuals"}{ to plot the log residuals}
#'                       \item{"trace"}{ to plot trace plots of a list of one or more parameter. Requires parameter names as additional arguments, see below.}
#'                       \item{"histogram"}{ to plot histograms of a list of one or more parameter. Requires parameter names as additional arguments, see below.}
#'                    }}
#' \item{\code{transformed}}{ a logical value indicating whether the quantity should be plotted on a transformed scale used during the Bayesian inference. Defaults to FALSE.}
#' \item{\code{title}}{ title of the plot. Defaults to NULL, i.e. no title}
#' \item{}{ additional characters or character vectors denoting the parameters to plot. Used when type is "trace" or "histogram". Allowed values are the parameter names found in the summary of the model object. See Examples.}
#' }
#' @return returns an object of class ggplot2 or Grob object
#' @importFrom ggplot2 ggplot aes geom_point geom_path geom_histogram geom_abline facet_wrap scale_color_manual scale_x_continuous scale_y_continuous label_parsed ggtitle xlab ylab
#' @importFrom rlang .data
#' @importFrom stats median
plot_fun <- function(x,...){
    args <- list(...)
    if('type' %in% names(args)){
        type <- args$type
        args <- args[names(args)!='type'] # remove type from args list
    }else{
        type <- 'rating_curve'
    }
    if('transformed' %in% names(args)){
        transformed <- args$transformed
        args <- args[names(args)!='transformed'] # remove transformed from args list
    }else{
        transformed <- F
    }
    if('title' %in% names(args)){
        title <- args$title
        args <- args[names(args)!='title'] # remove title from args list
    }else{
        title <- NULL
    }
    args <- unlist(args)
    legal_types <- c('rating_curve','rating_curve_mean','f','beta','sigma_eps','residuals','trace','histogram')
    if(!(type %in% legal_types)){
        stop(cat(paste('Type argument not recognized. Possible types are:\n -',paste(legal_types,collapse='\n - '))))
    }else if(type=='trace'){
        plot_dat <- gather_draws(x,args,transformed=transformed)
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
        plot_dat <- gather_draws(x,args,transformed=transformed)
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
            plot_dat$log_Q <- log(plot_dat$Q)
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
            stop('Plots of type "beta" are only for models with stage dependent power law exponent, s.a. "bgplm0" and "bgplm"')
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
plot_collage <- function(x,...){
    args <- list(...)
    args <- args[names(args)!='type']
    if('transformed' %in% names(args)){
        transformed <- args$transformed
        args <- args[names(args)!='transformed'] # remove transformed from args list
    }else{
        transformed <- F
    }
    args <- unique(unlist(args))
    if(length(args)==0){
        args <- c('rating_curve','residuals','f','sigma_eps')
    }
    legal_args <- c('rating_curve','rating_curve_mean','sigma_eps','f','beta','residuals')
    if(all(args %in% legal_args)){
        plot_list <- lapply(args,function(y){
            plot_fun(x,type=y,transformed=transformed)
        })
        p <- do.call(arrangeGrob,c(plot_list,ncol=round(sqrt(length(args)))))
    }else{
        stop(cat(paste('For type collage, arguments must be in the following list:\n -',paste(legal_args,collapse='\n - '))))
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
#' @param ... not used in this function
#' @seealso \code{\link{bplm0}} for fitting the bplm0 model, \code{\link{summary.bplm0}} for summaries, \code{\link{predict.bplm0}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{plot.bplm0}} to help visualize the full posterior distributions.
#' @examples
#' \dontrun{
#' data(V316_river)
#' f <- Q~W
#' bplm0.fit <- bplm0(f,V316_river)
#' print(bplm0.fit)
#' }
#' @export
print.bplm0 <- function(x,...){
    print_fun(x)
}

#' Summarizing bplm0 fit
#'
#' Summarize the results of a bplm0 object
#' @param object an object of class "bplm0"
#' @param ... Not used for this function
#' @seealso \code{\link{bplm0}} for fitting the bplm0 model, \code{\link{predict.bplm0}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{plot.bplm0}} to help visualize the full posterior distributions.
#' @examples
#' \dontrun{
#' data(V316_river)
#' f <- Q~W
#' bplm0.fit <- bplm0(f,V316_river)
#' summary(bplm0.fit)
#' }
#' @export
summary.bplm0 <- function(object,...){
    summary_fun(object)
}

#' Autoplot bplm0 fit
#'
#' Uses ggplot2 to plot bplm0 object
#' @param x an object of class "bplm0"
#' @param ... Arguments to be passed to other methods. The following arguments are supported:
#' \itemize{
#'  \item{\code{type}}{ a character denoting what type of plot should be drawn. Possible types are
#'                    \itemize{
#'                       \item{"rating_curve"}{ to plot the rating curve on original scale.}
#'                       \item{"rating_curve_mean"}{ to plot the rating curve on a log scale.}
#'                       \item{"f"}{ to plot the power-law exponent}
#'                       \item{"sigma_eps"}{ to plot the standard deviation on the data level}
#'                       \item{"residuals"}{ to plot the log residuals}
#'                       \item{"trace"}{ to plot trace plots of a list of one or more parameter. Requires parameter names as additional arguments, see below.}
#'                       \item{"histogram"}{ to plot histograms of a list of one or more parameter. Requires parameter names as additional arguments, see below.}
#'                    }}
#' \item{\code{transformed}}{ a logical value indicating whether the quantity should be plotted on a transformed scale used during the Bayesian inference. Defaults to FALSE.}
#' \item{\code{title}}{ title of the plot. Defaults to NULL, i.e. no title}
#' \item{additional}{ characters or character vectors denoting the parameters to plot. Used when type is "trace" or "histogram". Allowed values are the parameter names found in the summary of the model object. See Examples.}
#' }
#' @return returns an object of class ggplot2
#' @seealso \code{\link{bplm0}} for fitting the bplm0 model,\code{\link{summary.bplm0}} for summaries of model parameters, \code{\link{predict.bplm0}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{plot.bplm0}} to help visualize the full posterior distributions.
#' @examples
#' \dontrun{
#' data(V316_river)
#' f <- Q~W
#' bplm0.fit <- bplm0(f,V316_river)
#' autoplot(bplm0.fit)
#' }
#' @export
autoplot.bplm0 <- function(x,...){
    plot_fun(x,...)
}

#' Plot bplm0 fit
#'
#' Print the results of a  object
#' @param x an object of class "bplm0"
#' @param ... Arguments to be passed to other methods. The following arguments are supported:
#' \itemize{
#'  \item{\code{type}}{ a character denoting what type of plot should be drawn. Possible types are
#'                    \itemize{
#'                       \item{"rating_curve"}{ to plot the rating curve on original scale.}
#'                       \item{"rating_curve_mean"}{ to plot the rating curve on a log scale.}
#'                       \item{"f"}{ to plot the power-law exponent}
#'                       \item{"sigma_eps"}{ to plot the standard deviation on the data level}
#'                       \item{"residuals"}{ to plot the log residuals}
#'                       \item{"trace"}{ to plot trace plots of a list of one or more parameter. Requires parameter names as additional arguments, see below.}
#'                       \item{"histogram"}{ to plot histograms of a list of one or more parameter. Requires parameter names as additional arguments, see below.}
#'                    }}
#' \item{\code{transformed}}{ a logical value indicating whether the quantity should be plotted on a transformed scale used during the Bayesian inference. Defaults to FALSE.}
#' \item{\code{title}}{ title of the plot. Defaults to NULL, i.e. no title}
#' \item{additional}{ characters or character vectors denoting the parameters to plot. Used when type is "trace" or "histogram". Allowed values are the parameter names found in the summary of the model object. See Examples.}
#' }
#' @seealso \code{\link{bplm0}} for fitting the bplm0 model,\code{\link{summary.bplm0}} for summaries, \code{\link{predict.bplm0}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{plot.bplm0}} to help visualize the full posterior distributions.
#' @examples
#' \dontrun{
#' data(V316_river)
#' f <- Q~W
#' bplm0.fit <- bplm0(f,V316_river)
#' plot(bplm0.fit)
#' }
#' @export
#' @importFrom grid grid.draw
#' @importFrom ggplot2 autoplot
plot.bplm0 <- function(x,...){
    args <- list(...)
    if(is.null(args$type) || args$type!='collage'){
        p <- autoplot(x,...)
        print(p)
    }else{
        p <- plot_collage(x,...)
        grid.draw(p)
    }
}

#' Predict method for bplm0 fit
#'
#' Print the results of a  object
#' @param object an object of class "bplm0"
#' @param newdata a numeric vector of stage values for which to predict. If omitted, the stage values in the data are used.
#' @param ... not used in this function
#' @return numeric vector of discharge values for the stage values given in newdata
#' @seealso \code{\link{bplm0}} for fitting the bplm0 model,\code{\link{summary.bplm0}} for summaries, \code{\link{predict.bplm0}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{plot.bplm0}} to help visualize the full posterior distributions.
#' @examples
#' \dontrun{
#' data(V316_river)
#' f <- Q~W
#' bplm0.fit <- bplm0(f,V316_river,h_max=2)
#' #predict rating curve on a equally 1 cm spaced grid from 1 to 2 meters
#' predict(bplm0.fit,newdata=seq(1,2,by=0.01))
#' }
#' @export
predict.bplm0 <- function(object,newdata=NULL,...){
    predict_fun(object,newdata)
}

#' Print bplm object
#'
#' Print the results of a bplm object
#' @param x an object of class "bplm"
#' @param ... not used in this function
#' @return bgplm0 returns an object of class "bplm"\cr\cr
#' @seealso \code{\link{summary.bplm}} for summaries, \code{\link{predict.bplm}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{plot.bplm}} to help visualize the full posterior distributions.
#' @examples
#' \dontrun{
#' data(V316_river)
#' f <- Q~W
#' bplm.fit <- bplm(f,V316_river)
#' print(bplm.fit)
#' }
#' @export
#'
print.bplm <- function(x,...){
    print_fun(x)
}

#' Summarizing bplm fit
#'
#' Summarize the results of a bplm object
#' @param object an object of class "bplm"
#' @param ... not used in this function
#' @return bgplm0 returns an object of class "bplm"\cr\cr
#' @seealso \code{\link{bplm}} for fitting the bplm model,\code{\link{summary.bplm}} for summaries, \code{\link{predict.bplm}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{plot.bplm}} to help visualize the full posterior distributions.
#' @examples
#' \dontrun{
#' data(V316_river)
#' f <- Q~W
#' bplm.fit <- bplm(f,V316_river)
#' summary(bplm.fit)
#' }
#' @export
summary.bplm <- function(object,...){
    summary_fun(object)
}

#' Autoplot bplm fit
#'
#' Uses ggplot2 to plot bplm object
#' @param x an object of class "bplm"
#' @param ... Arguments to be passed to other methods. The following arguments are supported:
#' \itemize{
#'  \item{\code{type}}{ a character denoting what type of plot should be drawn. Possible types are
#'                    \itemize{
#'                       \item{"rating_curve"}{ to plot the rating curve on original scale.}
#'                       \item{"rating_curve_mean"}{ to plot the rating curve on a log scale.}
#'                       \item{"f"}{ to plot the power-law exponent}
#'                       \item{"sigma_eps"}{ to plot the standard deviation on the data level}
#'                       \item{"residuals"}{ to plot the log residuals}
#'                       \item{"trace"}{ to plot trace plots of a list of one or more parameter. Requires parameter names as additional arguments, see below.}
#'                       \item{"histogram"}{ to plot histograms of a list of one or more parameter. Requires parameter names as additional arguments, see below.}
#'                    }}
#' \item{\code{transformed}}{ a logical value indicating whether the quantity should be plotted on a transformed scale used during the Bayesian inference. Defaults to FALSE.}
#' \item{\code{title}}{ title of the plot. Defaults to NULL, i.e. no title}
#' \item{additional}{ characters or character vectors denoting the parameters to plot. Used when type is "trace" or "histogram". Allowed values are the parameter names found in the summary of the model object. See Examples.}
#' }
#' @return returns an object of class ggplot2
#' @seealso \code{\link{bplm}} for fitting the bplm model,\code{\link{summary.bplm}} for summaries of model parameters, \code{\link{predict.bplm}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{plot.bplm}} to help visualize the full posterior distributions.
#' @examples
#' \dontrun{
#' data(V316_river)
#' f <- Q~W
#' bplm.fit <- bplm(f,V316_river)
#' autoplot(bplm.fit)
#' }
#' @export
autoplot.bplm <- function(x,...){
    plot_fun(x,...)
}

#' Plot bplm fit
#'
#' Print the results of a  object
#' @param x an object of class "bplm"
#' @param ... Arguments to be passed to other methods. The following arguments are supported:
#' \itemize{
#'  \item{\code{type}}{ a character denoting what type of plot should be drawn. Possible types are
#'                    \itemize{
#'                       \item{"rating_curve"}{ to plot the rating curve on original scale.}
#'                       \item{"rating_curve_mean"}{ to plot the rating curve on a log scale.}
#'                       \item{"f"}{ to plot the power-law exponent}
#'                       \item{"sigma_eps"}{ to plot the standard deviation on the data level}
#'                       \item{"residuals"}{ to plot the log residuals}
#'                       \item{"trace"}{ to plot trace plots of a list of one or more parameter. Requires parameter names as additional arguments, see below.}
#'                       \item{"histogram"}{ to plot histograms of a list of one or more parameter. Requires parameter names as additional arguments, see below.}
#'                    }}
#' \item{\code{transformed}}{ a logical value indicating whether the quantity should be plotted on a transformed scale used during the Bayesian inference. Defaults to FALSE.}
#' \item{\code{title}}{ title of the plot. Defaults to NULL, i.e. no title}
#' \item{additional}{ characters or character vectors denoting the parameters to plot. Used when type is "trace" or "histogram". Allowed values are the parameter names found in the summary of the model object. See Examples.}
#' }
#' @seealso \code{\link{bplm}} for fitting the bplm model,\code{\link{summary.bplm}} for summaries, \code{\link{predict.bplm}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{plot.bplm}} to help visualize the full posterior distributions.
#' @examples
#' \dontrun{
#' data(V316_river)
#' f <- Q~W
#' bplm.fit <- bplm(f,V316_river)
#' plot(bplm.fit)
#' }
#' @export
#' @importFrom grid grid.draw
#' @importFrom ggplot2 autoplot
plot.bplm <- function(x,...){
    args <- list(...)
    if(is.null(args$type) || args$type!='collage'){
        p <- autoplot(x,...)
        print(p)
    }else{
        p <- plot_collage(x,...)
        grid.draw(p)
    }
}

#' Predict method for bplm fit
#'
#' Print the results of a  object
#' @param object an object of class "bplm"
#' @param newdata a numeric vector of stage values for which to predict. If omitted, the stage values in the data are used.
#' @param ... not used in this function
#' @return numeric vector of discharge values for the stage values given in newdata
#' @seealso \code{\link{bplm}} for fitting the bplm model,\code{\link{summary.bplm}} for summaries, \code{\link{predict.bplm}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{plot.bplm}} to help visualize the full posterior distributions.
#' @examples
#' \dontrun{
#' data(V316_river)
#' f <- Q~W
#' bplm.fit <- bplm(f,V316_river,h_max=2)
#' #predict rating curve on a equally 1 cm spaced grid from 1 to 2 meters
#' predict(bplm.fit,newdata=seq(1,2,by=0.01))
#' }
#' @export
predict.bplm <- function(object,newdata=NULL,...){
    predict_fun(object,newdata)
}

#' Print bgplm0 object
#'
#' Print the results of a bgplm0 object
#' @param x an object of class "bgplm0"
#' @param ... not used in this function
#' @seealso \code{\link{summary.bgplm0}} for summaries, \code{\link{predict.bgplm0}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{plot.bgplm0}} to help visualize the full posterior distributions.
#' @examples
#' \dontrun{
#' data(V316_river)
#' f <- Q~W
#' bgplm0.fit <- bgplm0(f,V316_river)
#' print(bgplm0.fit)
#' }
#' @export
print.bgplm0 <- function(x,...){
    print_fun(x)
}

#' Summarizing bgplm0 fit
#'
#' Summarize the results of a bgplm0 object
#' @param object an object of class "bgplm0"
#' @param ... not used in this function
#' @return bgplm0 returns an object of class "bplm"\cr\cr
#' @seealso \code{\link{bgplm0}} for fitting the bgplm0 model,\code{\link{summary.bgplm0}} for summaries, \code{\link{predict.bgplm0}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{plot.bgplm0}} to help visualize the full posterior distributions.
#' @examples
#' \dontrun{
#' data(V316_river)
#' f <- Q~W
#' bgplm0.fit <- bgplm0(f,V316_river)
#' summary(bgplm0.fit)
#' }
#' @export
summary.bgplm0 <- function(object,...){
    summary_fun(object)
}

#' Autoplot bgplm0 fit
#'
#' Uses ggplot2 to plot bgplm0 object
#' @param x an object of class "bgplm0"
#' @param ... Arguments to be passed to other methods. The following arguments are supported:
#' \itemize{
#'  \item{\code{type}}{ a character denoting what type of plot should be drawn. Possible types are
#'                    \itemize{
#'                       \item{"rating_curve"}{ to plot the rating curve on original scale.}
#'                       \item{"rating_curve_mean"}{ to plot the rating curve on a log scale.}
#'                       \item{"f"}{ to plot the power-law exponent}
#'                       \item{"sigma_eps"}{ to plot the standard deviation on the data level}
#'                       \item{"residuals"}{ to plot the log residuals}
#'                       \item{"trace"}{ to plot trace plots of a list of one or more parameter. Requires parameter names as additional arguments, see below.}
#'                       \item{"histogram"}{ to plot histograms of a list of one or more parameter. Requires parameter names as additional arguments, see below.}
#'                    }}
#' \item{\code{transformed}}{ a logical value indicating whether the quantity should be plotted on a transformed scale used during the Bayesian inference. Defaults to FALSE.}
#' \item{\code{title}}{ title of the plot. Defaults to NULL, i.e. no title}
#' \item{additional}{ characters or character vectors denoting the parameters to plot. Used when type is "trace" or "histogram". Allowed values are the parameter names found in the summary of the model object. See Examples.}
#' }
#' @return returns an object of class ggplot2
#' @seealso \code{\link{bgplm0}} for fitting the bgplm0 model,\code{\link{summary.bgplm0}} for summaries of model parameters, \code{\link{predict.bgplm0}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{plot.bgplm0}} to help visualize the full posterior distributions.
#' @examples
#' \dontrun{
#' data(V316_river)
#' f <- Q~W
#' bgplm0.fit <- bgplm0(f,V316_river)
#' autoplot(bgplm0.fit)
#' }
#' @export
autoplot.bgplm0 <- function(x,...){
    plot_fun(x,...)
}

#' Plot bgplm0 fit
#'
#' Print the results of a  object
#' @param x an object of class "bgplm0"
#' @param ... Arguments to be passed to other methods. The following arguments are supported:
#' \itemize{
#'  \item{\code{type}}{ a character denoting what type of plot should be drawn. Possible types are
#'                    \itemize{
#'                       \item{"rating_curve"}{ to plot the rating curve on original scale.}
#'                       \item{"rating_curve_mean"}{ to plot the rating curve on a log scale.}
#'                       \item{"f"}{ to plot the power-law exponent}
#'                       \item{"sigma_eps"}{ to plot the standard deviation on the data level}
#'                       \item{"residuals"}{ to plot the log residuals}
#'                       \item{"trace"}{ to plot trace plots of a list of one or more parameter. Requires parameter names as additional arguments, see below.}
#'                       \item{"histogram"}{ to plot histograms of a list of one or more parameter. Requires parameter names as additional arguments, see below.}
#'                    }}
#' \item{\code{transformed}}{ a logical value indicating whether the quantity should be plotted on a transformed scale used during the Bayesian inference. Defaults to FALSE.}
#' \item{\code{title}}{ title of the plot. Defaults to NULL, i.e. no title}
#' \item{additional}{ characters or character vectors denoting the parameters to plot. Used when type is "trace" or "histogram". Allowed values are the parameter names found in the summary of the model object. See Examples.}
#' }
#' @seealso \code{\link{bgplm0}} for fitting the bgplm0 model,\code{\link{summary.bgplm0}} for summaries, \code{\link{predict.bgplm0}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{plot.bgplm0}} to help visualize the full posterior distributions.
#' @examples
#' \dontrun{
#' data(V316_river)
#' f <- Q~W
#' bgplm0.fit <- bgplm0(f,V316_river)
#' plot(bgplm0.fit)
#' }
#' @export
#' @importFrom grid grid.draw
#' @importFrom ggplot2 autoplot
plot.bgplm0 <- function(x,...){
    args <- list(...)
    if(is.null(args$type) || args$type!='collage'){
        p <- autoplot(x,...)
        print(p)
    }else{
        p <- plot_collage(x,...)
        grid.draw(p)
    }
}

#' Predict method for bgplm0 fit
#'
#' Print the results of a  object
#' @param object an object of class "bgplm0"
#' @param newdata a numeric vector of stage values for which to predict. If omitted, the stage values in the data are used.
#' @param ... not used in this function
#' @return numeric vector of discharge values for the stage values given in newdata
#' @seealso \code{\link{bgplm0}} for fitting the bgplm0 model,\code{\link{summary.bgplm0}} for summaries, \code{\link{predict.bgplm0}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{plot.bgplm0}} to help visualize the full posterior distributions.
#' @examples
#' \dontrun{
#' data(V316_river)
#' f <- Q~W
#' bgplm0.fit <- bgplm0(f,V316_river,h_max=2)
#' #predict rating curve on a equally 1 cm spaced grid from 1 to 2 meters
#' predict(bgplm0.fit,newdata=seq(1,2,by=0.01))
#' }
#' @export
predict.bgplm0 <- function(object,newdata=NULL,...){
    predict_fun(object,newdata)
}

#' Print bgplm object
#'
#' Print the results of a bgplm object
#' @param x an object of class "bgplm"
#' @param ... not used in this function
#' @seealso \code{\link{summary.bgplm}} for summaries, \code{\link{predict.bgplm}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{plot.bgplm}} to help visualize the full posterior distributions.
#' @examples
#' \dontrun{
#' data(V316_river)
#' f <- Q~W
#' bgplm.fit <- bgplm(f,V316_river)
#' print(bgplm.fit)
#' }
#' @export
print.bgplm <- function(x,...){
    print_fun(x)
}

#' Summarizing bplm fit
#'
#' Summarize the results of a bgplm object
#' @param object an object of class "bgplm"
#' @param ... not used in this function
#' @return bgplm0 returns an object of class "bplm"\cr\cr
#' @seealso \code{\link{bgplm}} for fitting the bgplm model,\code{\link{summary.bgplm}} for summaries, \code{\link{predict.bgplm}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{plot.bgplm}} to help visualize the full posterior distributions.
#' @examples
#' \dontrun{
#' data(V316_river)
#' f <- Q~W
#' bgplm.fit <- bgplm(f,V316_river)
#' summary(bgplm.fit)
#' }
#' @export
summary.bgplm <- function(object,...){
    summary_fun(object)
}

#' Autoplot bgplm fit
#'
#' Uses ggplot2 to plot bgplm object
#'
#' @param x an object of class "bgplm"
#' @param ... Arguments to be passed to other methods. The following arguments are supported:
#' \itemize{
#'   \item{\code{type}}{ a character denoting what type of plot should be drawn. Possible types are
#'                    \itemize{
#'                       \item{"rating_curve"}{ to plot the rating curve on original scale.}
#'                       \item{"rating_curve_mean"}{ to plot the rating curve on a log scale.}
#'                       \item{"f"}{ to plot the power-law exponent}
#'                       \item{"sigma_eps"}{ to plot the standard deviation on the data level}
#'                       \item{"residuals"}{ to plot the log residuals}
#'                       \item{"trace"}{ to plot trace plots of a list of one or more parameter. Requires parameter names as additional arguments, see below.}
#'                       \item{"histogram"}{ to plot histograms of a list of one or more parameter. Requires parameter names as additional arguments, see below.}
#'                    }}
#'   \item{\code{transformed}}{ a logical value indicating whether the quantity should be plotted on a transformed scale used during the Bayesian inference. Defaults to FALSE.}
#'   \item{\code{title}}{ title of the plot. Defaults to NULL, i.e. no title}
#'   \item{additional}{ characters or character vectors denoting the parameters to plot. Used when type is "trace" or "histogram". Allowed values are the parameter names found in the summary of the model object. See Examples.}
#' }
#' @return returns an object of class ggplot2
#' @seealso \code{\link{bgplm}} for fitting the bgplm model,\code{\link{summary.bgplm}} for summaries of model parameters, \code{\link{predict.bgplm}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{plot.bgplm}} to help visualize the full posterior distributions.
#' @examples
#' \dontrun{
#' x=1
#' }
#' @export
autoplot.bgplm <- function(x,...){
    plot_fun(x,...)
}

#' Plot bgplm fit
#'
#' Print the results of a  object
#'
#' @param x an object of class "bgplm"
#' @param ... Arguments to be passed to other methods. The following arguments are supported:
#' \itemize{
#'  \item{\code{type}}{ a character denoting what type of plot should be drawn. Possible types are
#'                    \itemize{
#'                       \item{"rating_curve"}{ to plot the rating curve on original scale.}
#'                       \item{"rating_curve_mean"}{ to plot the rating curve on a log scale.}
#'                       \item{"f"}{ to plot the power-law exponent}
#'                       \item{"sigma_eps"}{ to plot the standard deviation on the data level}
#'                       \item{"residuals"}{ to plot the log residuals}
#'                       \item{"trace"}{ to plot trace plots of a list of one or more parameter. Requires parameter names as additional arguments, see below.}
#'                       \item{"histogram"}{ to plot histograms of a list of one or more parameter. Requires parameter names as additional arguments, see below.}
#'                    }}
#' \item{\code{transformed}}{ a logical value indicating whether the quantity should be plotted on a transformed scale used during the Bayesian inference. Defaults to FALSE.}
#' \item{\code{title}}{ title of the plot. Defaults to NULL, i.e. no title}
#' \item{additional }{characters or character vectors denoting the parameters to plot. Used when type is "trace" or "histogram". Allowed values are the parameter names found in the summary of the model object. See Examples.}
#' }
#' @seealso \code{\link{bgplm}} for fitting the bgplm model,\code{\link{summary.bgplm}} for summaries, \code{\link{predict.bgplm}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{plot.bgplm}} to help visualize the full posterior distributions.
#' @examples
#' \dontrun{
#' data(V316_river)
#' f <- Q~W
#' bgplm.fit <- bgplm(f,V316_river)
#' plot(bgplm.fit)
#' }
#' @export
#' @importFrom grid grid.draw
#' @importFrom ggplot2 autoplot
plot.bgplm <- function(x,...){
    args <- list(...)
    if(is.null(args$type) || args$type!='collage'){
        p <- autoplot(x,...)
        print(p)
    }else{
        p <- plot_collage(x,...)
        grid.draw(p)
    }
}

#' Predict method for bgplm fit
#'
#' Print the results of a  object
#' @param object an object of class "bgplm"
#' @param newdata a numeric vector of stage values for which to predict. If omitted, the stage values in the data are used.
#' @param ... not used in this function
#' @return numeric vector of discharge values for the stage values given in newdata
#' @seealso \code{\link{bgplm}} for fitting the bgplm model,\code{\link{summary.bgplm}} for summaries, \code{\link{predict.bgplm}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{plot.bgplm}} to help visualize the full posterior distributions.
#' @examples
#' \dontrun{
#' data(V316_river)
#' f <- Q~W
#' bgplm0.fit <- bgplm0(f,V316_river,h_max=2)
#' #predict rating curve on a equally 1 cm spaced grid from 1 to 2 meters
#' predict(bgplm0.fit,newdata=seq(1,2,by=0.01))
#' }
#' @export
predict.bgplm <- function(object,newdata=NULL,...){
    predict_fun(object,newdata)
}
