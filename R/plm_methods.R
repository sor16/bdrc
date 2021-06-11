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
#' @param ... not used in this function
#' @param scaling a numerical value which can be used to scale up or down the size of the text and titles of a plot that uses \code{theme_bdrc}. Defaults to 1.
#' @return returns a theme object for the package
#' @importFrom ggplot2 %+replace% theme_classic theme element_text element_blank element_rect
theme_bdrc <- function(...,scaling=1){
    title_size <- scaling*16
    text_size <- scaling*12
    plot_title_size=scaling*18
    theme_classic() %+replace%
        theme( #text = element_text(family="Times", face="plain"),
               strip.background = element_blank(),
               strip.text.x = element_text(size = title_size),
               axis.title.x = element_text(size=title_size),
               axis.title.y = element_text(size=title_size,angle=90),
               axis.text.x = element_text(size=text_size),
               axis.text.y = element_text(size=text_size),
               legend.text = element_text(size=text_size),
               legend.title = element_text(size=text_size),
               plot.title = element_text(size=plot_title_size),
               panel.border = element_rect(colour="black",fill=NA),
               ...)
}

#' @importFrom ggplot2 ggplot_gtable ggplot_build
extract_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
}

#' @importFrom scales extended_breaks
histogram_breaks <-function(x){
    default_breaks <- extended_breaks()(x)
    if((max(default_breaks)-min(default_breaks)) < 0.01){
        return(median(default_breaks))
    }else{
        return(default_breaks)
    }
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
#' @importFrom ggplot2 ggplot aes geom_point geom_path geom_histogram geom_abline geom_hline facet_wrap scale_color_manual scale_x_continuous scale_y_continuous label_parsed ggtitle xlab ylab
#' @importFrom rlang .data
#' @importFrom stats median
plot_fun <- function(x,type='rating_curve',param=NULL,transformed=F,...){
    args <- list(...)
    cbPalette <- c("green","red","slateblue1","hotpink","#56B4E9","#E69F00","#000000","#999999","#CC79A7","#D55E00","#0072B2","#009E73")
    legal_types <- c('rating_curve','rating_curve_mean','f','beta','sigma_eps','residuals','trace','histogram','r_hat','autocorrelation')
    if(!(type %in% legal_types)){
        stop(cat(paste('Type argument not recognized. Possible types are:\n -',paste(legal_types,collapse='\n - '))))
    }
    mod_params <- get_param_names(class(x),c_param=x$run_info$c_param)
    if(is.null(param)){
        param <- mod_params #defaults to all stage independent model parameters
    }else{
        param <- get_args_rollout(param,mod_params)
    }
    if(type=='trace'){
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
                geom_path(alpha=0.7) +
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
                geom_path(col="#0072B5FF",alpha=0.7) +
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
            facet_wrap(~name_expr,scales='free',labeller=label_parsed) +
            scale_x_continuous(breaks=histogram_breaks) +
            scale_y_continuous(expand=c(0,0,0.05,0)) +
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
            plot_dat$log_Q <- log(plot_dat[, all.vars(x$formula)[1]])

            plot_dat$log_lower <- log(plot_dat$lower)
            plot_dat$log_median <- log(plot_dat$median)
            plot_dat$log_upper <- log(plot_dat$upper)
            p <- ggplot(data=plot_dat) +
                geom_point(aes(x=.data$`log(h-c_hat)`,y=.data$log_Q), size=.9, shape=21, fill="gray60", color="black") +
                geom_path(aes(x=.data$`log(h-c_hat)`,y=.data$log_median)) +
                geom_path(aes(x=.data$`log(h-c_hat)`,y=.data$log_lower),linetype='dashed') +
                geom_path(aes(x=.data$`log(h-c_hat)`,y=.data$log_upper),linetype='dashed') +
                scale_x_continuous(limits=if(!is.null(args$xlim)) args$xlim else c(NA,NA),expand=c(0.01,0)) +
                scale_y_continuous(limits=if(!is.null(args$ylim)) args$ylim else c(NA,NA),expand=c(0.01,0)) +
                xlab(parse(text=x_lab)) +
                ylab(parse(text=y_lab)) +
                theme_bdrc()
        }else{
            #to generate label - latex2exp::TeX('$\\textit{Q}\\lbrack\\textit{m^3/s}\\rbrack$','character')
            x_lab <- "paste('','',italic(paste('Q')),paste('['),italic(paste('m',phantom() ^ {paste('3')},'/s')),paste(']'),'')"
            #to generate label - latex2exp::TeX('$\\textit{h}\\lbrack\\textit{m}\\rbrack$','character')
            y_lab <- "paste('','',italic(paste('h')),paste('['),italic(paste('m')),paste(']'),'')"
            p <- ggplot(data=x[[type]]) +
                geom_point(data=x$data,aes(.data[[all.vars(x$formula)[1]]],.data[[all.vars(x$formula)[2]]]), size=.9, shape=21, fill="gray60", color="black") +
                geom_path(aes(x=.data$median,y=.data$h)) +
                geom_path(aes(x=.data$lower,y=.data$h),linetype='dashed') +
                geom_path(aes(x=.data$upper,y=.data$h),linetype='dashed') +
                scale_x_continuous(limits=if(!is.null(args$xlim)) args$xlim else c(0,max(x$rating_curve$upper,x$data$Q)),expand=c(0.01,0)) +
                scale_y_continuous(limits=if(!is.null(args$ylim)) args$ylim else c(NA,NA),expand=c(0.01,0)) +
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
            scale_x_continuous(limits=if(!is.null(args$xlim)) args$xlim else c(NA,NA),expand=c(0,0)) +
            scale_y_continuous(limits=if(!is.null(args$ylim)) args$ylim else c(0,max(plot_dat$upper)*1.1),expand=c(0,0)) +
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
            scale_x_continuous(if(!is.null(args$xlim)) args$xlim else c(NA,NA),expand=c(0,0)) +
            scale_y_continuous(limits=if(!is.null(args$ylim)) args$ylim else c(NA,NA),expand=c(0,0)) +
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
            scale_x_continuous(limits=if(!is.null(args$xlim)) args$xlim else c(NA,NA),expand=c(0,0)) +
            scale_y_continuous(limits=if(!is.null(args$ylim)) args$ylim else c(min(1,0.9*min(plot_dat$lower)),max(3.5,1.1*max(plot_dat$upper))),expand=c(0,0)) +
            theme_bdrc()
    }else if(type=='residuals'){
        resid_dat <- get_residuals_dat(x)
        resid_lim <- max(abs(resid_dat$r_lower),resid_dat$r_upper,abs(resid_dat$r_median))
        #to generate label - latex2exp::TeX("$log(\\textit{Q})-log(\\textit{\\hat{Q}})$",'character')
        y_lab <- "paste('','log','(','',italic(paste('Q')),')','','-log','(','',italic(paste('',hat(paste('Q')))),')','','')"
        #to generate label - latex2exp::TeX("$log(\\textit{h - \\hat{c}})$",'character')
        x_lab <- "paste('','log','(','',italic(paste('h',phantom() - phantom(),'',hat(paste('c')))),')','','')"
        p <- ggplot(data=resid_dat) +
            geom_point(aes(.data$`log(h-c_hat)`,.data$r_median), size=.9, shape=21, fill="gray60", color="black") +
            geom_path(aes(.data$`log(h-c_hat)`,.data$r_lower),linetype='dashed') +
            geom_path(aes(.data$`log(h-c_hat)`,.data$r_upper),linetype='dashed') +
            geom_abline(intercept=0,slope=0,size=1.1) +
            xlab(parse(text=x_lab)) +
            ylab(parse(text=y_lab)) +
            scale_x_continuous(limits=if(!is.null(args$xlim)) args$xlim else c(NA,NA),expand=c(0.01,0.01)) +
            scale_y_continuous(limits=if(!is.null(args$ylim)) args$ylim else c(-resid_lim,resid_lim)) +
            theme_bdrc()
    }else if(type=='r_hat'){
        rhat_dat <- get_rhat_dat(x,param)
        rhat_dat$Rhat[rhat_dat$Rhat<1] <- 1
        rhat_dat$Rhat[rhat_dat$Rhat>2] <- 2
        param_expr <- parse(text=get_param_expression(param))
        #to generate label - latex2exp::TeX("$\\textit{\\hat{R}}$",'character')
        y_lab <- "paste('','',italic(paste('',hat(paste('R')))),'')"
        p <- ggplot(data=rhat_dat, aes(x=.data$iterations,y=.data$Rhat,color=.data$parameters)) +
             geom_hline(yintercept=1.1,linetype='dashed') +
             geom_line(na.rm=T) +
             scale_x_continuous(limits=c(4*x$run_info$thin+x$run_info$burnin,x$run_info$nr_iter),breaks=c(5000,10000,15000),expand=c(0,0)) +
             scale_y_continuous(limits=c(1,2),breaks=c(1,1.1,1.2,1.4,1.6,1.8,2),expand=c(0,0)) +
             scale_color_manual(values=cbPalette,name=class(x),labels=param_expr) +
             xlab('Iteration') +
             ylab(parse(text=y_lab)) +
             theme_bdrc()
    }else if(type=='autocorrelation'){
        auto_dat <- do.call('rbind',lapply(param,function(p) data.frame(lag=x$autocorrelation$lag,param=p,corr=x$autocorrelation[,p])))
        param_expr <- parse(text=get_param_expression(param))
        max_lag <- nrow(x$autocorrelation)
        p <- ggplot(data=auto_dat, aes(x=.data$lag,y=.data$corr,color=.data$param)) +
             geom_hline(yintercept=0) +
             geom_line() +
             geom_point(size=1) +
             scale_x_continuous(limits=c(1,max_lag),labels=c(1,seq(5,max_lag,5)),breaks=c(1,seq(5,max_lag,5)),expand=c(0,1)) +
             scale_y_continuous(limits=c(min(auto_dat$corr,-1/11),1),expand=c(0,0)) +
             scale_color_manual(values=cbPalette,name=class(x),labels=param_expr) +
             xlab('Lag') +
             ylab('Autocorrelation') +
             theme_bdrc()
    }
    if(!is.null(args$title)){
        p <- p + ggtitle(args$title)
    }
    return(p)
}


#' @importFrom gridExtra arrangeGrob
#' @importFrom grid textGrob gpar unit
#' @importFrom ggplot2 theme guides guide_legend
plot_grob <- function(x,type,transformed=F){
    if(type=='panel'){
        panel_types <- c('rating_curve','residuals','f','sigma_eps')
        plot_list <- lapply(panel_types,function(ty){
            plot_fun(x,type=ty,transformed=transformed)
        })
        p <- do.call(arrangeGrob,c(plot_list,ncol=round(sqrt(length(panel_types)))))
    }else if(type=='convergence_diagnostics'){
        autocorrelation_plot <- plot_fun(x,type='autocorrelation') +
                                theme_bdrc(legend.key.size = unit(0.8, "lines"),
                                           legend.justification = "top",
                                           legend.key.width = unit(0.5,"cm"),
                                           scaling=0.82)
        r_hat_plot <- plot_fun(x,type='r_hat') +
                      theme_bdrc(legend.key.size = unit(0.8, "lines"),
                                 legend.justification = "top",
                                 legend.key.width = unit(1,"cm"),
                                 scaling=0.82)
        legend <- extract_legend(autocorrelation_plot)
        p <- arrangeGrob(arrangeGrob(r_hat_plot+theme(legend.position="none"),
                                     autocorrelation_plot+theme(legend.position="none"),nrow=1),
                         legend,ncol=2,widths=c(4,1))
    }
    return(p)
}

#' @importFrom stats approx
predict_fun <- function(object,newdata=NULL,wide=FALSE){
    if(wide){
        if(!is.null(newdata)){
            stop('newdata must be NULL when wide is TRUE.')
        }else{
            newdata <- seq(ceiling(min(object$rating_curve$h)*100)/100,floor(max(object$rating_curve$h)*100)/100,by=0.01)
        }
    }else{
        newdata <- object$rating_curve$h
    }
    if(class(newdata) !='numeric'){
        stop('newdata must be a vector of type "numeric" or NULL')
    }
    if(any(is.na(newdata))){
        stop('newdata must not include NA')
    }
    if(any(newdata>max(object$rating_curve$h))){
        stop('newdata must contain values within the range of stage values used to fit the rating curve. See "h_max" option to extrapolate the rating curve to higher stages')
    }
    lower_pred <- approx(object$rating_curve$h,object$rating_curve$lower,xout=newdata)$y
    median_pred <- approx(object$rating_curve$h,object$rating_curve$median,xout=newdata)$y
    upper_pred <- approx(object$rating_curve$h,object$rating_curve$upper,xout=newdata)$y
    pred_dat <- data.frame(h=newdata,lower=lower_pred,median=median_pred,upper=upper_pred)
    pred_dat[is.na(pred_dat)] <- 0
    if(wide){
        pred_dat <- predict_wider(pred_dat)
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
#' @param ... further arguments passed to other methods. Currently supports:
#'                     \itemize{
#'                       \item{"title"}{ a character denoting the title of the plot}
#'                       \item{"xlim"}{ numeric vector of length 2, denoting the limits on the x axis of the plot. Only active for types "rating_curve","rating_curve_mean","f","beta","sigma_eps","residuals".}
#'                       \item{"ylim"}{  numeric vector of length 2, denoting the limits on the y axis of the plot. Only active for types "rating_curve","rating_curve_mean","f","beta","sigma_eps","residuals".}
#'                     }
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
autoplot.plm0 <- function(x,type='rating_curve',param=NULL,transformed=F,...){
    plot_fun(x,type=type,param=param,transformed=transformed,...)
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
#' @param ... further arguments passed to other methods. Currently supports:
#'                     \itemize{
#'                       \item{"title"}{ a character denoting the title of the plot}
#'                       \item{"xlim"}{ numeric vector of length 2, denoting the limits on the x axis of the plot. Only active for types "rating_curve","rating_curve_mean","f","beta","sigma_eps","residuals".}
#'                       \item{"ylim"}{  numeric vector of length 2, denoting the limits on the y axis of the plot. Only active for types "rating_curve","rating_curve_mean","f","beta","sigma_eps","residuals".}
#'                     }
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
plot.plm0 <- function(x,type='rating_curve',param=NULL,transformed=F,...){
    grob_types <- c('panel','convergence_diagnostics')
    if(is.null(type) || !(type%in%grob_types)){
        p <- autoplot(x,type=type,param=param,transformed=transformed,...)
        print(p)
    }else{
        p <- plot_grob(x,type=type,transformed=transformed)
        grid.draw(p)
    }
}

#' Predict method for plm0 fit
#'
#' Print the results of a  object
#' @param object an object of class "plm0"
#' @param newdata a numeric vector of stage values for which to predict. If omitted, the stage values in the data are used.
#' @param wide a logical statement determining weather to produce a wide prediction output. If TRUE, then only the predictions median values are presented as a tabular rating curve, with stage changing in decimeter increments with each row and centimeter increments with each column.
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
predict.plm0 <- function(object,newdata=NULL,wide=FALSE,...){
    predict_fun(object,newdata,wide)
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
#' @param ... further arguments passed to other methods. Currently supports:
#'                     \itemize{
#'                       \item{"title"}{ a character denoting the title of the plot}
#'                       \item{"xlim"}{ numeric vector of length 2, denoting the limits on the x axis of the plot. Only active for types "rating_curve","rating_curve_mean","f","beta","sigma_eps","residuals".}
#'                       \item{"ylim"}{  numeric vector of length 2, denoting the limits on the y axis of the plot. Only active for types "rating_curve","rating_curve_mean","f","beta","sigma_eps","residuals".}
#'                     }
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
autoplot.plm <- function(x,type='rating_curve',param=NULL,transformed=F,...){
    plot_fun(x,type=type,param=param,transformed=transformed,...)
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
#' @param ... further arguments passed to other methods. Currently supports:
#'                     \itemize{
#'                       \item{"title"}{ a character denoting the title of the plot}
#'                       \item{"xlim"}{ numeric vector of length 2, denoting the limits on the x axis of the plot. Only active for types "rating_curve","rating_curve_mean","f","beta","sigma_eps","residuals".}
#'                       \item{"ylim"}{  numeric vector of length 2, denoting the limits on the y axis of the plot. Only active for types "rating_curve","rating_curve_mean","f","beta","sigma_eps","residuals".}
#'                     }
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
plot.plm <- function(x,type='rating_curve',param=NULL,transformed=F,...){
    grob_types <- c('panel','convergence_diagnostics')
    if(is.null(type) || !(type%in%grob_types)){
        p <- autoplot(x,type=type,param=param,transformed=transformed,...)
        print(p)
    }else{
        p <- plot_grob(x,type=type,transformed=transformed)
        grid.draw(p)
    }
}

#' Predict method for plm fit
#'
#' Print the results of a  object
#' @param object an object of class "plm"
#' @param newdata a numeric vector of stage values for which to predict. If omitted, the stage values in the data are used.
#' @param wide a logical statement determining weather to produce a wide prediction output. If TRUE, then only the predictions median values are presented as a tabular rating curve, with stage changing in decimeter increments with each row and centimeter increments with each column.
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
predict.plm <- function(object,newdata=NULL,wide=FALSE,...){
    predict_fun(object,newdata,wide)
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
#' @param ... further arguments passed to other methods. Currently supports:
#'                     \itemize{
#'                       \item{"title"}{ a character denoting the title of the plot}
#'                       \item{"xlim"}{ numeric vector of length 2, denoting the limits on the x axis of the plot. Only active for types "rating_curve","rating_curve_mean","f","beta","sigma_eps","residuals".}
#'                       \item{"ylim"}{  numeric vector of length 2, denoting the limits on the y axis of the plot. Only active for types "rating_curve","rating_curve_mean","f","beta","sigma_eps","residuals".}
#'                     }
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
autoplot.gplm0 <- function(x,type='rating_curve',param=NULL,transformed=F,...){
    plot_fun(x,type=type,param=param,transformed=transformed,...)
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
#' @param ... further arguments passed to other methods. Currently supports:
#'                     \itemize{
#'                       \item{"title"}{ a character denoting the title of the plot}
#'                       \item{"xlim"}{ numeric vector of length 2, denoting the limits on the x axis of the plot. Only active for types "rating_curve","rating_curve_mean","f","beta","sigma_eps","residuals".}
#'                       \item{"ylim"}{  numeric vector of length 2, denoting the limits on the y axis of the plot. Only active for types "rating_curve","rating_curve_mean","f","beta","sigma_eps","residuals".}
#'                     }
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
plot.gplm0 <- function(x,type='rating_curve',param=NULL,transformed=F,...){
    grob_types <- c('panel','convergence_diagnostics')
    if(is.null(type) || !(type%in%grob_types)){
        p <- autoplot(x,type=type,param=param,transformed=transformed,...)
        print(p)
    }else{
        p <- plot_grob(x,type=type,transformed=transformed)
        grid.draw(p)
    }
}

#' Predict method for gplm0 fit
#'
#' Print the results of a  object
#' @param object an object of class "gplm0"
#' @param newdata a numeric vector of stage values for which to predict. If omitted, the stage values in the data are used.
#' @param wide a logical statement determining weather to produce a wide prediction output. If TRUE, then only the predictions median values are presented as a tabular rating curve, with stage changing in decimeter increments with each row and centimeter increments with each column.
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
predict.gplm0 <- function(object,newdata=NULL,wide=FALSE,...){
    predict_fun(object,newdata,wide)
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
#' @param ... further arguments passed to other methods. Currently supports:
#'                     \itemize{
#'                       \item{"title"}{ a character denoting the title of the plot}
#'                       \item{"xlim"}{ numeric vector of length 2, denoting the limits on the x axis of the plot. Only active for types "rating_curve","rating_curve_mean","f","beta","sigma_eps","residuals".}
#'                       \item{"ylim"}{  numeric vector of length 2, denoting the limits on the y axis of the plot. Only active for types "rating_curve","rating_curve_mean","f","beta","sigma_eps","residuals".}
#'                     }
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
autoplot.gplm <- function(x,type='rating_curve',param=NULL,transformed=F,...){
    plot_fun(x,type=type,param=param,transformed=transformed,...)
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
#' @param ... further arguments passed to other methods. Currently supports:
#'                     \itemize{
#'                       \item{"title"}{ a character denoting the title of the plot}
#'                       \item{"xlim"}{ numeric vector of length 2, denoting the limits on the x axis of the plot. Only active for types "rating_curve","rating_curve_mean","f","beta","sigma_eps","residuals".}
#'                       \item{"ylim"}{  numeric vector of length 2, denoting the limits on the y axis of the plot. Only active for types "rating_curve","rating_curve_mean","f","beta","sigma_eps","residuals".}
#'                     }
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
plot.gplm <- function(x,type='rating_curve',param=NULL,transformed=F,...){
    grob_types <- c('panel','convergence_diagnostics')
    if(is.null(type) || !(type%in%grob_types)){
        p <- autoplot(x,type=type,param=param,transformed=transformed,...)
        print(p)
    }else{
        p <- plot_grob(x,type=type,transformed=transformed)
        grid.draw(p)
    }
}

#' Predict method for gplm fit
#'
#' Print the results of a  object
#' @param object an object of class "gplm"
#' @param newdata a numeric vector of stage values for which to predict. If omitted, the stage values in the data are used.
#' @param wide a logical statement determining weather to produce a wide prediction output. If TRUE, then only the predictions median values are presented as a tabular rating curve, with stage changing in decimeter increments with each row and centimeter increments with each column.
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
predict.gplm <- function(object,newdata=NULL,wide=FALSE,...){
    predict_fun(object,newdata,wide)
}
