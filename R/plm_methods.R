print_fun <- function(x){
    cat(paste0(class(x)," - Call:\n"),
        paste(deparse(x$formula), collapse = "\n"), "\n")
}

summary_fun <- function(x){
    param_summary <- x$param_summary[,c('lower','median','upper')]
    names(param_summary) <- paste0(names(param_summary),c('-2.5%','-50%','-97.5%'))
    cat("\nFormula: \n",
        paste(deparse(x$formula), sep = "\n", collapse = "\n"))
    cat("\nLatent parameters:\n")
    print(param_summary[1:2,],row.names = TRUE,digits=3,right=FALSE)
    cat("\nHyperparameters:\n")
    print(param_summary[3:nrow(param_summary),],row.names = TRUE,digits=3,right=FALSE)
    cat("\nWAIC:",x$WAIC)
}

#' Custom bdrc theme
#'
#' @param ... not used in this function
#' @param scaling a numerical value which can be used to scale up or down the size of the text and titles of a plot that uses \code{theme_bdrc}. Defaults to 1.
#' @return returns a theme object for the package
#' @importFrom ggplot2 %+replace% theme_classic theme element_text element_blank element_rect
#' @keywords internal
theme_bdrc <- function(...,scaling=1){
    title_size <- scaling*12
    text_size <- scaling*12
    plot_title_size=scaling*12
    theme_classic() %+replace%
        theme( #text = element_text(family="Times", face="plain"),
               strip.background = element_blank(),
               strip.text.x = element_text(size = title_size),
               axis.title.x = element_text(size = title_size),
               axis.title.y = element_text(size = title_size,angle=90),
               axis.text.x = element_text(size = text_size),
               axis.text.y = element_text(size = text_size),
               legend.text = element_text(size = text_size),
               legend.title = element_text(size = text_size),
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
#' Visualize results from model objects in bdrc, plm0, plm, gplm0,gplm
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
#'                       \item{"panel"}{ to plot a 2x2 panel of plots: "rating curve", "residuals", "f" and "sigma_eps"}
#'                    }
#' @param param a character vector with the parameters to plot. Defaults to NULL and is only used if type is "trace" or "histogram". Allowed values are the parameters given in the model summary of x as well as "hyperparameters" or "latent_parameters" for specific groups of parameters.
#' @param transformed a logical value indicating whether the quantity should be plotted on a transformed scale used during the Bayesian inference. Defaults to FALSE.
#' @param title a character denoting the title of the plot. Defaults to NULL, i.e. no title.
#' @param xlim numeric vector of length 2, denoting the limits on the x axis of the plot. Applicable for types "rating_curve","rating_curve_mean","f","beta","sigma_eps","residuals".
#' @param ylim numeric vector of length 2, denoting the limits on the y axis of the plot. Applicable for types "rating_curve","rating_curve_mean","f","beta","sigma_eps","residuals".
#' @return returns an object of class ggplot2.
#' @importFrom ggplot2 ggplot aes geom_point geom_path geom_histogram geom_abline geom_hline geom_smooth facet_wrap scale_color_manual scale_x_continuous scale_y_continuous expansion label_parsed ggtitle xlab ylab geom_blank margin element_text theme
#' @importFrom rlang .data
#' @importFrom stats median
#' @keywords internal
plot_fun <- function(x,type='rating_curve',param=NULL,transformed=FALSE,title=NULL,xlim=NULL,ylim=NULL,...){
    color_palette <- c("green","red","slateblue1","hotpink","#56B4E9","#E69F00","#000000","#999999","#CC79A7","#D55E00","#0072B2","#009E73")
    legal_types <- c('rating_curve','rating_curve_mean','f','beta','sigma_eps','residuals','trace','histogram','r_hat','autocorrelation')
    if(!(type %in% legal_types)){
        stop(paste('Type argument not recognized. Possible types are:\n -',paste(legal_types,collapse='\n - ')))
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
                geom_path(alpha=0.7,linewidth=0.3) +
                facet_wrap(~name_expr,scales='free',labeller = label_parsed) +
                scale_color_manual(values=c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF"),
                                   name='Chain number') +
                xlab('Iteration') +
                ylab('') +
                ggtitle(if(!is.null(title)) title else "Traceplots") +
                theme_bdrc()
        }else{
            param_expr <- get_param_expression(params)
            plot_dat$chain_name <- paste0('Chain nr ',plot_dat$chain)
            p <- ggplot(plot_dat,aes(x=.data$iter,y=.data$value)) +
                geom_path(col="#0072B5FF",linewidth=0.3) +
                facet_wrap(~chain_name,scales='free') +
                xlab('Iteration') +
                ylab(parse(text=param_expr)) +
                ggtitle(if(!is.null(title)) title else "Traceplots") +
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
            facet_wrap(~name_expr,scales='free',labeller=label_parsed,strip.position = "bottom") +
            scale_x_continuous(breaks=histogram_breaks) +
            scale_y_continuous(expand=c(0,0,0.05,0)) +
            xlab('') +
            ylab('') +
            ggtitle(if(!is.null(title)) title else if( length(param)==1 ) "Histogram of posterior draws" else "Histograms of posterior draws" ) +
            theme_bdrc() +
            theme(plot.title = element_text( vjust = 2 ),
                  strip.placement = "outside",
                  plot.margin =  ggplot2::margin(t = 10, r = 10, b = -10, l = 0, unit = "pt"))
    }else if(type=='rating_curve' | type=='rating_curve_mean'){
        if(transformed){
            x_lab <- "paste('','',log,,,,'(','',italic(paste('h-',hat(paste('c')))),')','','')"
            y_lab <- "paste('','',log,,,,'(','',italic(paste('Q')),')','','')"
            c_hat <- if(is.null(x$run_info$c_param)) median(x$c_posterior) else x$run_info$c_param
            h_min <- min(x$data[[all.vars(x$formula)[2]]])
            plot_dat <- merge(x[[type]][x[[type]]$h>=h_min,],x$data,by.x='h',by.y=all.vars(x$formula)[2],all.x=TRUE)
            plot_dat[,'log(h-c_hat)'] <- log(plot_dat$h-c_hat)
            plot_dat$log_Q <- log(plot_dat[, all.vars(x$formula)[1]])

            plot_dat$log_lower <- log(plot_dat$lower)
            plot_dat$log_median <- log(plot_dat$median)
            plot_dat$log_upper <- log(plot_dat$upper)
            p <- ggplot(data=plot_dat) +
                geom_path(aes(x=.data$`log(h-c_hat)`,y=.data$log_median),alpha=0.95) +
                geom_path(aes(x=.data$`log(h-c_hat)`,y=.data$log_lower),linetype='dashed',alpha=0.95) +
                geom_path(aes(x=.data$`log(h-c_hat)`,y=.data$log_upper),linetype='dashed',alpha=0.95) +
                geom_point(data=plot_dat[!is.na(plot_dat$log_Q),],aes(x=.data$`log(h-c_hat)`,y=.data$log_Q), size=.9, shape=21, fill="gray60", color="black",alpha=0.95) +
                scale_x_continuous(limits=if(!is.null(xlim)) xlim else c(NA,NA),expand=c(0.01,0)) +
                scale_y_continuous(limits=if(!is.null(ylim)) ylim else c(NA,NA),expand=c(0.01,0)) +
                xlab(parse(text=x_lab)) +
                ylab(parse(text=y_lab)) +
                ggtitle(if(!is.null(title)) title else "Log-transformed rating curve") +
                theme_bdrc() +
                theme(plot.title = element_text( vjust = 2 ))
        }else{
            x_lab <- "paste('','',italic(paste('Q')),paste('['),italic(paste('m',phantom() ^ {paste('3')},'/s')),paste(']'),'')"
            y_lab <- "paste('','',italic(paste('h')),paste('['),italic(paste('m')),paste(']'),'')"
            p <- ggplot(data=x[[type]]) +
                geom_path(aes(x=.data$median,y=.data$h),alpha=0.95) +
                geom_path(aes(x=.data$lower,y=.data$h),linetype='dashed',alpha=0.95) +
                geom_path(aes(x=.data$upper,y=.data$h),linetype='dashed',alpha=0.95) +
                geom_point(data=x$data,aes(.data[[all.vars(x$formula)[1]]],.data[[all.vars(x$formula)[2]]]), size=.9, shape=21, fill="gray60", color="black",alpha=0.95) +
                scale_x_continuous(limits=if(!is.null(xlim)) xlim else c(NA,NA),expand=expansion(mult=0.01)) +
                #scale_x_continuous(limits=if(!is.null(xlim)) xlim else c(0,max(x$rating_curve$upper,x$data$Q)),expand=c(0.01,0)) +
                scale_y_continuous(limits=if(!is.null(ylim)) ylim else c(NA,NA),expand=c(0.01,0)) +
                xlab(parse(text=x_lab)) +
                ylab(parse(text=y_lab)) +
                ggtitle(if(!is.null(title)) title else "Rating curve") +
                theme_bdrc() +
                theme(plot.title = element_text( vjust = 2 ))
        }
    }else if(type=='sigma_eps'){
        x_lab <- "paste('','',italic(paste('h')),paste('['),italic(paste('m')),paste(']'),'')"
        h_in_data <- x$data[,all.vars(x$formula)[2],drop=TRUE]
        if('sigma_eps_summary' %in% names(x)){
            y_lab <- "paste('','',sigma,,,,phantom() [ {paste('',epsilon,,,)} ],'(','',italic(paste('h')),')','','')"
            plot_dat <- x$sigma_eps_summary[x$sigma_eps_summary$h>=min(h_in_data) & x$sigma_eps_summary$h<=max(h_in_data),]
        }else{
            y_lab <- "paste('','',sigma,,,,phantom() [ {paste('',epsilon,,,)} ],'')"
            plot_dat <- data.frame(h=x$data[,all.vars(x$formula)[2],drop=TRUE],
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
            scale_x_continuous(limits=if(!is.null(xlim)) xlim else c(NA,NA),expand=c(0,0)) +
            scale_y_continuous(limits=if(!is.null(ylim)) ylim else c(0,max(plot_dat$upper)*1.1),expand=c(0,0)) +
            ggtitle(if(!is.null(title)) title else "Std. dev. of the error terms") +
            theme_bdrc() +
            theme(plot.title = element_text( vjust = 2 ))
    }else if(type=='beta'){
        if(!('beta_summary' %in% names(x))){
            stop('Plots of type "beta" are only for models with stage dependent power law exponent, s.a. "gplm0" and "gplm"')
        }
        x_lab <- "paste('','',italic(paste('h')),paste('['),italic(paste('m')),paste(']'),'')"
        y_lab <- "paste('','',beta,,,,'(','',italic(paste('h')),')','','')"
        h_in_data <- x$data[,all.vars(x$formula)[2],drop=TRUE]
        p <- ggplot(data=x$beta_summary[x$beta_summary$h>=min(h_in_data) & x$beta_summary$h<=max(h_in_data),]) +
            geom_path(aes(.data$h,.data$median)) +
            geom_path(aes(.data$h,.data$lower),linetype='dashed') +
            geom_path(aes(.data$h,.data$upper),linetype='dashed') +
            xlab(parse(text=x_lab)) +
            ylab(parse(text=y_lab)) +
            scale_x_continuous(if(!is.null(xlim)) xlim else c(NA,NA),expand=c(0,0)) +
            scale_y_continuous(limits=if(!is.null(ylim)) ylim else c(NA,NA),expand=expansion(mult=rep(.05,2))) +
            ggtitle(if(!is.null(title)) title else "Power-law exponent deviations") +
            theme_bdrc() +
            theme(plot.title = element_text( vjust = 2 ))
    }else if(type=='f'){
        x_lab <- "paste('','',italic(paste('h')),paste('['),italic(paste('m')),paste(']'),'')"
        h_in_data <- x$data[,all.vars(x$formula)[2],drop=TRUE]
        if('f_summary' %in% names(x)){
            y_lab <- "paste('','',italic(paste('b+',beta,,,,'(','h',')','')),'')"
            plot_dat <- x$f_summary[x$f_summary$h>=min(h_in_data) & x$f_summary$h<=max(h_in_data),]
        }else{
            y_lab <- "paste('','',italic(paste('b')),'')"
            plot_dat <- data.frame(h=x$data[,all.vars(x$formula)[2],drop=TRUE],
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
            scale_x_continuous(limits=if(!is.null(xlim)) xlim else c(NA,NA),expand=c(0,0)) +
            scale_y_continuous(limits=if(!is.null(ylim)) ylim else c(min(1,0.9*min(plot_dat$lower)),max(3.5,1.1*max(plot_dat$upper))),expand=c(0,0)) +
            ggtitle(if(!is.null(title)) title else "Power-law exponent") +
            theme_bdrc() +
            theme(plot.title = element_text( vjust = 2 ))
    }else if(type=='residuals'){
        resid_dat <- get_residuals_dat(x)
        y_lab <- "paste('','log','(','',italic(paste('Q')),')','','-log','(','',italic(paste('',hat(paste('Q')))),')','','')"
        x_lab <- "paste('','log','(','',italic(paste('h',phantom() - phantom(),'',hat(paste('c')))),')','','')"
        method <- 'loess'
        span <- 0.3
        p <- ggplot(data=resid_dat) +
            geom_hline(yintercept=0,linewidth=0.8,alpha=.95) +
            geom_point(data=resid_dat[!is.na(resid_dat$Q),],aes(.data$`log(h-c_hat)`,.data$r_median), size=.9, shape=21, fill="gray60", color="black",alpha=0.95) +
            geom_blank(aes(y=-.data$r_median)) +
            geom_blank(aes(y=-.data$r_upper)) +
            geom_blank(aes(y=-.data$r_lower)) +
            geom_blank(aes(y=-.data$m_upper)) +
            geom_blank(aes(y=-.data$m_lower)) +
            geom_smooth(aes(x=.data$`log(h-c_hat)`,y=.data$r_upper),span=span,se=FALSE,stat = "smooth",color='black',linetype='dashed',linewidth=0.5,alpha=0.95,method=method,formula='y~x') +
            geom_smooth(aes(x=.data$`log(h-c_hat)`,y=.data$r_lower),span=span,se=FALSE,stat = "smooth",color='black',linetype='dashed',linewidth=0.5,alpha=0.95,method=method,formula='y~x') +
            geom_smooth(aes(x=.data$`log(h-c_hat)`,y=.data$m_upper),span=span,se=FALSE,stat = "smooth",color='black',linetype='solid',linewidth=0.3,alpha=0.95,method=method,formula='y~x') +
            geom_smooth(aes(x=.data$`log(h-c_hat)`,y=.data$m_lower),span=span,se=FALSE,stat = "smooth",color='black',linetype='solid',linewidth=0.3,alpha=0.95,method=method,formula='y~x') +
            xlab(parse(text=x_lab)) +
            ylab(parse(text=y_lab)) +
            scale_x_continuous(limits=if(!is.null(xlim)) xlim else c(NA,NA),expand=expansion(mult=rep(.01,2))) +
            scale_y_continuous(limits=if(!is.null(ylim)) ylim else c(NA,NA),expand=expansion(mult=rep(.05,2))) +
            ggtitle(if(!is.null(title)) title else "Residuals") +
            theme_bdrc() +
            theme(plot.title = element_text( vjust = 2 ))
    }else if(type=='r_hat'){
        rhat_dat <- get_rhat_dat(x,param)
        rhat_dat$Rhat[rhat_dat$Rhat<1] <- 1
        rhat_dat$Rhat[rhat_dat$Rhat>2] <- 2
        param_expr <- parse(text=get_param_expression(param))
        y_lab <- "paste('','',italic(paste('',hat(paste('R')))),'')"
        p <- ggplot(data=rhat_dat, aes(x=.data$iterations,y=.data$Rhat,color=.data$parameters)) +
             geom_hline(yintercept=1.1,linetype='dashed') +
             geom_line(na.rm=TRUE) +
             scale_x_continuous(limits=c(4*x$run_info$thin+x$run_info$burnin,x$run_info$nr_iter),breaks=c(5000,10000,15000),expand=c(0,0)) +
             scale_y_continuous(limits=c(1,2),breaks=c(1,1.1,1.2,1.4,1.6,1.8,2),expand=c(0,0)) +
             scale_color_manual(values=color_palette,name=class(x),labels=param_expr) +
             xlab('Iteration') +
             ylab(parse(text=y_lab)) +
             ggtitle(if(!is.null(title)) title else "Gelman-Rubin statistic") +
             theme_bdrc() +
             theme(plot.title = element_text( vjust = 2 ),
                   axis.title.y = element_text( vjust = 3 ),
                   plot.margin = ggplot2::margin( t=7, r=7, b=7, l=12, unit = "pt" ) )
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
             scale_color_manual(values=color_palette,name=class(x),labels=param_expr) +
             xlab('Lag') +
             ylab('Sample autocorrelation') +
             ggtitle(if(!is.null(title)) title else "Autocorrelation in posterior draws") +
             theme_bdrc() +
             theme(plot.title = element_text( vjust = 2 ),
                   axis.title.y = element_text( vjust = 3 ),
                   plot.margin = ggplot2::margin( t=7, r=7, b=7, l=12, unit = "pt" ) )
    }
    return(p)
}


#' @importFrom gridExtra arrangeGrob
#' @importFrom grid textGrob gpar unit unit.pmax
#' @importFrom ggplot2 theme guides guide_legend ggplotGrob
plot_grob <- function(x,type,transformed=FALSE){
    if(type=='panel'){
        panel_types <- c('rating_curve','residuals','f','sigma_eps')
        grob_list <- lapply(panel_types,function(ty){
            ggplotGrob(plot_fun(x,type=ty,transformed=transformed))
        })
        maxHeight <-  unit.pmax( grob_list[[1]]$heights[2:9], grob_list[[2]]$heights[2:9],
                                 grob_list[[3]]$heights[2:9], grob_list[[4]]$heights[2:9])
        maxWidth <-  unit.pmax( grob_list[[1]]$widths[2:5], grob_list[[2]]$widths[2:5],
                                grob_list[[3]]$widths[2:5], grob_list[[4]]$widths[2:5])
        for(j in 1:4){
            grob_list[[j]]$heights[2:9] <- as.list(maxHeight)
            grob_list[[j]]$widths[2:5] <- as.list(maxWidth)
        }
        p <- do.call(arrangeGrob,c(grob_list,ncol=round(sqrt(length(panel_types)))))
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
        if(is.null(newdata)){
            newdata <- object$rating_curve$h
        }
    }
    if(!inherits(newdata,'numeric')){
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

#' Print method for discharge rating curves
#'
#' Print a discharge rating curve model object
#' @param x an object of class "plm0", "plm", "gplm0" or "gplm".
#' @param ... not used in this function
#' @seealso \code{\link{plm0}}, \code{\link{plm}}, \code{\link{gplm0}}, \code{\link{gplm}} for fitting a discharge rating curve and \code{\link{summary.plm0}}, \code{\link{summary.plm}}, \code{\link{summary.gplm0}} and \code{\link{summary.gplm}} for summaries. It is also useful to look at \code{\link{plot.plm0}}, \code{\link{plot.plm}}, \code{\link{plot.gplm0}} and \code{\link{plot.gplm}} to help visualize all aspects of the fitted discharge rating curve. Additionally, \code{\link{spread_draws}} and \code{\link{spread_draws}} help working directly with the MCMC samples.
#' @describeIn print.plm0 Print method for plm0
#' @export
print.plm0 <- function(x,...){
    print_fun(x)
}

#' Summary method for discharge rating curves
#'
#' Summarize a discharge rating curve model object
#' @param object an object of class "plm0", "plm", "gplm0" or "gplm".
#' @param ... Not used for this function
#' @seealso \code{\link{plm0}}, \code{\link{plm}}, \code{\link{gplm0}} and \code{\link{gplm}} for fitting a discharge rating curve. It is also useful to look at \code{\link{plot.plm0}}, \code{\link{plot.plm}}, \code{\link{plot.gplm0}} and \code{\link{plot.gplm}} to help visualize all aspects of the fitted discharge rating curve. Additionally, \code{\link{spread_draws}} and \code{\link{spread_draws}} help working directly with the MCMC samples.
#' @examples
#' \donttest{
#' data(krokfors)
#' set.seed(1)
#' plm0.fit <- plm0(formula=Q~W,data=krokfors,num_cores=2)
#' summary(plm0.fit)
#' }
#' @describeIn summary.plm0 Summary method for plm0
#' @export
summary.plm0 <- function(object,...){
    summary_fun(object)
}

#' Autoplot method for discharge rating curves
#'
#' Visualize discharge rating curve model objects
#' @param object an object of class "plm0","plm","gplm0" or "gplm".
#' @param ... other plotting parameters (not used in this function)
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
#' @param title a character denoting the title of the plot
#' @param xlim numeric vector of length 2, denoting the limits on the x axis of the plot. Applicable for types "rating_curve","rating_curve_mean","f","beta","sigma_eps","residuals".
#' @param ylim numeric vector of length 2, denoting the limits on the y axis of the plot. Applicable for types "rating_curve","rating_curve_mean","f","beta","sigma_eps","residuals".
#'
#' @return returns an object of class "ggplot2".
#' @seealso \code{\link{plm0}}, \code{\link{plm}}, \code{\link{gplm0}} and \code{\link{gplm}} for fitting a discharge rating curve and \code{\link{summary.plm0}}, \code{\link{summary.plm}}, \code{\link{summary.gplm0}} and \code{\link{summary.gplm}} for summaries. It is also useful to look at \code{\link{spread_draws}} and \code{\link{gather_draws}} to work directly with the MCMC samples.
#' @examples
#' \donttest{
#' library(ggplot2)
#' data(krokfors)
#' set.seed(1)
#' plm0.fit <- plm0(Q~W,krokfors,num_cores=2)
#' autoplot(plm0.fit)
#' autoplot(plm0.fit,transformed=TRUE)
#' autoplot(plm0.fit,type='histogram',param='c')
#' autoplot(plm0.fit,type='histogram',param='c',transformed=TRUE)
#' autoplot(plm0.fit,type='histogram',param='hyperparameters')
#' autoplot(plm0.fit,type='histogram',param='latent_parameters')
#' autoplot(plm0.fit,type='residuals')
#' autoplot(plm0.fit,type='f')
#' autoplot(plm0.fit,type='sigma_eps')
#' }
#' @describeIn autoplot.plm0 Autoplot method for plm0
#' @importFrom ggplot2 autoplot
#' @export
autoplot.plm0 <- function(object,...,type='rating_curve',param=NULL,transformed=FALSE,title=NULL,xlim=NULL,ylim=NULL){
    plot_fun(object,type=type,param=param,transformed=transformed,title=title,xlim=xlim,ylim=ylim)
}

#' Plot method for discharge rating curves
#'
#' Visualize discharge rating curve model objects
#' @param x object of class "plm0", "plm", "gplm0" or "gplm".
#' @param ... other plotting parameters (not used in this function)
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
#'                       \item{"panel"}{ to plot a 2x2 panel of plots: "rating curve", "residuals", "f" and "sigma_eps"}
#'                    }
#' @param param a character vector with the parameters to plot. Defaults to NULL and is only used if type is "trace" or "histogram". Allowed values are the parameters given in the model summary of x as well as "hyperparameters" or "latent_parameters" for specific groups of parameters.
#' @param transformed a logical value indicating whether the quantity should be plotted on a transformed scale used during the Bayesian inference. Defaults to FALSE.
#' @param title a character denoting the title of the plot
#' @param xlim numeric vector of length 2, denoting the limits on the x axis of the plot. Applicable for types "rating_curve","rating_curve_mean","f","beta","sigma_eps","residuals".
#' @param ylim numeric vector of length 2, denoting the limits on the y axis of the plot. Applicable for types "rating_curve","rating_curve_mean","f","beta","sigma_eps","residuals".
#'
#' @return No return value, called for side effects.
#' @seealso \code{\link{plm0}}, \code{\link{plm}}, \code{\link{gplm0}} and \code{\link{gplm}} for fitting a discharge rating curve and \code{\link{summary.plm0}}, \code{\link{summary.plm}}, \code{\link{summary.gplm0}} and \code{\link{summary.gplm}} for summaries. It is also useful to look at \code{\link{spread_draws}} and \code{\link{gather_draws}} to work directly with the MCMC samples.
#' @examples
#' \donttest{
#' data(krokfors)
#' set.seed(1)
#' plm0.fit <- plm0(formula=Q~W,data=krokfors,num_cores=2)
#'
#' plot(plm0.fit)
#' plot(plm0.fit,transformed=TRUE)
#' plot(plm0.fit,type='histogram',param='c')
#' plot(plm0.fit,type='histogram',param='c',transformed=TRUE)
#' plot(plm0.fit,type='histogram',param='hyperparameters')
#' plot(plm0.fit,type='histogram',param='latent_parameters')
#' plot(plm0.fit,type='residuals')
#' plot(plm0.fit,type='f')
#' plot(plm0.fit,type='sigma_eps')
#' }
#' @describeIn plot.plm0 Plot method for plm0
#' @export
#' @importFrom grid grid.draw
#' @importFrom ggplot2 autoplot
plot.plm0 <- function(x,...,type='rating_curve',param=NULL,transformed=FALSE,title=NULL,xlim=NULL,ylim=NULL){
    grob_types <- c('panel','convergence_diagnostics')
    if(is.null(type) || !(type%in%grob_types)){
        p <- autoplot(x,type=type,param=param,transformed=transformed,title=title,xlim=xlim,ylim=ylim)
        print(p)
    }else{
        p <- plot_grob(x,type=type,transformed=transformed)
        grid.draw(p)
    }
}

#' Predict method for discharge rating curves
#'
#' Predict the discharge for given stage values based on a discharge rating curve model object.
#' @param object an object of class "plm0", "plm", "gplm0" or "gplm".
#' @param newdata a numeric vector of stage values for which to predict. If omitted, the stage values in the data are used.
#' @param wide a logical value denoting whether to produce a wide prediction output.If TRUE, then the output is a table with median prediction values for an equally spaced grid of stages with 1 cm increments, each row containing predictions in a decimeter range of stages.
#' @param ... not used in this function
#' @return an object of class "data.frame" with four columns, h (stage), lower (2.5\% posterior predictive quantile), median (50\% posterior predictive quantile), upper (97.5\% posterior predictive quantile). If wide=TRUE, a matrix as described above (see wide parameter) is returned.
#' @seealso \code{\link{plm0}}, \code{\link{plm}}, \code{\link{gplm0}} and \code{\link{gplm}} for fitting a discharge rating curve and \code{\link{summary.plm0}}, \code{\link{summary.plm}}, \code{\link{summary.gplm0}} and \code{\link{summary.gplm}} for summaries. It is also useful to look at \code{\link{plot.plm0}}, \code{\link{plot.plm}}, \code{\link{plot.gplm0}} and \code{\link{plot.gplm}} to help visualize all aspects of the fitted discharge rating curve. Additionally, \code{\link{spread_draws}} and \code{\link{spread_draws}} help working directly with the MCMC samples.
#' @examples
#' \donttest{
#' data(krokfors)
#' set.seed(1)
#' plm0.fit <- plm0(formula=Q~W,data=krokfors,h_max=10,num_cores=2)
#' #predict rating curve on a equally 10 cm spaced grid from 9 to 10 meters
#' predict(plm0.fit,newdata=seq(9,10,by=0.1))
#' }
#' @describeIn predict.plm0 Predict method for plm0
#' @export
predict.plm0 <- function(object,...,newdata=NULL,wide=FALSE){
    predict_fun(object,newdata=newdata,wide=wide)
}

#' @describeIn print.plm0 Print method for plm
#' @export
#'
print.plm <- function(x,...){
    print_fun(x)
}

#' @describeIn summary.plm0 Summary method for plm
#' @export
summary.plm <- function(object,...){
    summary_fun(object)
}

#' @describeIn autoplot.plm0 Autoplot method for plm
#' @export
autoplot.plm <- function(object,...,type='rating_curve',param=NULL,transformed=FALSE,title=NULL,xlim=NULL,ylim=NULL){
    plot_fun(object,type=type,param=param,transformed=transformed,title=title,xlim=xlim,ylim=ylim)
}

#' @describeIn plot.plm0 Plot method for plm
#' @export
#' @importFrom grid grid.draw
#' @importFrom ggplot2 autoplot
plot.plm <- function(x,...,type='rating_curve',param=NULL,transformed=FALSE,title=NULL,xlim=NULL,ylim=NULL){
    grob_types <- c('panel','convergence_diagnostics')
    if(is.null(type) || !(type%in%grob_types)){
        p <- autoplot(x,type=type,param=param,transformed=transformed,title=title,xlim=xlim,ylim=ylim)
        print(p)
    }else{
        p <- plot_grob(x,type=type,transformed=transformed)
        grid.draw(p)
    }
}

#' @describeIn predict.plm0 Predict method for plm
#' @export
predict.plm <- function(object,...,newdata=NULL,wide=FALSE){
    predict_fun(object,newdata=newdata,wide=wide)
}

#' @describeIn print.plm0 Print method for gplm0
#' @export
print.gplm0 <- function(x,...){
    print_fun(x)
}

#' @describeIn summary.plm0 Summary method for gplm0
#' @export
summary.gplm0 <- function(object,...){
    summary_fun(object)
}

#' @describeIn autoplot.plm0 Autoplot method for gplm0
#' @export
autoplot.gplm0 <- function(object,...,type='rating_curve',param=NULL,transformed=FALSE,title=NULL,xlim=NULL,ylim=NULL){
    plot_fun(object,type=type,param=param,transformed=transformed,title=title,xlim=xlim,ylim=ylim)
}

#' @describeIn plot.plm0 Plot method for gplm0
#' @export
#' @importFrom grid grid.draw
#' @importFrom ggplot2 autoplot
plot.gplm0 <- function(x,...,type='rating_curve',param=NULL,transformed=FALSE,title=NULL,xlim=NULL,ylim=NULL){
    grob_types <- c('panel','convergence_diagnostics')
    if(is.null(type) || !(type%in%grob_types)){
        p <- autoplot(x,type=type,param=param,transformed=transformed,title=title,xlim=xlim,ylim=ylim)
        print(p)
    }else{
        p <- plot_grob(x,type=type,transformed=transformed)
        grid.draw(p)
    }
}

#' @describeIn predict.plm0 Predict method for gplm0
#' @export
predict.gplm0 <- function(object,...,newdata=NULL,wide=FALSE){
    predict_fun(object,newdata=newdata,wide=wide)
}

#' @describeIn print.plm0 Print method for gplm
#' @export
print.gplm <- function(x,...){
    print_fun(x)
}

#' @describeIn summary.plm0 Summary method for gplm
#' @export
summary.gplm <- function(object,...){
    summary_fun(object)
}

#' @describeIn autoplot.plm0 Autoplot method for gplm
#' @export
autoplot.gplm <- function(object,...,type='rating_curve',param=NULL,transformed=FALSE,title=NULL,xlim=NULL,ylim=NULL){
    plot_fun(object,type=type,param=param,transformed=transformed,title=title,xlim=xlim,ylim=ylim)
}

#' @describeIn plot.plm0 Plot method for gplm
#' @export
#' @importFrom grid grid.draw
#' @importFrom ggplot2 autoplot
plot.gplm <- function(x,...,type='rating_curve',param=NULL,transformed=FALSE,title=NULL,xlim=NULL,ylim=NULL){
    grob_types <- c('panel','convergence_diagnostics')
    if(is.null(type) || !(type%in%grob_types)){
        p <- autoplot(x,type=type,param=param,transformed=transformed,title=title,xlim=xlim,ylim=ylim)
        print(p)
    }else{
        p <- plot_grob(x,type=type,transformed=transformed)
        grid.draw(p)
    }
}

#' @describeIn predict.plm0 Predict method for gplm
#' @export
predict.gplm <- function(object,...,newdata=NULL,wide=FALSE){
    predict_fun(object,newdata=newdata,wide=wide)
}
