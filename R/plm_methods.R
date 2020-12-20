print_fun <- function(x){
    cat("\nCall:\n",
        paste(deparse(x$formula), sep = "\n", collapse = "\n"), "\n\n", sep = "")
}

summary_fun <- function(x){
    names(x$param_summary)=paste0(names(x$param_summary),c('-2.5%','-50%','-97.5%'))
    cat("\nFormula: \n",
        paste(deparse(x$formula), sep = "\n", collapse = "\n"))
    cat("\nLatent parameters:\n")
    print(x$param_summary[1:2,],row.names = T,digits=3,right=F)
    cat("\nHyperparameters:\n")
    print(x$param_summary[3:nrow(x$param_summary),],row.names = T,digits=3,right=F)
    cat("\nDIC:",x$DIC_summary[,2])
}


plot_fun <- function(x,type=NULL,...,transformed=F){
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
                      strip.text.x = element_text(size = 16))
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
                      strip.text.x = element_text(size = 12),
                      axis.title.x = element_text(size = 12),
                      axis.title.y = element_text(size = 12))
        }
    }else if(type=='histogram'){
        plot_dat <- gather_draws(x,...,transformed=transformed)
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
                  strip.text.x = element_text(size = 16))
    }else if(type=='rating_curve' | type=='rating_curve_mean'){
        if(transformed){
            x_lab <- 'log(h-hat(c))'
            y_lab <- 'log(Q)'
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
            x_lab <- 'Discharge~(m^3/s)'
            p <- ggplot(data=x[[type]]) +
                geom_point(data=x$data,aes_string(all.vars(x$formula)[1],all.vars(x$formula)[2])) +
                geom_line(aes(median,h)) +
                geom_line(aes(lower,h),linetype='dashed') +
                geom_line(aes(upper,h),linetype='dashed') +
                xlab(parse(text=x_lab)) +
                ylab('Stage (m)') +
                theme_classic() +
                theme(axis.text.x = element_text(size = 12),
                      axis.text.y = element_text(size = 12),
                      axis.title.x = element_text(size = 16),
                      axis.title.y = element_text(size = 16))
        }
    }else if(type=='sigma_eps'){
        if(!('sigma_eps_summary' %in% names(x))){
            stop('Plots of type "sigma_eps" are only for models with stage dependent variance, s.a. "bgplm" and "bplm"')
        }
        y_lab <- 'sigma[epsilon](h)'
        p <- ggplot(data=x$sigma_eps_summary) +
            geom_line(aes(h,median)) +
            geom_line(aes(h,lower),linetype='dashed') +
            geom_line(aes(h,upper),linetype='dashed') +
            xlab('Stage (m)') +
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
        y_lab <- 'beta(h)'
        p <- ggplot(data=x$beta_summary) +
            geom_line(aes(h,median)) +
            geom_line(aes(h,lower),linetype='dashed') +
            geom_line(aes(h,upper),linetype='dashed') +
            xlab('Stage (m)') +
            ylab(parse(text=y_lab)) +
            theme_classic() +
            theme(axis.text.x = element_text(size = 12),
                  axis.text.y = element_text(size = 12),
                  axis.title.x = element_text(size = 16),
                  axis.title.y = element_text(size = 16))
    }else if(type=='f'){
        if(!('f_summary' %in% names(x))){
            stop('Plots of type "f" are only for models with stage dependent power law exponent, s.a. "bgplm0" and "bgplm"')
        }
        y_lab <- 'f(h)'
        p <- ggplot(data=x$f_summary) +
            geom_line(aes(median,h)) +
            geom_line(aes(lower,h),linetype='dashed') +
            geom_line(aes(upper,h),linetype='dashed') +
            xlab('Stage (m)') +
            ylab(parse(text=y_lab)) +
            theme_classic() +
            theme(axis.text.x = element_text(size = 12),
                  axis.text.y = element_text(size = 12),
                  axis.title.x = element_text(size = 16),
                  axis.title.y = element_text(size = 16))
    }else if(type=='residuals'){
        resid_dat <- merge(x$rating_curve,x$data,by.x='h',by.y=all.vars(x$formula)[2])
        c_hat <- if(is.null(x$run_info_c_param)) median(x$c_posterior) else x$run_info$c_param
        resid_dat[,'log(h-c_hat)'] <- log(resid_dat$h-c_hat)
        resid_dat$r_median <- log(resid_dat$Q)-log(resid_dat$median)
        resid_dat$r_lower <- log(resid_dat$lower)-log(resid_dat$median)
        resid_dat$r_upper <- log(resid_dat$upper)-log(resid_dat$median)
        y_lab <- 'log(Q)-log(hat(Q))'
        x_lab <- 'log(h-hat(c))'
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


print.bplm0 <- function(x,...){
    print_fun(x)
}

summary.bplm0 <- function(object,...){
    summary_fun(object)
}

plot.bplm0 <- function(x,type='rating_curve',...,transformed=F){
    plot_fun(x,type,...,transformed=transformed)
}

predict.bplm0 <- function(object,newdata=NULL){
    predict_fun(object,newdata)
}

print.bplm <- function(x,...){
    print_fun(x)
}

summary.bplm <- function(object,...){
    summary_fun(object)
}

plot.bplm <- function(x,type='rating_curve',...,transformed=F){
    plot_fun(x,type,...,transformed=transformed)
}

predict.bplm <- function(object,newdata=NULL){
    predict_fun(object,newdata)
}

print.bgplm0 <- function(x,...){
    print_fun(x)
}

summary.bgplm0 <- function(object,...){
    summary_fun(object)
}

plot.bgplm0 <- function(x,type='rating_curve',...,transformed=F){
    plot_fun(x,type,...,transformed=transformed)
}

predict.bgplm0 <- function(object,newdata=NULL){
    predict_fun(object,newdata)
}

print.bgplm <- function(x,...){
    print_fun(x)
}

summary.bgplm <- function(object,...){
    summary_fun(object)
}

plot.bgplm <- function(x,type='rating_curve',...,transformed=F){
    plot_fun(x,type,...,transformed=transformed)
}

predict.bgplm <- function(object,newdata=NULL){
    predict_fun(object,newdata)
}
