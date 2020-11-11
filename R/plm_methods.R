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


plot_fun <- function(x,type=NULL,param=NULL,transformed=F){
    if(type=='trace'){
        plot_dat <- spread_draws(x,param,transformed)
        if('w' %in% names(plot_dat)){
            stop('Plots of type "trace" can only be of stage-independent parameters')
        }
        if('param' %in% names(plot_dat)){
            params <- unique(plot_dat$param)
            param_levels <- get_parameter_levels(params)
            plot_dat$param_expr <- factor(plot_dat$param,levels=param_levels,labels=sapply(param_levels,get_param_expression))
            plot_dat$chain <- factor(as.character(plot_dat$chain),levels=1:max(plot_dat$chain))
            p <- ggplot(plot_dat,aes(iter,value,col=chain)) +
                geom_line() +
                facet_wrap(~param_expr,scales='free',labeller = label_parsed) +
                scale_color_manual(values=c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF"),
                                   name='Chain number') +
                xlab('Iteration') +
                ylab('') +
                theme_classic() +
                theme(strip.background = element_blank(),
                      strip.text.x = element_text(size = 16))
        }else{
            param_name <- names(plot_dat)[3]
            param_expr <- get_param_expression(param_name)
            plot_dat$chain_name <- paste0('Chain nr ',plot_dat$chain)
            p <- ggplot(plot_dat,aes_string('iter',paste0('`',param_name,'`'))) +
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
        plot_dat <- spread_draws(x,param,transformed)
        if('w' %in% names(plot_dat)){
            stop('Plots of type "histogram" can only be of stage-independent parameters')
        }
        if('param' %in% names(plot_dat)){
            params <- unique(plot_dat$param)
            param_levels <- get_parameter_levels(params)
            plot_dat$param_expr <- factor(plot_dat$param,levels=param_levels,labels=sapply(param_levels,get_param_expression))
            plot_dat$chain <- factor(as.character(plot_dat$chain),levels=1:max(plot_dat$chain))
            p <- ggplot(plot_dat,aes(value)) +
                geom_histogram(bins=100,fill="#0072B5FF") +
                facet_wrap(~param_expr,scales='free',labeller = label_parsed) +
                xlab('') +
                ylab('') +
                theme_classic() +
                theme(strip.background = element_blank(),
                      strip.text.x = element_text(size = 16))
        }else{
            param_name <- names(plot_dat)[3]
            param_expr <- get_param_expression(param_name)
            plot_dat$chain_name <- paste0('Chain nr ',plot_dat$chain)
            p <- ggplot(plot_dat,aes_string(paste0('`',param_name,'`')),col="#BC3C29FF") +
                geom_histogram(bins=100,fill="#0072B5FF") +
                xlab(parse(text=param_expr)) +
                theme_classic() +
                theme(axis.title.x = element_text(size = 16),
                      axis.title.y = element_text(size = 12))
        }
    }else if(type=='rating_curve' | type=='rating_curve_mean'){
        if(transformed){
            x_lab <- 'log(W-hat(c))'
            y_lab <- 'log(Q)'
            c_hat <- if(is.null(x$run_info_c_param)) median(x$c_posterior) else x$run_info$c_param
            plot_dat <- merge(x[[type]],x$data,by.x='w',by.y=all.vars(x$formula)[2])
            plot_dat[,'log(w-c_hat)'] <- log(plot_dat$w-c_hat)
            p <- ggplot(data=plot_dat) +
                geom_point(aes(`log(w-c_hat)`,log(Q))) +
                geom_line(aes(`log(w-c_hat)`,log(median))) +
                geom_line(aes(`log(w-c_hat)`,log(lower)),linetype='dashed') +
                geom_line(aes(`log(w-c_hat)`,log(upper)),linetype='dashed') +
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
                geom_line(aes(median,w)) +
                geom_line(aes(lower,w),linetype='dashed') +
                geom_line(aes(upper,w),linetype='dashed') +
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
            geom_line(aes(w,median)) +
            geom_line(aes(w,lower),linetype='dashed') +
            geom_line(aes(w,upper),linetype='dashed') +
            xlab('Stage (m)') +
            ylab(parse(text=y_lab)) +
            theme_classic() +
            theme(axis.text.x = element_text(size = 12),
                  axis.text.y = element_text(size = 12),
                  axis.title.x = element_text(size = 16),
                  axis.title.y = element_text(size = 16))
    }else if(type=='beta' | type=='f'){
        if(!('beta_summary' %in% names(x))){
            stop('Plots of type "beta" are only for models with stage dependent power law exponent, s.a. "bgplm0" and "bgplm"')
        }
        y_lab <- 'beta(h)'
        p <- ggplot(data=x$beta_summary) +
            geom_line(aes(w,median)) +
            geom_line(aes(w,lower),linetype='dashed') +
            geom_line(aes(w,upper),linetype='dashed') +
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
            geom_line(aes(median,w)) +
            geom_line(aes(lower,w),linetype='dashed') +
            geom_line(aes(upper,w),linetype='dashed') +
            xlab('Stage (m)') +
            ylab(parse(text=y_lab)) +
            theme_classic() +
            theme(axis.text.x = element_text(size = 12),
                  axis.text.y = element_text(size = 12),
                  axis.title.x = element_text(size = 16),
                  axis.title.y = element_text(size = 16))
    }else if(type=='residuals'){
        resid_dat <- merge(mod$rating_curve,mod$data,by.x='w',by.y=all.vars(mod$formula)[2])
        c_hat <- if(is.null(mod$run_info_c_param)) median(mod$c_posterior) else mod$run_info$c_param
        resid_dat[,'log(w-c_hat)'] <- log(resid_dat$w-c_hat)
        resid_dat$r_median <- log(resid_dat$Q)-log(resid_dat$median)
        resid_dat$r_lower <- log(resid_dat$lower)-log(resid_dat$median)
        resid_dat$r_upper <- log(resid_dat$upper)-log(resid_dat$median)
        y_lab <- 'log(Q)-log(hat(Q))'
        x_lab <- 'log(W-hat(c))'
        p <- ggplot(data=resid_dat) +
            geom_point(aes(`log(w-c_hat)`,r_median),size=2) +
            geom_line(aes(`log(w-c_hat)`,r_lower),linetype='dashed') +
            geom_line(aes(`log(w-c_hat)`,r_upper),linetype='dashed') +
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
        merged_data <- merge(object$rating_curve,object$data,by.x='w',by.y=all.vars(object$formula)[2])
        pred_dat <- merged_data[,c('w','lower','median','upper')]
    }else{
        if(class(newdata) !='numeric'){
            stop('newdata must be a vector of type "numeric" or NULL')
        }
        if(any(is.na(newdata))){
            stop('newdata must include NA')
        }
        if(any(newdata<min(object$rating_curve$w) | newdata>max(object$rating_curve$w))){
            stop('newdata must contain values within the range of stage values used to fit the rating curve. See "w_max" option to extrapolate the rating curve to higher stages')
        }
        lower_pred <- stats::approx(object$rating_curve$w,object$rating_curve$lower,xout=newdata)$y
        median_pred <- stats::approx(object$rating_curve$w,object$rating_curve$median,xout=newdata)$y
        upper_pred <- stats::approx(object$rating_curve$w,object$rating_curve$upper,xout=newdata)$y
        pred_dat <- data.frame(w=newdata,lower=lower_pred,median=median_pred,upper=upper_pred)
    }
    return(pred_dat)
}


print.bplm0 <- function(x,...){
    print_fun(x)
}

summary.bplm0 <- function(object,...){
    summary_fun(object)
}

plot.bplm0 <- function(x,type='rating_curve',param=NULL,transformed=F,...){
    plot_fun(x,type,param,transformed)
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

plot.bplm <- function(x,type='rating_curve',param=NULL,transformed=F,...){
    plot_fun(x,type,param,transformed)
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

plot.bgplm0 <- function(x,type='rating_curve',param=NULL,transformed=F,...){
    plot_fun(x,type,param,transformed)
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

plot.bgplm <- function(x,type='rating_curve',param=NULL,transformed=F,...){
    plot_fun(x,type,param,transformed)
}

predict.bgplm <- function(object,newdata=NULL){
    predict_fun(object,newdata)
}
