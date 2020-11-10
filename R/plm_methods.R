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
        y_lab <- 'Discharge~(m^3/s)'
        p <- ggplot(data=x[[type]]) +
            geom_point(data=x$data,aes_string(all.vars(x$formula)[1],all.vars(x$formula)[2])) +
            geom_line(aes(median,w)) +
            geom_line(aes(lower,w),linetype='dashed') +
            geom_line(aes(upper,w),linetype='dashed') +
            xlab(parse(text=y_lab)) +
            ylab('Stage (m)') +
            theme_classic() +
            theme(axis.text.x = element_text(size = 12),
                  axis.text.y = element_text(size = 12),
                  axis.title.x = element_text(size = 16),
                  axis.title.y = element_text(size = 16))
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
        resid_dat <- merge(x$rating_curve,x$data,by.x='w',by.y=all.vars(x$formula)[2])
        resid_dat$r_median <- log(resid_dat$Q)-log(resid_dat$median)
        resid_dat$r_lower <- log(resid_dat$lower)-log(resid_dat$median)
        resid_dat$r_upper <- log(resid_dat$upper)-log(resid_dat$median)
        y_lab <- 'log(Q)-log(hat(Q))'
        p <- ggplot(data=resid_dat) +
            geom_point(aes(w,r_median),size=2) +
            geom_line(aes(w,r_lower),linetype='dashed') +
            geom_line(aes(w,r_upper),linetype='dashed') +
            geom_abline(intercept=0,slope=0,size=1.1) +
            xlab('Stage (m)') +
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
        merged_data <- merge(object$rating_curve,object$data,by.x='w',by.y=all.vars(x$formula)[2])
        pred_dat <- merged_data[,c('w','lower','median','upper')]
    }else{

    }
    return(pred_dat)
}


print.bplm0 <- function(x,...){
    print_fun(x)
}

summary.bplm0 <- function(x,...){
    summary_fun(x)
}

plot.bplm0 <- function(x,type='rating_curve',...){
}

predict.bplm0 <- function(x,...){
}

print.bplm <- function(x,...){
    print_fun(x)
}

summary.bplm <- function(x,...){
    summary_fun(x)
}

plot.bplm <- function(x,type='rating_curve',...){
}

predict.bplm <- function(x,...){
}

print.bgplm0 <- function(x,...){
    print_fun(x)
}

summary.bgplm0 <- function(x,...){
    summary_fun(x)
}

plot.bgplm0 <- function(x,type='rating_curve',...){
}

predict.bgplm0 <- function(x,...){
}

print.bgplm <- function(x,...){
    print_fun(x)
}

summary.bgplm <- function(x,...){
    summary_fun(x)
}

plot.bgplm <- function(x,type='rating_curve',param=NULL,transformed=F,...){
    plot_fun(x,type,param,transformed)
}

predict.bgplm <- function(x,...){
}
