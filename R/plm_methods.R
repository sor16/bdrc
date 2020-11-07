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
    plot_dat <- spread_draws(x,param,transformed)
    if(type=='trace'){
        if('param' %in% names(plot_dat)){
            params <- unique(plot_dat$param)
            plot_dat$param <- factor(sapply(plot_dat$param,get_parameter_levels,transformed),levels=get_parameter_levels(params),labels=params)
            plot_dat$chain <- factor(plot_dat$chain,levels=1:max(plot_dat$chain))
            plot_dat$param_expr <- sapply(plot_dat$param,get_param_expression,transformed)
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
            param_expr <- get_param_expression(param,transformed)
            p <- ggplot(plot_dat,aes_('',param),col="#BC3C29FF") +
                geom_line() +
                facet_wrap(~chain,scales='free') +
                xlab('Iteration') +
                ylab() +
                theme_classic() +
                theme(strip.background = element_blank(),
                      strip.text.x = element_text(size = 12))
        }

    }else if(type=='histogram'){

    }else if (type=='graph'){
        plot_dat <- spread_draws(x,param,transformed)
    }
    if(type %in% c('rating_curve',))
    spread_draws(x,)



}

predict_fun <- function(x,type){

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
    summary_fun()
}

plot.bgplm <- function(x,type='rating_curve',...){
}

predict.bgplm <- function(x,...){
}
