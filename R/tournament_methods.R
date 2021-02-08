#' Print tournament object
#'
#' Print the results of a tournament of model comparisons
#' @param x an object of class "tournament"
#' @seealso \code{\link{summary.torunament}} for summaries
#' @export
print.tournament <- function(x){
    cat(paste0('Tournament with winner ',class(x$winner)))
}

#' Print summary of tournament object
#'
#' Print the summary of a tournament of model comparisons
#' @param x an object of class "tournament"
#' @export
summary.tournament <- function(x){
    x$summary
}

#' Plot comparison of models in tournament
#'
#' @param x an object of class "tournament"
#' @param type character denoting the type of plot desired. One of "tree","DIC","residuals","f","sigma_eps","rating_curve","rating_curve_log"
#' @seealso \code{\link{summary.torunament}} for summaries
#' @importFrom gridExtra grid.arrange
#' @export
plot.tournament <- function(x,type="rating_curve"){
    if(type=="Deviance"){
        Deviance_post_dat <- lapply(x$contestants,function(m){
            data.frame(model=class(m),D=c(m$Deviance_posterior))
        })
        Deviance_post_dat <- do.call(rbind,Deviance_post_dat)
        p <- ggplot() +
             geom_boxplot(data=Deviance_post_dat,aes(x=model,y=D),width=0.4) +
             theme_classic() +
             xlab('Model') +
             ylab('Deviance')
    }else if(type=="residuals"){
        plot_list <- lapply(x$contestants,function(m){
                        plot(m,type='residuals',title=class(m))
                    })
        p <- do.call(gridExtra::grid.arrange,c(plot_list,ncol=2))
    }else if(type=="sigma_eps"){
        ylim_dat <- sapply(x$contestants,function(m){
                        if(grepl('0',class(m))){
                            c(m$param_summary['sigma_eps','lower'],m$param_summary['sigma_eps','upper'])
                        }else{
                            c(min(m$sigma_eps_summary$lower),max(m$sigma_eps_summary$upper))
                        }
                    })
        ylim_min <- max(c(0,min(ylim_dat[1,])))
        ylim_max <- max(ylim_dat[2,])
        plot_list <- lapply(x$contestants,function(m){
            plot(m,type="sigma_eps",title=class(m)) + ylim(ylim_min,ylim_max)
        })
        p <- do.call(gridExtra::grid.arrange,c(plot_list,ncol=2))
    }else if(type=="f"){
        ylim_dat <- sapply(x$contestants,function(m){
            if(!grepl('g',class(m))){
                c(m$param_summary['b','lower'],m$param_summary['b','upper'])
            }else{
                c(min(m$f_summary$lower),max(m$f_summary$upper))
            }
        })
        ylim_min <- min(ylim_dat[1,])
        ylim_max <- max(ylim_dat[2,])
        plot_list <- lapply(x$contestants,function(m){
            plot(m,type="f",title=class(m)) + ylim(ylim_min,ylim_max)
        })
        p <- do.call(gridExtra::grid.arrange,c(plot_list,ncol=2))
    }else if(type %in% c("rating_curve","rating_curve_log","rating_curve_mean","rating_curve_mean")){
        rc_type <- type
        transformed <- grepl('log',rc_type)
        rc_type <- gsub('_log','',type)
        plot_list <- lapply(x$contestants,function(m){
            plot(m,type=rc_type,transformed=transformed,title=class(m))
        })
        p <- do.call(gridExtra::grid.arrange,c(plot_list,ncol=2))
    }else{
        stop('type not recognized. Must be one of "Deviance","residuals","sigma_eps","f","rating_curve","rating_curve_log","rating_curve_mean","rating_curve_mean"')
    }
    return(p)
}
