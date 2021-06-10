#' Print tournament object
#'
#' Print the results of a tournament of model comparisons
#' @param x an object of class "tournament"
#' @param ... not used in this function
#' @seealso \code{\link{summary.tournament}} for summaries
#' @examples
#' \dontrun{
#' data(V316_river)
#' f <- Q~W
#' t_obj <- tournament(f,V316_river)
#' print(t_obj)
#' }
#' @export
print.tournament <- function(x,...){
    cat(paste0('Tournament with winner ',class(x$winner)))
}

#' Print summary of tournament object
#'
#' Print the summary of a tournament of model comparisons
#' @param object an object of class "tournament"
#' @param ... not used in this function
#' @examples
#' \dontrun{
#' data(V316_river)
#' f <- Q~W
#' t_obj <- tournament(f,V316_river)
#' summary(t_obj)
#' }
#' @export
summary.tournament <- function(object,...){
    object$summary
}

#' Autoplot - Comparison of models in tournament
#'
#' Compare the four models from the tournament object in different ways
#'
#' @param x an object of class "tournament"
#' @param type a character denoting what type of plot should be drawn. Possible types are
#' \itemize{
#'  \item{"deviance"}{ to plot the rating curve on original scale.}
#' }
#' @param ... not used in this function
#' @seealso \code{\link{summary.tournament}} for summaries
#' @examples
#' \dontrun{
#' data(V316_river)
#' f <- Q~W
#' t_obj <- tournament(f,V316_river)
#' autoplot(t_obj)
#' }
#' @importFrom ggplot2 ggplot geom_boxplot stat_boxplot geom_line geom_point xlab ylab
#' @importFrom rlang .data
#' @export
autoplot.tournament <- function(x,type='deviance',...){
    args <- list(...)
    legal_types <- c('deviance')
    if(!(type %in% legal_types)){
        stop(cat(paste('Type argument not recognized. Possible types are:\n - ',paste(legal_types,collapse='\n - '))))
    }else if(type=="deviance"){
        deviance_post_dat <- lapply(x$contestants,function(m){
            data.frame(model=class(m),D=c(m$Deviance_posterior))
        })
        deviance_post_dat <- do.call(rbind,deviance_post_dat)
        DIC_dat <- lapply(x$contestants,function(m){
            data.frame(model=class(m),DIC=c(m$DIC))
        })
        DIC_dat <- do.call(rbind,DIC_dat)
        p <- ggplot(data=deviance_post_dat,aes(x=.data$model,y=.data$D)) +
             geom_boxplot(size=.4,color="black",outlier.size=0.1,outlier.shape=21,outlier.fill="gray90",fill="gray90") +
             stat_boxplot(geom='errorbar',width=0.4) +
             geom_line(data=DIC_dat,aes(x=.data$model,y=.data$DIC,group=1),color='gray30') +
             geom_point(data=DIC_dat,aes(x=.data$model,y=.data$DIC),size=3,shape=23,fill='red2',color='black') +
             theme_bdrc() +
             xlab('') +
             ylab('Deviance')
    }
    return(p)
}


#' Plot comparison of models in tournament
#'
#' Compare the four models from the tournament object in multiple ways
#'
#' @param x an object of class "tournament"
#' @param type a character denoting what type of plot should be drawn. Possible types are
#' \itemize{
#'   \item{"deviance"}{ to plot the rating curve on original scale.}
#'   \item{"rating_curve"}{ to plot the rating curve on original scale.}
#'   \item{"rating_curve_mean"}{ to plot the rating curve on a log scale.}
#'   \item{"f"}{ to plot the power-law exponent}
#'   \item{"sigma_eps"}{ to plot the standard deviation on the data level}
#'   \item{"residuals"}{ to plot the log residuals}
#'  }
#' @param transformed a logical value indicating whether the quantity should be plotted on a transformed scale used during the Bayesian inference. Defaults to FALSE.
#' @param ... further arguments passed to other methods (currently unused).
#' @seealso \code{\link{summary.tournament}} for summaries
#' @examples
#' \dontrun{
#' data(V316_river)
#' f <- Q~W
#' t_obj <- tournament(f,V316_river)
#' plot(t_obj)
#' plot(t_obj,type='rating_curve_log')
#' plot(t_obj,type='deviance')
#' plot(t_obj,type='f')
#' plot(t_obj,type='sigma_eps')
#' plot(t_obj,type='residuals')
#' }
#' @importFrom ggplot2 ylim scale_y_continuous
#' @importFrom gridExtra grid.arrange
#' @importFrom grid grid.draw
#' @export
plot.tournament <- function(x,type='deviance',transformed=F,...){
    args <- list(...)
    legal_types <- c("deviance","rating_curve","rating_curve_mean","sigma_eps","f","residuals",'convergence_diagnostics')
    if(is.null(type) || type=='deviance'){
        p <- autoplot(x,transformed=transformed)
    }else if(type=="residuals"){
        res_dat <- sapply(names(x$contestants),function(m){
            max(abs(get_residuals_dat(x$contestants[[m]])[,c('r_median','r_lower','r_upper')]))
        })
        max_res <- max(res_dat)
        ylim_max <- 1.1*max_res
        plot_list <- lapply(x$contestants,function(m){
            suppressMessages(autoplot(m,type=type,title=class(m)) + ylim(-ylim_max,ylim_max))
                    })
        p <- do.call(arrangeGrob,c(plot_list,ncol=2))
    }else if(type=="sigma_eps"){
        ylim_dat <- sapply(x$contestants,function(m){
                        if(grepl('0',class(m))){
                            m$param_summary['sigma_eps','upper']
                        }else{
                            max(m$sigma_eps_summary$upper)
                        }
                    })
        ylim_min <- 0
        ylim_max <- 1.1*max(ylim_dat)
        plot_list <- lapply(x$contestants,function(m){
            suppressMessages(autoplot(m,type=type,title=class(m)) + scale_y_continuous(limits=c(ylim_min,ylim_max),expand=c(0,0)))
        })
        p <- do.call(arrangeGrob,c(plot_list,ncol=2))
    }else if(type=="f"){
        ylim_dat <- sapply(x$contestants,function(m){
            if(!grepl('g',class(m))){
                c(m$param_summary['b','lower'],m$param_summary['b','upper'])
            }else{
                c(min(m$f_summary$lower),max(m$f_summary$upper))
            }
        })
        ylim_min <- min(1,0.9*ylim_dat[1,])
        ylim_max <- max(3.5,1.1*ylim_dat[2,])
        plot_list <- lapply(x$contestants,function(m){
            supressMessages(autoplot(m,type=type,title=class(m)) + ylim(ylim_min,ylim_max))
        })
        p <- do.call(arrangeGrob,c(plot_list,ncol=2))
    }else if(type %in% c("rating_curve","rating_curve_mean")){
        plot_list <- lapply(x$contestants,function(m){
            autoplot(m,type=type,transformed=transformed,title=class(m))
        })
        p <- do.call(arrangeGrob,c(plot_list,ncol=2))
    }else if(type=='convergence_diagnostics'){
        plot_list <- lapply(x$contestants,function(m){
            plot_grob(m,type=type)
        })
        p <- do.call(arrangeGrob,c(plot_list,nrow=4))
    }else{
        stop(cat(paste0('type not recognized. Possible types are:',paste(legal_types,collapse='\n - '))))
    }
    if('ggplot' %in% class(p)){
        print(p)
    }else{
        grid.draw(p)
    }
}

#Idea, create graphical plot for tournament games
#df <- data.frame(round=c(1,1,1,1,2,2,3),
#                 game=c(1,1,2,2,1,1,1),
#                 model=c(c('plm0','plm','gplm0','gplm'),c('plm0','gplm'),'plm0'),
#                 xloc=c(rep(0,4),rep(1,2),2),
#                 yloc=c(seq(0,3,by=1),c(0.5,2.5),1.5))
# seg_dat <- data.frame(x1=c(c(xloc[1:4]),y))
#ggplot(df) + geom_text(aes(x=xloc,y=yloc,label=model),size=10) + theme_classic() + theme(line=element_blank(),text=element_blank())





