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
#' @param ... Arguments to be passed to other methods. The following arguments are supported:
#' \itemize{
#'   \item{\code{type}}{ a character denoting what type of plot should be drawn. Possible types are
#'                    \itemize{
#'                       \item{"Deviance"}{ to plot the rating curve on original scale.}
#'                    }}
#'   \item{\code{title}}{ title of the plot. Defaults to NULL, i.e. no title}
#' }
#' @seealso \code{\link{summary.tournament}} for summaries
#' @examples
#' \dontrun{
#' data(V316_river)
#' f <- Q~W
#' t_obj <- tournament(f,V316_river)
#' autoplot(t_obj)
#' }
#' @importFrom ggplot2 ggplot geom_boxplot xlab ylab
#' @importFrom rlang .data
#' @export
autoplot.tournament <- function(x,...){
    args <- list(...)
    if('type' %in% names(args)){
        type <- args$type
    }else{
        type <- 'Deviance'
    }
    if('title' %in% names(args)){
        title <- args$title
    }else{
        title <- NULL
    }
    legal_types <- c('Deviance')
    if(!(type %in% legal_types)){
        stop(cat(paste('Type argument not recognized. Possible types are:\n - ',paste(legal_types,collapse='\n - '))))
    }else if(type=="Deviance"){
        Deviance_post_dat <- lapply(x$contestants,function(m){
            data.frame(model=class(m),D=c(m$Deviance_posterior))
        })
        Deviance_post_dat <- do.call(rbind,Deviance_post_dat)
        DIC_dat <- lapply(x$contestants,function(m){
            data.frame(model=class(m),DIC=c(m$DIC))
        })
        DIC_dat <- do.call(rbind,DIC_dat)
        p <- ggplot(data=Deviance_post_dat,aes(x=.data$model,y=.data$D)) +
             geom_boxplot(size=.4,color="black",outlier.size=0.1,outlier.shape=21,outlier.fill="gray90",fill="gray90") +
             stat_boxplot(geom='errorbar') +
             geom_line(data=DIC_dat,aes(x=.data$model,y=.data$DIC,group=1),color='gray30') +
             geom_point(data=DIC_dat,aes(x=.data$model,y=.data$DIC),size=3,shape=23,fill='red2',color='black') +
             theme_bdrc() +
             xlab('Model') +
             ylab('Deviance')
    }
    #TODO - add num effective parameters plot
    return(p)
}


#' Plot comparison of models in tournament
#'
#' Compare the four models from the tournament object in multiple ways
#'
#' @param x an object of class "tournament"
#' @param ... Arguments to be passed to other methods. The following arguments are supported:
#' \itemize{
#'   \item{\code{type}}{ a character denoting what type of plot should be drawn. Possible types are
#'                    \itemize{
#'                       \item{"Deviance"}{ to plot the rating curve on original scale.}
#'                       \item{"rating_curve"}{ to plot the rating curve on original scale.}
#'                       \item{"rating_curve_mean"}{ to plot the rating curve on a log scale.}
#'                       \item{"f"}{ to plot the power-law exponent}
#'                       \item{"sigma_eps"}{ to plot the standard deviation on the data level}
#'                       \item{"residuals"}{ to plot the log residuals}
#'                    }}
#'   \item{\code{transformed}}{ a logical value indicating whether the quantity should be plotted on a transformed scale used during the Bayesian inference. Defaults to FALSE.}
#'   \item{\code{title}}{ title of the plot. Defaults to NULL, i.e. no title}
#' }
#' @seealso \code{\link{summary.tournament}} for summaries
#' @examples
#' \dontrun{
#' data(V316_river)
#' f <- Q~W
#' t_obj <- tournament(f,V316_river)
#' plot(t_obj)
#' plot(t_obj,type='rating_curve_log')
#' plot(t_obj,type='Deviance')
#' plot(t_obj,type='f')
#' plot(t_obj,type='sigma_eps')
#' plot(t_obj,type='residuals')
#' }
#' @importFrom ggplot2 ylim
#' @importFrom gridExtra grid.arrange
#' @importFrom grid grid.draw
#' @export
plot.tournament <- function(x,...){
    args <- list(...)
    if('transformed' %in% names(args)){
        transformed <- args$transformed
    }else{
        transformed <- F
    }
    legal_types <- c("Deviance","rating_curve","rating_curve_mean","sigma_eps","f","residuals")
    if(is.null(args$type) || args$type=='Deviance'){
        p <- autoplot(x,...)
    }else if(args$type=="residuals"){
        plot_list <- lapply(x$contestants,function(m){
                        autoplot(m,type=args$type,title=class(m))
                    })
        p <- do.call(grid.arrange,c(plot_list,ncol=2))
    }else if(args$type=="sigma_eps"){
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
            autoplot(m,type=args$type,title=class(m)) + ylim(ylim_min,ylim_max)
        })
        p <- do.call(grid.arrange,c(plot_list,ncol=2))
    }else if(args$type=="f"){
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
            autoplot(m,type=args$type,title=class(m)) + ylim(ylim_min,ylim_max)
        })
        p <- do.call(grid.arrange,c(plot_list,ncol=2))
    }else if(args$type %in% c("rating_curve","rating_curve_mean")){
        plot_list <- lapply(x$contestants,function(m){
            autoplot(m,type=args$type,transformed=args$transformed,title=class(m))
        })
        p <- do.call(grid.arrange,c(plot_list,ncol=2))
    }else{
        stop(cat(paste0('type not recognized. Possible types are:',paste(legal_types,collapse='\n - '))))
    }
    if('ggplot' %in% class(p)){
        print(p)
    }else{
        grid.draw(p)
    }
}
