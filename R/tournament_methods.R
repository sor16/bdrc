#' @importFrom ggplot2 autoplot ggplot geom_segment geom_text scale_colour_manual theme_classic geom_label theme unit element_blank
#' @importFrom gridExtra arrangeGrob
plot_tournament_fun <- function(x,type='deviance'){
    if(type=='deviance'){
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
    }else if(type=='tournament_results'){
        loc_pts <- data.frame(x=c(seq(0,3),0.5,2.5,1.5),
                              y=c(rep(0,4),1,1,2),
                              xend=c(0.5,0.5,2.5,2.5,1.5,1.5,NA),
                              yend=c(rep(1,4),2,2,NA),
                              model=c(sapply(x$contestants,function(m)class(m)),
                                      x$summary$model[5:6],
                                      paste0('Tournament winner  =>  ',class(x$winner),
                                             paste0(rep(' ',20),collapse=' '))))
        prob_dat <- data.frame(P=round(x$summary$P,digits=3),
                               winner=x$summary$winner,
                               x=c(loc_pts$x[1:4],0.8*(loc_pts$x[5:6]-1.5)+1.5),
                               y=loc_pts$y[1:6]+0.5)
        game_results <- ggplot() +
            geom_segment(data=loc_pts[1:6,],aes(x=x,y=y,xend=xend,yend=yend)) +
            geom_text(data=prob_dat, aes(x=x,y=y,label=P,color=winner,size=7)) +
            scale_colour_manual(values = c("red", "green3")) +
            theme_classic() +
            geom_label(data=loc_pts[5:7,],
                       label=loc_pts$model[5:7],
                       x=loc_pts$x[5:7],
                       y=loc_pts$y[5:7],
                       label.padding=unit(0.5,"lines"),
                       label.size=0,
                       color="Black",
                       fill="white",
                       size=6) +
            theme(line=element_blank(),
                  text=element_blank(),
                  plot.margin=unit(c(0,1,-0.5,3),"cm"),
                  legend.position="none")
        grob_list <- plot_tournament_grob(t_obj,type='residuals')
        p <- arrangeGrob(game_results,arrangeGrob(grobs=grob_list$grobs,ncol=4),nrow=2,heights=c(1,1))
    }
    return(p)
}


#' @importFrom ggplot2 autoplot
#' @importFrom gridExtra arrangeGrob
plot_tournament_grob <- function(x,type='panel',transformed=F){
    ylim_types <- c('sigma_eps','f','panel','tournament_results','residuals')
    if(type%in%ylim_types){
        ylim_dat <- lapply(x$contestants,function(m){
            if(grepl('0',class(m))){
                sig_ylim <- c( min=0, max=1.1*m$param_summary['sigma_eps','upper'] )
            }else{
                sig_ylim <- c( min=0, max=1.1*max(m$sigma_eps_summary$upper) )
            }
            if(grepl('g',class(m))){
                f_ylim <- c( min=min(1,0.9*min(m$f_summary$lower)), max=max(3.5,1.1*max(m$f_summary$upper)) )
            }else{
                f_ylim <- c( min=min(1,0.9*m$param_summary['b','lower']), max=max(3.5,1.1*m$param_summary['b','upper']) )
            }
            max_res <- 1.1*max(abs(get_residuals_dat(m)[,c('r_median','r_lower','r_upper')]))
            data.frame(sigma_eps_min=sig_ylim[1],sigma_eps_max=sig_ylim[2],f_min=f_ylim[1],f_max=f_ylim[2],residuals_min=-max_res,residuals_max=max_res)
        })
        ylim_dat <- do.call('rbind',ylim_dat)
        ylim_dat <- sapply(colnames(ylim_dat),function(col){
            if(grepl('min',col)){
                min(ylim_dat[,col])
            }else{
                max(ylim_dat[,col])
            }
        })
        ylim_dat <- c(ylim_dat,rating_curve_min=NA,rating_curve_max=NA)
    }
    grob_list <- list()
    if(type=="residuals"){
        plot_list <- lapply(x$contestants,function(m){
            autoplot(m,type=type,title=class(m),ylim=ylim_dat[c('residuals_min','residuals_max')])
        })
        grob_list <- do.call(arrangeGrob,c(plot_list,ncol=2))
    }else if(type=="sigma_eps"){
        plot_list <- lapply(x$contestants,function(m){
            autoplot(m,type=type,title=class(m),ylim=ylim_dat[c('sigma_eps_min','sigma_eps_max')])
        })
        grob_list <- do.call(arrangeGrob,c(plot_list,ncol=2))
    }else if(type=="f"){
        plot_list <- lapply(x$contestants,function(m){
            autoplot(m,type=type,title=class(m),ylim=ylim_dat[c('f_min','f_max')])
        })
        grob_list <- do.call(arrangeGrob,c(plot_list,ncol=2))
    }else if(type %in% c("rating_curve","rating_curve_mean")){
        plot_list <- lapply(x$contestants,function(m){
            autoplot(m,type=type,transformed=transformed,title=class(m))
        })
        grob_list <- do.call(arrangeGrob,c(plot_list,ncol=2))
    }else if(type=='convergence_diagnostics'){
        plot_list <- lapply(x$contestants,function(m){
            plot_grob(m,type=type)
        })
        grob_list <- do.call(arrangeGrob,c(plot_list,nrow=4))
    }else if(type=='panel'){
        panel_types <- c('rating_curve','residuals','f','sigma_eps')
        grob_list <- lapply(x$contestants,function(m){
            plot_list <- lapply(panel_types,function(ty){
                ylim_vec <- ylim_dat[c(paste0(ty,'_min'),paste0(ty,'_max'))]
                plot_fun(m,type=ty,transformed=transformed,param=NULL,ylim=ylim_vec)
            })
            p <- do.call('arrangeGrob',c(plot_list,ncol=round(sqrt(length(panel_types)))))
        })
        names(grob_list) <- sapply(x$contestants,class)
    }else if(type=='tournament_results'){
        grob_list <- plot_tournament_fun(x,type='tournament_results')
    }
    return(grob_list)
}

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
        p <- plot_tournament_fun(x,type=type)
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
#' @importFrom grid grid.draw
#' @importFrom gridExtra grid.arrange
#' @export
plot.tournament <- function(x,type='deviance',transformed=F,...){
    args <- list(...)
    legal_types <- c("deviance","rating_curve","rating_curve_mean","sigma_eps","f","residuals",'convergence_diagnostics','panel','tournament_results')
    if(is.null(type) || type=='deviance'){
        p <- autoplot(x,type=type)
    }else if(type%in%legal_types){
        p <- plot_tournament_grob(x,type=type,transformed=transformed,...)
    }else{
        stop(cat(paste0('type not recognized. Possible types are:',paste(legal_types,collapse='\n - '))))
    }
    if('ggplot' %in% class(p)){
        print(p)
    }else{
        if(type=='panel'){
            grid.draw(p[[class(x$winner)]])
        }else{
            grid.draw(p)
        }
    }
}
