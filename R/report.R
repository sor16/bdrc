#' @importFrom stats quantile predict
#' @importFrom ggplot2 autoplot scale_x_continuous scale_y_continuous
#' @importFrom gridExtra tableGrob ttheme_minimal
#' @importFrom grid textGrob
get_report_components <- function(x,type=1){
    error_msg1 <- 'Please provide a single object of types tournament, gplm, gplm0, plm or plm0.'
    error_msg2 <- 'It is only possible to produce a type 1 report for a single model object of type gplm, gplm0, plm or plm0.'
    legal_types <- c('gplm','gplm0','plm','plm0','tournament')
    if(!(type %in% c(1,2))){
        stop('Please input an integer value of 1 or 2 to indicate which type of report is to be produce.')
    }
    if(!(class(x) %in% legal_types)){
        stop(error_msg1)
    }
    if(inherits(x,'tournament')){
        if(type==1){
            m_obj <- list(x$winner)
            names(m_obj) <- class(x$winner)
        }else{
            m_obj <- x$contestants
            t_obj <- x
        }
    }else{
        if(type==2){
            stop('It is only possible to produce a type 1 report for a single model object of type gplm, gplm0, plm or plm0.')
        }
        m_obj <- list(x)
        names(m_obj) <- class(x)
    }
    output_list <- list()
    main_plot_types <- c('rating_curve','residuals','sigma_eps','f')
    if(type==1){
        output_list$main_page_plots <- plot_grob(m_obj[[1]],type='panel')
        output_list$main_page_table <- lapply(m_obj,function(m){
                                            param <- get_param_names(class(m),m$run_info$c_param)
                                            table <- rbind(m$param_summary[,c('lower','median','upper')],c(m$Deviance_summary))
                                            names(table) <- paste0(names(table),c('-2.5%','-50%','-97.5%'))
                                            row.names(table) <- c(sapply(1:length(param),get_param_expression),"Deviance")
                                            table <- format(round(table,digits=3),nsmall=3)
                                            tableGrob(table,theme=ttheme_minimal(rowhead=list(fg_params=list(parse=TRUE))))
                                        })
        p_mat <- predict(m_obj[[1]],wide=TRUE)
        num_pages <- floor(nrow(p_mat)/41)+1
        output_list$p_mat_list <- lapply(1:num_pages,function(i){
            idx <- if(num_pages==1) 1:nrow(p_mat) else if(i==num_pages) ((num_pages-1)*40+1):nrow(p_mat) else ((i-1)*40+1):((i)*40)
            tableGrob(format(round(p_mat[idx,],digits=3),nsmall=3),
                      theme=ttheme_minimal(core=list(bg_params = list(fill = c("#F7FBFF","#DEEBF7"), col=NA),fg_params=list(fontface=3)),
                                           colhead=list(fg_params=list(col="black",fontface=2L)),
                                           rowhead=list(fg_params=list(col="black",fontface=2L)),
                                           base_size = 10))
        })
        output_list$obj_class <- class(m_obj[[1]])
    }else{
        output_list$main_page_plots <- plot_tournament_grob(t_obj)
        output_list$main_page_table <- lapply(m_obj,function(m){
                                            param <- get_param_names(class(m),m$run_info$c_param)
                                            table <- rbind(m$param_summary,data.frame(m$Deviance_summary,r_hat=NA,eff_n_samples=NA))
                                            table[,c('lower','median','upper','r_hat')] <- format(round(table[,c('lower','median','upper','r_hat')],digits=3),nsmall=3)
                                            row.names(table) <- c(sapply(1:length(param),get_param_expression),"Deviance")
                                            table['Deviance',c('eff_n_samples','r_hat')] <- ''
                                            names(table) <- c(paste0(names(table[,c('lower','median','upper')]),c('-2.5%','-50%','-97.5%')),'eff_n_samples','r_hat')
                                            tableGrob(table,theme=ttheme_minimal(rowhead=list(fg_params=list(parse=TRUE))))
                                        })
        output_list$mcmc_hist_list <- lapply(m_obj,function(m){
            params <- get_param_names(class(m),m$run_info$c_param)
            hist_plot_list <- lapply(1:length(params), function(j){
                autoplot(m,type='histogram',param=params[j],transformed = TRUE,title='')
            })
        })
        mcmc_table <- lapply(m_obj,function(m){
            data.frame(nr_iter=m$run_info$nr_iter,
                       nr_chains=m$run_info$num_chains,
                       burnin=m$run_info$burnin,
                       thin=m$run_info$thin,
                       eff_num_param=format(m$effective_num_param_WAIC,digits=2),
                       acceptance_rate=format(m$acceptance_rate,digits=2))
        })
        mcmc_table <- t(do.call('rbind',mcmc_table))
        output_list$mcmc_table <- tableGrob(mcmc_table,theme=ttheme_minimal())
        tour_table <- t_obj$summary
        col_idx1 <- (names(tour_table) %in% c('round','game','model','WAIC','DIC','Post_prob','winner'))
        col_idx2 <- !(names(tour_table) %in% c('round','game','model','winner'))
        tour_table[,col_idx2] <- round(tour_table[,col_idx2],digits=2)
        output_list$tour_table <- tableGrob(tour_table[,col_idx1],theme=ttheme_minimal(),rows=NULL)
        output_list$dev_boxplot <- autoplot(t_obj,type='deviance')
        output_list$tournament_results <- plot_tournament_grob(t_obj,type='tournament_results')
        output_list$conv_diag_plots <- lapply(t_obj$contestants,function(x){
            plot_grob(x,type='convergence_diagnostics')
        })
    }
    return(output_list)
}

#' @importFrom gridExtra arrangeGrob
#' @importFrom grid textGrob
get_report_pages_fun <- function(x,type=1){
    report_components <- get_report_components(x,type=type)
    if(type==1){
        main_page_plot <- arrangeGrob(report_components$main_page_plots,
                                      report_components$main_page_table[[1]],
                                      nrow=2,
                                      as.table=TRUE,
                                      heights=c(5,3),
                                      top=textGrob(report_components$obj_class,gp=gpar(fontsize=22,facetype='bold')))
        predict_mat <- lapply(report_components$p_mat_list,function(p){
            arrangeGrob(p,
                        nrow=1,
                        as.table=TRUE,
                        heights=c(1),
                        top=textGrob(paste0('Tabular Rating Curve - ',report_components$obj_class),gp=gpar(fontsize=22,facetype='bold')))
        })
        report_pages <- c(list(main_page_plot),predict_mat)
    }else{
        main_page_plots <- lapply(names(report_components$main_page_plots),function(m){
                                  arrangeGrob(report_components$main_page_plots[[m]],
                                              report_components$main_page_table[[m]],
                                              nrow=2,
                                              as.table=TRUE,
                                              heights=c(5,3),
                                              top=textGrob(m,gp=gpar(fontsize=22,facetype='bold')))
                           })
        tournament_page <- arrangeGrob(arrangeGrob(report_components$dev_boxplot,
                                                   arrangeGrob(report_components$mcmc_table,
                                                               report_components$tour_table,nrow=2),
                                                   ncol=2),
                                       report_components$tournament_results,
                                       nrow=2,
                                       as.table=TRUE,
                                       heights=c(1,1),
                                       top=textGrob('Model comparison and MCMC diagnostics',gp=gpar(fontsize=22,facetype='bold')))
        convergence_page <- arrangeGrob(grobs=report_components$conv_diag_plots,
                                        nrow=4,
                                        as.table=TRUE)
        histogram_pages <- lapply(names(report_components$main_page_plots), function(m) {
                                  arrangeGrob(arrangeGrob(grobs=report_components$mcmc_hist_list[[m]],nrow=4,ncol=3),
                                              top=textGrob(paste0('Estimated parameters of ',m),gp=gpar(fontsize=20,facetype='bold')))
                            })
        report_pages <- c(main_page_plots,list(tournament_page),list(convergence_page),histogram_pages)
    }
    return(report_pages)
}

#' @importFrom utils askYesNo
#' @importFrom grDevices pdf dev.off
save_report <- function(report_pages,path=NULL,paper='a4',width=9,height=11){
    if(is.null(path)){
        path <- paste0(getwd(),'/report.pdf')
    }
    if(interactive()){
        answer <- askYesNo(paste0('Do you wish to save the report as a pdf file at the following location:\n ',path,'?'))
    }else{
        stop('Unable to ask permission for writing the report to the file system. get_report() must be used in an interactive R session')
    }
    if(answer==FALSE | is.na(answer)) stop('Permission to store a pdf file was not granted')
    pdf(file=path,paper=paper,width=width,height=height)
    for(i in 1:length(report_pages)){
        grid.arrange(report_pages[[i]],as.table=TRUE)
    }
    invisible(dev.off())
    message(paste0('The report was saved at the following location:\n',path))
}

#### S3 methods
#' Report pages for a discharge rating curve or tournament
#'
#' Get a list of the pages of a report on a discharge rating curve model or tournament
#' @param x an object of class "tournament", "plm0", "plm", "gplm0" or "gplm".
#' @param type an integer denoting what type of report is to be produced. Defaults to type 1. Possible types are
#'                    \itemize{
#'                       \item{1 - }{produces a report displaying the results of the model (winning model if a tournament provided). The first page contains a panel of four plots and a summary of the posterior distributions of the parameters. On the second page a tabular prediction of discharge on an equally spaced grid of stages is displayed. This prediction table can span multiple pages.}
#'                       \item{2 - }{produces a ten page report and is only permissible for objects of class "tournament". The first four pages contain a panel of four plots and a summary of the posterior distributions of the parameters for each of the four models in the tournament, the fifth page shows model comparison plots and tables, the sixth page convergence diagnostics plots, and the final four pages shows the histograms of the parameters in each of the four models.}
#'                    }
#' @param ... further arguments passed to other methods (currently unused).
#' @return A list of objects of type "grob" that correspond to the pages in a rating curve report.
#' @seealso \code{\link{tournament}} for running a tournament,\code{\link{summary.tournament}} for summaries and \code{\link{get_report}} for generating and saving a report of a tournament object.
#' @examples
#' \donttest{
#' data(krokfors)
#' set.seed(1)
#' plm0.fit <- plm0(formula=Q~W,data=krokfors,num_cores=2)
#' plm0_pages <- get_report_pages(plm0.fit)
#' }
#' @export
get_report_pages <- function(x,type=1,...) UseMethod("get_report_pages")

#' @describeIn get_report_pages Get report pages for plm0 model object
#' @export
get_report_pages.plm0 <- function(x,type=1,...){
    get_report_pages_fun(x,type=type)
}

#' @describeIn get_report_pages Get report pages for plm model object
#' @export
get_report_pages.plm <- function(x,type=1,...){
    get_report_pages_fun(x,type=type)
}

#' @describeIn get_report_pages Get report pages for gplm0 model object
#' @export
get_report_pages.gplm0 <- function(x,type=1,...){
    get_report_pages_fun(x,type=type)
}

#' @describeIn get_report_pages Get report pages for gplm model object
#' @export
get_report_pages.gplm <- function(x,type=1,...){
    get_report_pages_fun(x,type=type)
}

#' @describeIn get_report_pages Get report pages for discharge rating curve tournament model object
#' @export
get_report_pages.tournament <- function(x,type=1,...){
    get_report_pages_fun(x,type=type)
}

#' Report for a discharge rating curve or tournament
#'
#' Save a pdf file with a report of a discharge rating curve object or tournament.
#' @param x an object of class "tournament", "plm0", "plm", "gplm0" or "gplm".
#' @param path file path to which the pdf file of the report is saved. If NULL, the current working directory is used.
#' @param type an integer denoting what type of report is to be produced. Defaults to type 1. Only type 1 is permissible for an object of class "plm0", "plm", "gplm0" or "gplm". Possible types are
#'                    \itemize{
#'                       \item{1 - }{produces a report displaying the results of the model (winning model if a tournament provided). The first page contains a panel of four plots and a summary of the posterior distributions of the parameters. On the second page a tabular prediction of discharge on an equally spaced grid of stages is displayed. This prediction table can span multiple pages.}
#'                       \item{2 - }{produces a ten page report and is only permissible for objects of class "tournament". The first four pages contain a panel of four plots and a summary of the posterior distributions of the parameters for each of the four models in the tournament, the fifth page shows model comparison plots and tables, the sixth page convergence diagnostics plots, and the final four pages shows the histograms of the parameters in each of the four models.}
#'                    }
#' @param ... further arguments passed to other methods (currently unused).
#' @details  This function can only be used in an interactive R session as it asks permission from the user to write to their file system.
#' @return No return value, called for side effects
#' @seealso \code{\link{get_report}} for generating and saving a report.
#' @examples
#' \donttest{
#' data(krokfors)
#' set.seed(1)
#' plm0.fit <- plm0(formula=Q~W,data=krokfors,num_cores=2)
#' }
#' \dontrun{
#' get_report(plm0.fit)
#' }
#' @export
get_report <- function(x,path=NULL,type=1,...) UseMethod("get_report")

#' @describeIn get_report Get report for plm0 model object
#' @export
get_report.plm0 <- function(x,path=NULL,type=1,...){
    report_pages <- get_report_pages_fun(x,type=type)
    save_report(report_pages,path=path,...)
}

#' @describeIn get_report Get report for plm model object
#' @export
get_report.plm <- function(x,path=NULL,type=1,...){
    report_pages <- get_report_pages_fun(x,type=type)
    save_report(report_pages,path=path,...)
}

#' @describeIn get_report Get report for gplm0 model object
#' @export
get_report.gplm0 <- function(x,path=NULL,type=1,...){
    report_pages <- get_report_pages_fun(x,type=type)
    save_report(report_pages,path=path,...)
}

#' @describeIn get_report Get report for gplm
#' @export
get_report.gplm <- function(x,path=NULL,type=1,...){
    report_pages <- get_report_pages_fun(x,type=type)
    save_report(report_pages,path=path,...)
}

#' @describeIn get_report Get report for discharge rating curve tournament
#' @export
get_report.tournament <- function(x,path=NULL,type=1,...){
    report_pages <- get_report_pages_fun(x,type=type)
    save_report(report_pages,path=path,...)
}
