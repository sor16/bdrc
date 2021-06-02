### TODO: Option in the predict function?? I.e. an help function for predict method. Ath aukastafi. Titill á töflu
#' @importFrom gridExtra tableGrob ttheme_minimal
#' @importFrom stats predict
predict_matrix <- function(x){
    # c_param <- if(is.null(x$run_info$c_param)) median(x$c_posterior) else x$run_info$c_param
    # c_param <- ceiling(c_param*100)/100
    min_rc_h <- min(x$rating_curve[,'h'])
    min_rc_h <- ceiling(min_rc_h*100)/100
    grid_max <- x$run_info$h_max
    p_mat <- predict(x,newdata=seq(min_rc_h,grid_max,by=0.01),wide=TRUE)
    p_mat <- format(round(p_mat,digits=3),nsmall=3)
    p_mat_grob <- tableGrob(p_mat,theme=ttheme_minimal(core=list(bg_params = list(fill = c("#F7FBFF","#DEEBF7"), col=NA),
                                                       fg_params=list(fontface=3)),
                                                       colhead=list(fg_params=list(col="black",fontface=2L)),
                                                       rowhead=list(fg_params=list(col="black",fontface=2L)),
                                                       base_size = 10))
    # split_row <- 6
    # n_pages <- ceiling(nrow(p_mat)/split_row)
    # if(n_pages>1){
    #     split_points <- seq(split_row,floor(nrow(p_mat)/split_row)*split_row,by=split_row)
    #     idx_vectors <- unname(split(seq(1,nrow(p_mat)), cumsum(seq_along(seq(1,nrow(p_mat))) %in% split_points)))
    #     p_mat_list <- lapply(idx_vectors,function(x){
    #                       format(p_mat[x,],nsmall=3)
    #                   })
    # }
    # p_mat_list <- lapply(p_mat_list,function(x){
    #                   tableGrob(x,theme=ttheme_minimal(core=list(bg_params = list(fill = blues9[1:2], col=NA),fg_params=list(fontface=3)),
    #                                                        colhead=list(fg_params=list(col="black",fontface=2L)),
    #                                                        rowhead=list(fg_params=list(col="black",fontface=2L))))
    #                })
    return(p_mat_grob)
}


#' @importFrom stats quantile
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
    if(class(x)=='tournament'){
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
    h_dat <- m_obj[[1]]$data[[all.vars(m_obj[[1]]$formula)[2]]]
    q_dat <- m_obj[[1]]$data[[all.vars(m_obj[[1]]$formula)[1]]]
    posterior_list <- lapply(m_obj,function(m){
                            if(is.null(m$run_info$c_param)){
                                post_dat <- spread_draws(m,'rating_curve','f','sigma_eps','c')
                            }else{
                                post_dat <- spread_draws(m,'rating_curve','f','sigma_eps')
                                post_dat$c <- m$run_info$c_param
                            }
                            post_dat[post_dat$h>=min(h_dat) & post_dat$h<=max(h_dat),]
                      })
    res_dat <- sapply(names(m_obj),function(x){
                    max(abs(get_residuals_dat(m_obj[[x]])[,c('r_median','r_lower','r_upper')]))
               })
    max_res <- max(res_dat)
    lim_list <- lapply(posterior_list,function(df){
                    data.frame(rating_curve_x_min=0,rating_curve_x_max=1.01*max(quantile(df[df$h==max(df$h),]$rating_curve,0.975),max(q_dat)),
                               rating_curve_y_min=median(df$c),rating_curve_y_max=1.01*max(df$h)-0.01*min(df$h),
                               residuals_y_min=1.1*(-max_res),residuals_y_max=1.1*max_res,
                               residuals_x_min=NA,residuals_x_max=NA,
                               sigma_eps_x_min=min(df$h),sigma_eps_x_max=max(df$h),
                               sigma_eps_y_min=0,sigma_eps_y_max=max(df$sigma_eps),
                               f_x_min=min(df$h),f_x_max=max(df$h),
                               f_y_min=min(df$f,1),f_y_max=max(df$f,3.5))
                })
    lim_dat <- do.call('rbind',lim_list)
    output_list <- list()
    main_plot_types <- c('rating_curve','residuals','sigma_eps','f')
    output_list$main_page_plots <- lapply(m_obj,function(m){
                                        pt_plot_list <- lapply(main_plot_types,function(pt){
                                                            autoplot(m,type=pt) +
                                                                scale_x_continuous(limits = c(min(lim_dat[[paste0(pt,'_x_min')]]),max(lim_dat[[paste0(pt,'_x_max')]]))) +
                                                                scale_y_continuous(limits = c(min(lim_dat[[paste0(pt,'_y_min')]]),max(lim_dat[[paste0(pt,'_y_max')]])))
                                                        })
                                        do.call('arrangeGrob',pt_plot_list)
                                    })
    if(type==1){
        output_list$main_page_table <- lapply(m_obj,function(m){
                                            param <- get_param_names(class(m),m$run_info$c_param)
                                            table <- rbind(m$param_summary[,c('lower','median','upper')],c(m$Deviance_summary))
                                            names(table) <- paste0(names(table),c('-2.5%','-50%','-97.5%'))
                                            row.names(table) <- c(sapply(1:length(param),get_param_expression),"Deviance")
                                            table <- format(round(table,digits=3),nsmall=3)
                                            tableGrob(table,theme=ttheme_minimal(rowhead=list(fg_params=list(parse=TRUE))))
                                        })
        output_list$p_mat_list <- predict_matrix(x)

        output_list$obj_class <- class(m_obj[[1]])
    }else{
        output_list$main_page_table <- lapply(m_obj,function(m){
                                            param <- get_param_names(class(m),m$run_info$c_param)
                                            table <- rbind(m$param_summary,data.frame(m$Deviance_summary,r_hat=NA,n_eff_samples=NA))
                                            table[,c('lower','median','upper','r_hat')] <- format(round(table[,c('lower','median','upper','r_hat')],digits=3),nsmall=3)
                                            row.names(table) <- c(sapply(1:length(param),get_param_expression),"Deviance")
                                            table['Deviance',c('n_eff_samples','r_hat')] <- ''
                                            names(table) <- c(paste0(names(table[,c('lower','median','upper')]),c('-2.5%','-50%','-97.5%')),'num_eff_samples','r_hat')
                                            tableGrob(table,theme=ttheme_minimal(rowhead=list(fg_params=list(parse=TRUE))))
                                        })
        output_list$mcmc_hist_list <- lapply(m_obj,function(m){
            params <- get_param_names(class(m),m$run_info$c_param)
            hist_plot_list <- lapply(1:length(params), function(j){
                autoplot(m,type='histogram',param=params[j],transformed = T)
            })
        })
        mcmc_table <- lapply(m_obj,function(m){
            data.frame(nr_iter=m$run_info$nr_iter,
                       nr_chains=m$run_info$num_chains,
                       burnin=m$run_info$burnin,
                       thin=m$run_info$thin,
                       nr_eff_param=format(m$num_effective_param,digits=2),
                       acceptance_rate=format(m$acceptance_rate,digits=2))
        })
        mcmc_table <- t(do.call('rbind',mcmc_table))
        output_list$mcmc_table <- tableGrob(mcmc_table,theme=ttheme_minimal())
        tour_table <- t_obj$summary[c("round","game","model","DIC","P","winner")]
        tour_table[c('DIC','P')] <- round(tour_table[c('DIC','P')],digits=2)
        output_list$tour_table <- tableGrob(tour_table,theme=ttheme_minimal(),rows=NULL)
        output_list$dev_boxplot <- autoplot(t_obj,type='deviance')
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
        main_page_plot <- arrangeGrob(report_components$main_page_plots[[1]],
                                      report_components$main_page_table[[1]],
                                      nrow=2,
                                      as.table=TRUE,
                                      heights=c(5,3),
                                      top=textGrob(report_components$obj_class,gp=gpar(fontsize=22,facetype='bold')))
        predict_mat <- arrangeGrob(report_components$p_mat_list,     # thu ert her !!!
                                   nrow=1,
                                   as.table=TRUE,
                                   heights=c(1),
                                   #top=textGrob(paste0('Rating curve predictions in cubic meters per second as a function of stage in meters.\n The stage increments by a decimeter per row and the columns represent centimeter increments.',report_components$obj_class),gp=gpar(fontsize=12,facetype='bold')))
                                   top=textGrob(paste0('Tabular Rating Curve - ',report_components$obj_class),gp=gpar(fontsize=22,facetype='bold')))
        report_pages <- list(main_page_plot,predict_mat)
    }else{
        main_page_plots <- lapply(names(report_components$main_page_plots),function(m){
                                  arrangeGrob(report_components$main_page_plots[[m]],
                                              report_components$main_page_table[[m]],
                                              nrow=2,
                                              as.table=TRUE,
                                              heights=c(5,3),
                                              top=textGrob(m,gp=gpar(fontsize=22,facetype='bold')))
                           })
        tournament_page <- arrangeGrob(arrangeGrob(report_components$mcmc_table,report_components$tour_table,nrow=1),
                                       arrangeGrob(report_components$dev_boxplot,ncol=2),
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

#' @importFrom grDevices pdf dev.off
save_report <- function(report_pages,path=NULL,paper='a4',width=9,height=11){
    if(is.null(path)){
        path <- 'report.pdf'
    }
    pdf(file=path,paper=paper,width=width,height=height)
    for(i in 1:length(report_pages)){
        grid.arrange(report_pages[[i]],as.table=T)
    }
    dev.off()
}

#### S3 methods


#' Get report pages
#'
#' Get a list of the grob objects used by \code{\link{get_report}} to create a pdf report.
#' @param x an object of class "tournament", "plm0", "plm", "gplm0" or "gplm".
#' @param type an integer denoting what type of report is to be produced. Defaults to type 1. Only type 1 is permissible for an object of class "plm0", "plm", "gplm0" or "gplm". Possible types are
#'                    \itemize{
#'                       \item{1}{ produces a list of two grob objects, the first grob object consists of a plot collage and parameter summary table, the second a tabular rating curve from the results of a model object, or if input class is "tournament", then the most approriate model as determined by \code{\link{tournament}}.}
#'                       \item{2}{ produces a list of ten grob objects, the first four grob objects consist of plot collage and parameter summary table for the four models, the fifth object consists of model comparison plots and tables, the sixth consists of convergence diagnostics plots, and the final four grob objects consist of histogram plots of the estmated parameters for the four models.}
#'                    }
#' @param ... further arguments passed to other methods (currently unused).
#' @seealso \code{\link{tournament}} for running a tournament,\code{\link{summary.tournament}} for summaries.
#' @examples
#' \dontrun{
#' data(bunnerviken)
#' f <- Q~W
#' t_obj <- tournament(f,bunnerviken)
#' grob_object_list <- get_report_pages(t_obj)
#' }
#' @export
get_report_pages <- function(x,type=1,...) UseMethod("get_report_pages")


#' Get plm0 fit report pages
#'
#' Get a list of the grob objects used by \code{\link{get_report}} to create a pdf report.
#' @param x an object of class "plm0".
#' @param type an integer to determine whether to get the pages from report type 1 or 2. Only type 1 is permissible for an object of class "plm0" which produces a list of two grob objects, the first grob object consists of a plot collage and parameter summary table, the second a tabular rating curve.
#' @param ... further arguments passed to other methods (currently unused).
#' @seealso \code{\link{plm0}} for fitting the plm0 model,\code{\link{summary.plm0}} for summaries, \code{\link{predict.plm0}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{plot.plm0}} to help visualize the full posterior distributions.
#' @examples
#' \dontrun{
#' data(halla)
#' f <- Q~W
#' plm0.fit <- plm0(f,halla)
#' grob_object_list <- get_report_pages(plm0.fit)
#' }
#' @export
get_report_pages.plm0 <- function(x,type=1,...){
    get_report_pages_fun(x,type=type)
}


#' Get plm fit report pages
#'
#' Get a list of the grob objects used to by \code{\link{get_report}} create a pdf report.
#' @param x an object of class "plm".
#' @param type an integer to determine whether to get the pages from report type 1 or 2. Only type 1 is permissible for an object of class "plm" which produces a list of two grob objects, the first grob object consists of a plot collage and parameter summary table, the second a tabular rating curve.
#' @param ... further arguments passed to other methods (currently unused).
#' @seealso \code{\link{plm}} for fitting the plm model,\code{\link{summary.plm}} for summaries, \code{\link{predict.plm}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{plot.plm}} to help visualize the full posterior distributions.
#' @examples
#' \dontrun{
#' data(lisjobacken)
#' f <- Q~W
#' plm.fit <- plm(f,lisjobacken)
#' grob_object_list <- get_report_pages(plm.fit)
#' }
#' @export
get_report_pages.plm <- function(x,type=1,...){
    get_report_pages_fun(x,type=type)
}


#' Get gplm0 fit report pages
#'
#' Get a list of the grob objects used by \code{\link{get_report}} to create a pdf report.
#' @param x an object of class "gplm0".
#' @param type an integer to determine whether to get the pages from report type 1 or 2. Only type 1 is permissible for an object of class "gplm0" which produces a list of two grob objects, the first grob object consists of a plot collage and parameter summary table, the second a tabular rating curve.
#' @param ... further arguments passed to other methods (currently unused).
#' @seealso \code{\link{gplm0}} for fitting the gplm0 model,\code{\link{summary.gplm0}} for summaries, \code{\link{predict.gplm0}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{plot.gplm0}} to help visualize the full posterior distributions.
#' @examples
#' \dontrun{
#' data(flotemarken)
#' f <- Q~W
#' gplm0.fit <- gplm0(f,flotemarken)
#' grob_object_list <- get_report_pages(gplm0.fit)
#' }
#' @export
get_report_pages.gplm0 <- function(x,type=1,...){
    get_report_pages_fun(x,type=type)
}


#' Get gplm fit report pages
#'
#' Get a list of the grob objects used by \code{\link{get_report}} to create a pdf report.
#' @param x an object of class "gplm".
#' @param type an integer to determine whether to get the pages from report type 1 or 2. Only type 1 is permissible for an object of class "gplm" which produces a list of two grob objects, the first grob object consists of a plot collage and parameter summary table, the second a tabular rating curve.
#' @param ... further arguments passed to other methods (currently unused).
#' @seealso \code{\link{gplm}} for fitting the gplm model,\code{\link{summary.gplm}} for summaries, \code{\link{predict.gplm}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{plot.gplm}} to help visualize the full posterior distributions.
#' @examples
#' \dontrun{
#' data(bunnerviken)
#' f <- Q~W
#' gplm.fit <- gplm(f,bunnerviken)
#' grob_object_list <- get_report_pages(gplm.fit)
#' }
#' @export
get_report_pages.gplm <- function(x,type=1,...){
    get_report_pages_fun(x,type=type)
}


#' Get tournament report pages
#'
#' Get a list of the grob objects used by \code{\link{get_report}} to create a pdf report.
#' @param x an object of class "tournament".
#' @param type an integer denoting what type of report is to be produced. Defaults to type 1. Possible types are
#'                    \itemize{
#'                       \item{1}{ produces a list of two grob objects, the first grob object consists of a plot collage and parameter summary table, the second a tabular rating curve from the results of the most approriate model for a given data set, as determined by \code{\link{tournament}}.}
#'                       \item{2}{ produces a list of ten grob objects, the first four grob objects consist of plot collage and parameter summary table for the four models, the fifth object consists of model comparison plots and tables, the sixth consists of convergence diagnostics plots, and the final four grob objects consist of histogram plots of the estmated parameters for the four models.}
#'                    }
#' @param ... further arguments passed to other methods (currently unused).
#' @seealso \code{\link{tournament}} for running a tournament,\code{\link{summary.tournament}} for summaries.
#' @examples
#' \dontrun{
#' data(bunnerviken)
#' f <- Q~W
#' t_obj <- tournament(f,bunnerviken)
#' grob_object_list <- get_report_pages(t_obj)
#' }
#' @export
get_report_pages.tournament <- function(x,type=1,...){
    get_report_pages_fun(x,type=type)
}


#' Get report
#'
#' Save a pdf report of an object
#' @param x an object of class  "tournament", "plm0", "plm", "gplm0" or "gplm".
#' @param path a valid path to store the pdf report. Defaults to save report at working directory as report.pdf.
#' @param type an integer denoting what type of report is to be produced. Defaults to type 1. Only type 1 is permissible for an object of class "plm0", "plm", "gplm0" or "gplm". Possible types are
#'                    \itemize{
#'                       \item{1}{ produces a two page report consisting of a plot collage, a parameter summary table and a tabular rating curve, from the results of a model object, or if input class is "tournament", to the most approriate model, as determined by \code{\link{tournament}}.}
#'                       \item{2}{ produces a 10 page report from the results of all four models; \code{\link{plm0}}, \code{\link{plm}}, \code{\link{gplm0}} and \code{\link{gplm}}. The first four pages show a plot collage and the parameter summary table for the four models. The fifth page is a model comparison page, the sixth is a covergence diagnostics page and the last four pages show the histograms of the estimated parameters in the four models.}
#'                    }
#' @param ... further arguments passed to other methods (currently unused).
#' @seealso \code{\link{tournament}} for running a tournament,\code{\link{summary.tournament}} for summaries.
#' @examples
#' \dontrun{
#' data(bunnerviken)
#' f <- Q~W
#' t_obj <- tournament(f,bunnerviken)
#' get_report(t_obj)
#' }
#' @export
get_report <- function(x,path=NULL,type=1,...) UseMethod("get_report")


#' Get plm0 fit report
#'
#' Save a pdf report of an object
#' @param x an object of class "plm0".
#' @param path a valid path to store the pdf report. Defaults to save report at working directory as report.pdf.
#' @param type an integer denoting what type of report is to be produced. Only type 1 is permissible for an object of class "plm0" which produces a two page report consisting of a plot collage, a parameter summary table and a tabular rating curve.
#' @param ... further arguments passed to other methods (currently unused).
#' @seealso \code{\link{plm0}} for fitting the plm0 model,\code{\link{summary.plm0}} for summaries, \code{\link{predict.plm0}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{plot.plm0}} to help visualize the full posterior distributions.
#' @examples
#' \dontrun{
#' data(halla)
#' f <- Q~W
#' plm0.fit <- plm0(f,halla)
#' get_report(plm0.fit)
#' }
#' @export
get_report.plm0 <- function(x,path=NULL,type=1,...){
    report_pages <- get_report_pages_fun(x,type=type)
    save_report(report_pages,path=path,...)
}


#' Get plm fit report as pdf file
#'
#' Save a pdf report of an object
#' @param x an object of class "plm".
#' @param path a valid path to store the pdf report. Defaults to save report at working directory as report.pdf.
#' @param type an integer denoting what type of report is to be produced. Only type 1 is permissible for an object of class "plm" which produces a two page report consisting of a plot collage, a parameter summary table and a tabular rating curve.
#' @param ... further arguments passed to other methods (currently unused).
#' @seealso \code{\link{plm}} for fitting the plm model,\code{\link{summary.plm}} for summaries, \code{\link{predict.plm}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{plot.plm}} to help visualize the full posterior distributions.
#' @examples
#' \dontrun{
#' data(lisjobacken)
#' f <- Q~W
#' plm.fit <- plm(f,lisjobacken)
#' get_report(plm.fit)
#' }
#' @export
get_report.plm <- function(x,path=NULL,type=1,...){
    report_pages <- get_report_pages_fun(x,type=type)
    save_report(report_pages,path=path,...)
}


#' Get gplm0 fit report as pdf file
#'
#' Save a pdf report of an object
#' @param x an object of class "gplm0".
#' @param path a valid path to store the pdf report. Defaults to save report at working directory as report.pdf.
#' @param type an integer denoting what type of report is to be produced. Only type 1 is permissible for an object of class "gplm0" which produces a two page report consisting of a plot collage, a parameter summary table and a tabular rating curve.
#' @param ... further arguments passed to other methods (currently unused).
#' @seealso \code{\link{gplm0}} for fitting the gplm0 model,\code{\link{summary.gplm0}} for summaries, \code{\link{predict.gplm0}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{plot.gplm0}} to help visualize the full posterior distributions.
#' @examples
#' \dontrun{
#' data(flotemarken)
#' f <- Q~W
#' gplm0.fit <- gplm0(f,flotemarken)
#' get_report(gplm0.fit)
#' }
#' @export
get_report.gplm0 <- function(x,path=NULL,type=1,...){
    report_pages <- get_report_pages_fun(x,type=type)
    save_report(report_pages,path=path,...)
}


#' Get gplm fit report as pdf file
#'
#' Save a pdf report of an object
#' @param x an object of class "gplm".
#' @param path a valid path to store the pdf report. Defaults to save report at working directory as report.pdf.
#' @param type an integer denoting what type of report is to be produced. Only type 1 is permissible for an object of class "gplm" which produces a two page report consisting of a plot collage, a parameter summary table and a tabular rating curve.
#' @param ... further arguments passed to other methods (currently unused).
#' @seealso \code{\link{gplm}} for fitting the gplm model,\code{\link{summary.gplm}} for summaries, \code{\link{predict.gplm}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{plot.gplm}} to help visualize the full posterior distributions.
#' @examples
#' \dontrun{
#' data(bunnerviken)
#' f <- Q~W
#' gplm.fit <- gplm(f,bunnerviken)
#' get_report(gplm.fit)
#' }
#' @export
get_report.gplm <- function(x,path=NULL,type=1,...){
    report_pages <- get_report_pages_fun(x,type=type)
    save_report(report_pages,path=path,...)
}


#' Get tournament report
#'
#' Save a pdf report of a tournament object
#' @param x an object of class "tournament".
#' @param path a valid path to store the pdf report. Defaults to save report at working directory as report.pdf.
#' @param type an integer denoting what type of report is to be produced. Defaults to type 1. Possible types are
#'                    \itemize{
#'                       \item{1}{ produces a two page report consisting of a plot collage, a parameter summary table and a tabular rating curve, from the results of the most approriate model for a given data set, as determined by \code{\link{tournament}}.}
#'                       \item{2}{ produces a 10 page report from the results of all four models; \code{\link{plm0}}, \code{\link{plm}}, \code{\link{gplm0}} and \code{\link{gplm}}. The first four pages show a plot collage and the parameter summary table for the four models. The fifth page is a model comparison page, the sixth is a covergence diagnostics page and the last four pages show the histograms of the estimated parameters in the four models.}
#'                    }
#' @param ... further arguments passed to other methods (currently unused).
#' @seealso \code{\link{tournament}} for running a tournament,\code{\link{summary.tournament}} for summaries.
#' @examples
#' \dontrun{
#' data(bunnerviken)
#' f <- Q~W
#' t_obj <- tournament(f,bunnerviken)
#' get_report(t_obj)
#' }
#' @export
get_report.tournament <- function(x,path=NULL,type=1,...){
    report_pages <- get_report_pages_fun(x,type=type)
    save_report(report_pages,path=path,...)
}
