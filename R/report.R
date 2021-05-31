### TODO: Option in the predict function?? I.e. an help function for predict method. Ath aukastafi. Titill á töflu
#' @importFrom gridExtra tableGrob ttheme_minimal
predict_matrix <- function(x){
    # c_param <- if(is.null(x$run_info$c_param)) median(x$c_posterior) else x$run_info$c_param
    # c_param <- ceiling(c_param*100)/100
    min_rc_h <- min(x$rating_curve[,'h'])
    min_rc_h <- ceiling(min_rc_h*100)/100
    grid_max <- x$run_info$h_max
    p_dat <- predict(x,newdata=seq(min_rc_h,grid_max,by=0.01))[c('h','median')]  # c_param (breyta)
    p_dat$decimal <- floor(p_dat$h*10)/10
    first_decimal <- length(p_dat$decimal[p_dat$decimal==p_dat$decimal[1]])
    if(first_decimal!=10) {
        n <- 10-first_decimal
        top_rows <- data.frame(h=sapply(n:1,function(x) p_dat$h[1]-0.01*x),median=rep(0,n),decimal=rep(p_dat$decimal[1],n))
        p_dat <- rbind(top_rows,p_dat)
    }
    last_decimal <- length(p_dat$decimal[p_dat$decimal==p_dat$decimal[length(p_dat$decimal)]])
    if(last_decimal!=10){
        m <- 10-last_decimal
        bot_rows <- data.frame(h=sapply(1:m,function(x) p_dat$h[length(p_dat$h)]+0.01*x),median=rep(NA,m),decimal=rep(p_dat$decimal[length(p_dat$decimal)],m))
        p_dat <- rbind(p_dat,bot_rows)
    }
    p_mat <- lapply(unique(p_dat$decimal),function(d) p_dat$median[p_dat$decimal==d])
    p_mat <- do.call('rbind',p_mat)
    rownames(p_mat) <- unique(p_dat$decimal)
    colnames(p_mat) <- seq(0,0.09,by=0.01)
    p_mat <- round(p_mat,digits=3)
    p_mat <- tableGrob(p_mat,theme=ttheme_minimal(base_family = "Times",
                                                  core=list(bg_params = list(fill = blues9[1:2], col=NA),fg_params=list(fontface=3)),
                                                  colhead=list(fg_params=list(col="black",fontface=2L)),
                                                  rowhead=list(fg_params=list(col="black",fontface=2L))))
    return(p_mat)
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
                    data.frame(rating_curve_x_min=quantile(df$rating_curve,0.025),rating_curve_x_max=1.01*max(quantile(df[df$h==max(df$h),]$rating_curve,0.975),max(q_dat)),
                               rating_curve_y_min=min(df$h),rating_curve_y_max=1.01*max(df$h)-0.01*min(df$h),
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
                                            table <- rbind(m$param_summary,m$Deviance_summary)
                                            row.names(table) <- c(sapply(1:length(param), function(x) get_param_expression(x)[[1]]),"Deviance")
                                            table <- format(round(table,digits=3),nsmall=3)
                                            tableGrob(table,theme=ttheme_minimal(base_family="Times",rowhead=list(fg_params=list(parse=TRUE))))
                                        })
        output_list$p_mat <- predict_matrix(m_obj[[1]])
        output_list$obj_class <- class(m_obj[[1]])
    }else{
        output_list$main_page_table <- lapply(m_obj,function(m){
                                            param <- get_param_names(class(m),m$run_info$c_param)
                                            table <- rbind(m$param_summary,m$Deviance_summary)
                                            row.names(table) <- c(sapply(1:length(param), function(x) get_param_expression(x)[[1]]),"Deviance")
                                            rhat_col <- c(m$r_hat$r_hat,NA)
                                            neff_col <- c(m$num_effective_samples$num_effective_samples,NA)
                                            table <- cbind(table,'R_hat'=rhat_col,'n_eff'=neff_col)
                                            table[,c('lower','median','upper','R_hat')] <- format(round(table[,c('lower','median','upper','R_hat')],digits=3),nsmall=3)
                                            table['Deviance',c('R_hat','n_eff')] <- ''
                                            tableGrob(table,theme=ttheme_minimal(base_family="Times",rowhead=list(fg_params=list(parse=TRUE))))
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
        output_list$mcmc_table <- tableGrob(mcmc_table,theme=ttheme_minimal(base_family="Times"))
        tour_table <- t_obj$summary[c("round","game","model","DIC","P","winner")]
        tour_table[c('DIC','P')] <- round(tour_table[c('DIC','P')],digits=2)
        output_list$tour_table <- tableGrob(tour_table,theme=ttheme_minimal(base_family="Times"),rows=NULL)
        output_list$dev_boxplot <- autoplot(t_obj,type='deviance')
        output_list$conv_diag_plots <- plot(t_obj,type='convergence_diagnostics')
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
                                      top=textGrob(report_components$obj_class,gp=gpar(fontfamily="Times",fontsize=22,facetype='bold')))
        predict_mat <- arrangeGrob(report_components$p_mat,
                                   nrow=1,
                                   as.table=TRUE,
                                   heights=c(1),
                                   top=textGrob(paste0('Rating curve predictions for ',report_components$obj_class),gp=gpar(fontsize=22,facetype='bold',fontfamily="Times")))
        report_pages <- list(main_page_plot,predict_mat)
    }else{
        main_page_plots <- lapply(names(report_components$main_page_plots),function(m){
                                  arrangeGrob(report_components$main_page_plots[[m]],
                                              report_components$main_page_table[[m]],
                                              nrow=2,
                                              as.table=TRUE,
                                              heights=c(5,3),
                                              top=textGrob(m,gp=gpar(fontfamily="Times",fontsize=22,facetype='bold')))
                           })
        tournament_page <- arrangeGrob(arrangeGrob(report_components$mcmc_table,report_components$tour_table,nrow=1),
                                       arrangeGrob(report_components$dev_boxplot,ncol=2),
                                       nrow=2,
                                       as.table=TRUE,
                                       heights=c(1,1),
                                       top=textGrob('Model comparison and MCMC diagnostics',gp=gpar(fontsize=22,facetype='bold',fontfamily="Times")))
        convergence_page <- arrangeGrob(report_components$conv_diag_plots,
                                        as.table=TRUE,
                                        top=textGrob('',gp=gpar(fontsize=22,facetype='bold',fontfamily="Times")))
        histogram_pages <- lapply(names(report_components$main_page_plots), function(m) {
                                  arrangeGrob(arrangeGrob(grobs=report_components$mcmc_hist_list[[m]],nrow=4,ncol=3),
                                              top=textGrob(paste0('Estimated parameters of ',m),gp=gpar(fontsize=20,facetype='bold',fontfamily="Times")))
                            })
        report_pages <- list(main_page_plots,tournament_page,convergence_page,histogram_pages)
    }
    return(report_pages)
}

#' @importFrom grDevices pdf dev.off
save_report <- function(report_pages,path=NULL,paper='a4',width=8,height=11){
    if(is.null(path)){
        path <- 'report.pdf'
    }
    pdf(file=path,paper=paper,width=width,height=height)
    for(i in 1:length(report_pages)){
        if(any(class(report_pages[[i]])=='list')){
            for (j in 1:length(report_pages[[i]])) {
                grid.arrange(report_pages[[i]][[j]],as.table=T)
            }
        }else{
            grid.arrange(report_pages[[i]],as.table=T)
        }
    }
    dev.off()
}

#### S3 methods
get_report_pages <- function(x,type=1,...) UseMethod("get_report_pages")

#' @export
get_report_pages.plm0 <- function(x,type=1,...){
    get_report_pages_fun(x,type=type)
}

#' @export
get_report_pages.plm <- function(x,type=1,...){
    get_report_pages_fun(x,type=type)
}

#' @export
get_report_pages.gplm0 <- function(x,type=1,...){
    get_report_pages_fun(x,type=type)
}

#' @export
get_report_pages.gplm <- function(x,type=1,...){
    get_report_pages_fun(x,type=type)
}

#' @export
get_report_pages.tournament <- function(x,type=1,...){
    get_report_pages_fun(x,type=type)
}

get_report <- function(x,path=NULL,type=1,...) UseMethod("get_report")

#' @export
get_report.plm0 <- function(x,path=NULL,type=1,...){
    report_pages <- get_report_pages_fun(x,type=type)
    save_report(report_pages,path=path,...)
}

#' @export
get_report.plm <- function(x,path=NULL,type=1,...){
    report_pages <- get_report_pages_fun(x,type=type)
    save_report(report_pages,path=path,...)
}

#' @export
get_report.gplm0 <- function(x,path=NULL,type=1,...){
    report_pages <- get_report_pages_fun(x,type=type)
    save_report(report_pages,path=path,...)
}

#' @export
get_report.gplm <- function(x,path=NULL,type=1,...){
    report_pages <- get_report_pages_fun(x,type=type)
    save_report(report_pages,path=path,...)
}

#' @export
get_report.tournament <- function(x,path=NULL,type=1,...){
    report_pages <- get_report_pages_fun(x,type=type)
    save_report(report_pages,path=path,...)
}
