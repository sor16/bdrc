#' Compare two models using Bayes factor
#'
#' evaluate_game uses the Bayes factor of two models to determine whether one model favors the other
#'
#' @param m a list of two model objects fit on the same dataset. Model objects allowed are "gplm", "gplm0", "plm" and "plm0"
#' @return
#' A data.frame with the summary of the results of the game
#' @references B. Hrafnkelsson, H. Sigurdarson, S.M. Gardarsson, 2020, Generalization of the power-law rating curve using hydrodynamic theory and Bayesian hierarchical modeling. arXiv preprint 2010.04769.
#'
#' @seealso \code{\link{tournament}}
evaluate_game <- function(m){
    B_vec <- sapply(m,function(x) 1/mean(exp(0.5*x$Deviance_posterior)))
    BF <- B_vec[1]/B_vec[2]
    PR_m1 <- 1/(1+(1/BF))
    DIC_vec <- sapply(m,function(x) x$DIC)
    num_eff_param <- sapply(m,function(x) x$num_effective_param)
    winner <- ifelse(PR_m1>=0.75,1,2)
    data.frame(model=sapply(m,class),
               B=B_vec,
               DIC=DIC_vec,
               num_eff_param=num_eff_param,
               P=c(PR_m1,1-PR_m1),
               winner=1:2==winner)
}


#' Determine the most adequate rating curve model
#'
#' tournament compares four rating curve models of different complexities and determines which model is the most adequate
#'
#' @param ... if data and formula are set to NULL, one can add four model objects of types "gplm", "gplm0", "plm" and "plm0". This prevents the function from running all four models explicitly.
#' @param formula an object of class "formula", with discharge column name as response and stage column name as a covariate.
#' @param data data.frame containing the variables specified in formula.
#' @details  TODO
#' @return
#' An object of type "tournament" with the following elements
#' \describe{
#'  \item{\code{summary}}{a data frame with information on reults of the different games in the tournament.}
#'  \item{\code{contestants}}{model objects of types "plm0","plm","bpglm0" and "gplm" being compared.}
#'  \item{\code{winner}}{model object of the tournament winner.}
#' }
#'
#' @references B. Hrafnkelsson, H. Sigurdarson, S.M. Gardarsson, 2020, Generalization of the power-law rating curve using hydrodynamic theory and Bayesian hierarchical modeling. arXiv preprint 2010.04769.
#'
#' @seealso \code{\link{summary.tournament}}
#' @examples
#' \dontrun{
#' data(V316_river)
#' f <- Q~W
#' t_obj <- tournament(f,V316_river)
#' plot(t_obj)
#' }
#' @export
tournament <- function(formula=NULL,data=NULL,...) {
    args <- list(...)
    error_msg <- 'Please provide either formula and data (name arguments explicitly) or four model objects of types gplm, gplm0, plm and plm0.'
    if(!inherits(formula,'formula') | !is.data.frame(data)){
        args <- c(list(formula,data),args)
        if(length(args)!=4){
            stop(error_msg)
        }else{
            args_class <- unlist(lapply(args,class))
            if(!all(sort(args_class)==c('gplm','gplm0','plm','plm0'))){
                stop(error_msg)
            }else{
                names(args) <- args_class
                args <- args[order(names(args))]
                if(length(unique(lapply(args,function(x) x$data)))!=1){
                    stop('The four models added have to be fit on the same data set')
                }
                if(length(unique(lapply(args,function(x) x$c_param)))!=1){
                    stop('The four models added have to be fit either all with the same stage of zero discharge (c), or all with unknown c')
                }
            }
        }
    }else{
        args <- list()
        args$gplm <- gplm(formula, data, ...)
        args$gplm0 <- gplm0(formula, data, ...)
        args$plm <- plm(formula, data, ...)
        args$plm0 <- plm0(formula, data, ...)
    }
    round1 <- list(list(args$gplm,args$gplm0),list(args$plm,args$plm0))
    round1_res <- lapply(1:length(round1),function(i){
                    game_df <- evaluate_game(round1[[i]])
                    round_df <- data.frame(round=1,game=i)
                    cbind(round_df,game_df)
                  })
    round1_res <- do.call(rbind,round1_res)
    round1_winners <- round1_res$model[round1_res$winner]

    round2 <- lapply(1:length(round1),function(i){
        round1[[i]][[which(round1_res$winner[round1_res$game==i])]]
    })
    round2_res <- cbind(data.frame(round=2,game=3),evaluate_game(round2))
    round2_winner <- round2_res$model[round2_res$winner]
    out_obj <- list()
    attr(out_obj, "class") <- "tournament"
    out_obj$contestants <- args[!(names(args) %in% c('formula','data'))]
    out_obj$winner <- round2[[which(round2_res$winner)]]
    out_obj$summary <- rbind(round1_res,round2_res)
    return(out_obj)
}


