#' Compare two models using a specified model selection criteria
#'
#' evaluate_game uses the Bayes factor, DIC or WAIC of two models to determine whether one model is more appropriate than the other
#'
#' @param m a list of two model objects fit on the same dataset. The allowed model objects are "gplm", "gplm0", "plm" and "plm0"
#' @param method a string specifying the method used to declare a winner. The allowed methods are "Delta_WAIC", "Delta_DIC" and "Bayes_factor".
#' @param winning_criteria a numerical value which sets the threshold for which the first model in the list must exceed for it to be declared the more appropriate model. This value defaults to 0.75 when Bayes factor is chosen as the model selection criteria but defaults to 1.5 when either DIC or WAIC are chosen.
#' @return
#' A data.frame with the summary of the results of the game
#' @references Hrafnkelsson, B., Sigurdarson, H., and Gardarsson, S. M. (2022). Generalization of the power-law rating curve using hydrodynamic theory and Bayesian hierarchical modeling, Environmetrics, 33(2):e2711.
#'
#' @seealso \code{\link{tournament}}
#' @keywords internal
evaluate_game <- function(m,method,winning_criteria){
    B_vec <- sapply(m,function(x) 1/mean(exp(0.5*x$Deviance_posterior)))
    BF <- B_vec[1]/B_vec[2]
    PR_m1 <- 1/(1+(1/BF))
    DIC_vec <- sapply(m,function(x) x$DIC)
    WAIC_vec <- sapply(m,function(x) x$WAIC)
    eff_num_param_DIC <- sapply(m,function(x) x$effective_num_param_DIC)
    eff_num_param_WAIC <- sapply(m,function(x) x$effective_num_param_WAIC)
    DDIC <- DIC_vec[2]-DIC_vec[1]
    DWAIC <- WAIC_vec[2]-WAIC_vec[1]
    delta_vec <- c('Bayes_factor'=PR_m1,'Delta_DIC'=DDIC,'Delta_WAIC'=DWAIC)
    winner <- if(delta_vec[[method]]>=winning_criteria) 1 else 2
    data.frame(model=sapply(m,class),
               BF=B_vec,
               P=c(PR_m1,1-PR_m1),
               eff_num_param_DIC=eff_num_param_DIC,
               DIC=DIC_vec,
               Delta_DIC=c(DDIC,NA),
               eff_num_param_WAIC=eff_num_param_WAIC,
               WAIC=WAIC_vec,
               Delta_WAIC=c(DWAIC,NA),
               winner=1:2==winner)
}


#' Tournament - Model comparison
#'
#' tournament compares four rating curve models of different complexities and determines the model that provides the best fit of the data at hand.
#'
#' @param formula an object of class "formula", with discharge column name as response and stage column name as a covariate.
#' @param data data.frame containing the variables specified in formula.
#' @param ... optional arguments passed to the model functions. Also, if data and formula are NULL, one can add four model objects of types "gplm", "gplm0", "plm" and "plm0". This runs the tournament for the input models and prevents running all four models explicitly.
#' @param method a string specifying the method used to declare the winner. The allowed methods are "Delta_WAIC", "Delta_DIC" and "Bayes_factor".
#' @param winning_criteria a numerical value which sets a threshold which the more complex model in each model comparison must exceed to be deemed the more appropriate model. When "Bayes_factor" is used as the model selection criterion this value should be between 0 and 1 as this sets the threshold for which the probability of the more complex model given the data in each model comparison, must exceed for it to be declared the more appropriate model. This value defaults to 0.75 to favor the less complex models when the superiority of the more complex model is somewhat ambiguous. When "Delta_WAIC" or "Delta_DIC" is used as the model selection criterion this value should be a real number. In this case this value defaults to 1.5 to favor the less complex models when the superiority of the more complex model is somewhat ambiguous. See the Details section.
#' @details Tournament is a comparison method that uses Bayes factor to compute the posterior probabilities of the models, or the differences in either DIC or WAIC, and select the most appropriate of the four models given the data. The first round of model comparisons sets up model types "gplm" vs. "gplm0" and "plm" vs. "plm0". If the more complex model ("gplm" and "plm", respectively) exceeds the "winning_criteria" (default value = 1.5, when using "Delta_WAIC") then it is chosen as the more appropriate model and moves on to the second and final round, where the winners from the first round will be compared in the same way. In the second round, if the more complex model (now the generalized power-law model) exceeds the same "winning_criteria" then it is chosen as the overall tournament winner and deemed the most appropriate model given the data.
#' @return
#' An object of type "tournament" with the following elements
#' \describe{
#'  \item{\code{contestants}}{model objects of types "plm0","plm","gplm0" and "gplm" being compared.}
#'  \item{\code{winner}}{model object of the tournament winner.}
#'  \item{\code{summary}}{a data frame with information on results of the different games in the tournament.}
#'  \item{\code{info}}{information about the tournament; the overall winner of the tournament; the method used to compare the models; the numerical threshold value for which a more complex model in a model comparison must have exceeded for it to have been deemed the appropriate model.}
#' }
#'
#' @references B. Hrafnkelsson, H. Sigurdarson, S.M. Gardarsson. (2020). Generalization of the power-law rating curve using hydrodynamic theory and Bayesian hierarchical modeling. arXiv preprint 2010.04769.
#' @references Jeffreys, H. (1961). Theory of Probability, Third Edition. Oxford University Press.
#' @references Kass, R., and A. Raftery, A. (1995). Bayes Factors. Journal of the American Statistical Association, 90, 773-795.
#' @references Spiegelhalter, D., Best, N., Carlin, B., Van Der Linde, A. (2002). Bayesian measures of model complexity and fit. Journal of the Royal Statistical Society: Series B (Statistical Methodology) 64(4), 583–639.
#' @references Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and widely applicable information criterion in singular learning theory. J. Mach. Learn. Res. 11, 3571–3594.
#'
#' @seealso \code{\link{summary.tournament}} and \code{\link{plot.tournament}}
#' @examples
#' \donttest{
#' data(krokfors)
#' set.seed(1)
#' t_obj <- tournament(formula=Q~W,data=krokfors,num_cores=2)
#' t_obj
#' summary(t_obj)
#' }
#' @export
tournament <- function(formula=NULL,data=NULL,...,method='Delta_WAIC',winning_criteria=NULL) {
    args <- list(...)
    default_win_crit <- c('Delta_WAIC'=1.5,'Delta_DIC'=1.5,'Bayes_factor'=0.75)
    error_msg <- "The method input must contain a string indicating the method to be used for comparing the models. The methods are 'Delta_WAIC' (default), 'Delta_DIC' and 'Bayes_factor'."
    if( is.null(method) ){
        stop(error_msg)
    }else{
        if( !(method%in%c('Delta_WAIC','Delta_DIC','Bayes_factor')) ){
            stop(error_msg)
        }
    }
    if(grepl("Delta",method)){
        error_msg <- "The winning_criteria when the method is set to 'Delta_DIC' or 'Delta_WAIC' must be a numerical real value. It is a threshold which a more complex model needs to surpass to be declared the more appropriate model when compared with a less complex model. This value defaults to 1.5 when the method is set to 'Delta_DIC' or 'Delta_WAIC'."
    }else{
        error_msg <- "The winning_criteria when the method is set to 'Bayes_factor' must be a numerical value between 0 and 1. It is a threshold which a more complex model needs to surpass to be declared the more appropriate model when compared with a less complex model. This value defaults to 0.75 when the method is set to 'Bayes_factor'."
    }
    if(!is.null(winning_criteria)){
        if(class(winning_criteria)!="numeric"){
            stop(error_msg)
        }else if( method=='Bayes_factor' & abs(winning_criteria-0.5)>0.5 ){
            stop(error_msg)
        }
    }
    if(is.null(winning_criteria)){
        winning_criteria <- default_win_crit[[method]]
    }
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
        message('Running tournament:')
        args$gplm <- gplm(formula, data, ...)
        message('25% - gplm finished')
        args$gplm0 <- gplm0(formula, data, ...)
        message('50% - gplm0 finished')
        args$plm <- plm(formula, data, ...)
        message('75% - plm finished')
        args$plm0 <- plm0(formula, data, ...)
        message('100% - plm0 finished')
    }
    round1 <- list(list(args$gplm,args$gplm0),list(args$plm,args$plm0))
    round1_res <- lapply(1:length(round1),function(i){
        game_df <- evaluate_game(round1[[i]],method,winning_criteria)
        round_df <- data.frame(round=1,game=i)
        cbind(round_df,game_df)
    })
    round1_res <- do.call(rbind,round1_res)
    round1_winners <- round1_res$model[round1_res$winner]

    round2 <- lapply(1:length(round1),function(i){
        round1[[i]][[which(round1_res$winner[round1_res$game==i])]]
    })
    round2_res <- cbind(data.frame(round=2,game=3),evaluate_game(round2,method,winning_criteria))
    round2_winner <- round2_res$model[round2_res$winner]
    out_obj <- list()
    attr(out_obj, "class") <- "tournament"
    out_obj$contestants <- args[!(names(args) %in% c('formula','data'))]
    out_obj$winner <- round2[[which(round2_res$winner)]]
    out_obj$summary <- rbind(round1_res,round2_res)
    out_obj$info <- list("winner"=class(round2[[which(round2_res$winner)]]),"method"=method,"winning_criteria"=winning_criteria)
    return(out_obj)
}


