#' Compare two models using a specified model selection criteria
#'
#' evaluate_game uses WAIC, DIC or the posterior probabilities of the models, calculated with Bayes factor, to determine whether one model is more appropriate than the other model given the data at hand.
#'
#' @param m a list of two model objects fit on the same dataset. The allowed model objects are "gplm", "gplm0", "plm" and "plm0"
#' @param method a string specifying the method used to estimate the predictive performance of the models. The allowed methods are "WAIC", "DIC" and "Posterior_probability".
#' @param winning_criteria a numerical value which sets the threshold which the first model in the list must exceed for it to be declared the more appropriate model. This value defaults to 2.2 for methods "WAIC" and "DIC", but defaults to 0.75 for method "Posterior_probability".
#' @return
#' A data.frame with the summary of the results of the game
#' @references Hrafnkelsson, B., Sigurdarson, H., and Gardarsson, S. M. (2022). Generalization of the power-law rating curve using hydrodynamic theory and Bayesian hierarchical modeling, Environmetrics, 33(2):e2711.
#'
#' @seealso \code{\link{tournament}}
#' @keywords internal
evaluate_game <- function(m,method,winning_criteria){
    if(method=="Posterior_probability"){
        ml <- sapply(m,function(x) 1/mean(exp(0.5*x$Deviance_posterior)))
        BF <- ml[1]/ml[2]
        PR_m1 <- 1/(1+(1/BF))
        df <- data.frame("marg_lik"=ml,"Post_prob"=c(PR_m1,1-PR_m1))
    }else if(method=="DIC"){
        D_hat <- sapply(m,function(x) x$D_hat)
        eff_num_param_DIC <- sapply(m,function(x) x$effective_num_param_DIC)
        DIC_vec <- sapply(m,function(x) x$DIC)
        DDIC <- DIC_vec[2]-DIC_vec[1]
        df <- data.frame("D_hat"=D_hat,"eff_num_param"=eff_num_param_DIC,"DIC"=DIC_vec,"Delta_DIC"=c(DDIC,NA))
    }else if(method=="WAIC"){
        lppd <- sapply(m,function(x) x$lppd)
        eff_num_param_WAIC <- sapply(m,function(x) x$effective_num_param_WAIC)
        WAIC_vec <- sapply(m,function(x) x$WAIC)
        DWAIC <- WAIC_vec[2]-WAIC_vec[1]
        df <- data.frame("lppd"=lppd,"eff_num_param"=eff_num_param_WAIC,"WAIC"=WAIC_vec,"Delta_WAIC"=c(DWAIC,NA))
    }
    winner <- if(df[1,ncol(df)]>=winning_criteria) 1 else 2
    data.frame(model=sapply(m,class),df,winner=1:2==winner)
}


#' Tournament - Model comparison
#'
#' tournament compares four rating curve models of different complexities and determines the model that provides the best fit of the data at hand.
#'
#' @param formula an object of class "formula", with discharge column name as response and stage column name as a covariate.
#' @param data data.frame containing the variables specified in formula.
#' @param model_list list of exactly four model objects of types "plm0","plm","gplm0" and "gplm" to be used in the tournament. Note that all of the model objects are required to be run with the same data and same c_param.
#' @param method a string specifying the method used to estimate the predictive performance of the models. The allowed methods are "WAIC", "DIC" and "Posterior_probability".
#' @param winning_criteria a numerical value which sets a threshold the more complex model in each model comparison must exceed to be deemed the more appropriate model. See the Details section.
#' @param ... optional arguments passed to the model functions.
#' @details Tournament is a model comparison method that uses WAIC to estimate the predictive performance of the four models and select the most appropriate model given the data. The first round of model comparisons sets up two games between model types, "gplm" vs. "gplm0" and "plm" vs. "plm0". The two comparisons are conducted such that if the WAIC of the more complex model ("gplm" and "plm", respectively) is smaller than the WAIC of the simpler models ("gplm0" and "plm0", respectively) by an input argument called the \code{winning_criteria} (default value = 2.2), then it is chosen as the more appropriate model. If not, the simpler model is chosen. The more appropriate models move on to the second round and are compared in the same way. The winner of the second round is chosen as the overall tournament winner and deemed the most appropriate model given the data.
#'
#' The default method "WAIC", or the Widely Applicable Information Criterion (see Watanabe (2010)), is used to estimate the predictive performance of the models. This method is a fully Bayesian method that uses the full set of posterior draws to estimate of the expected log pointwise predictive density.
#'
#' Method "DIC", or Deviance Information Criterion (see Spiegelhalter (2002)), is similar to the "WAIC" but instead of using the full set of posterior draws to compute the estimate of the expected log pointwise predictive density, it uses a point estimate of the posterior distribution.
#'
#' Method "Posterior_probability" uses the posterior probabilities of the models, calculated with Bayes factor (see Jeffreys (1961) and Kass and Raftery (1995)), to compare the models, where all the models are assumed a priori to be equally likely. This method is not chosen as the default method because the Bayes factor calculations can be quite unstable.
#'
#' When methods "WAIC" or "DIC" are used, the \code{winning_criteria} should be a real number. The winning criteria is a threshold value which the more complex model in each model comparison must exceed for it to be declared the more appropriate model. Setting the winning criteria slightly above 0 (default value = 2.2 for both "WAIC" and "DIC") gives the less complex model in each comparison a slight advantage. When method "Posterior_probability" is used, the winning criteria should be a real value between 0 and 1 (default value = 0.75). This sets the threshold value for which the posterior probability of the more complex model, given the data, in each model comparison must exceed for it to be declared the more appropriate model. In all three cases, the default value is selected so as to give the less complex models a slight advantage, and should give more or less consistent results when applying the tournament to real world data.
#' @return
#' An object of type "tournament" with the following elements
#' \describe{
#'  \item{\code{contestants}}{model objects of types "plm0","plm","gplm0" and "gplm" being compared.}
#'  \item{\code{winner}}{model object of the tournament winner.}
#'  \item{\code{summary}}{a data frame with information on results of the different games in the tournament.}
#'  \item{\code{info}}{specifics about the tournament; the overall winner; the method used; and the winning criteria.}
#' }
#'
#' @references Hrafnkelsson, B., Sigurdarson, H., and Gardarsson, S. M. (2022). Generalization of the power-law rating curve using hydrodynamic theory and Bayesian hierarchical modeling, Environmetrics, 33(2):e2711.
#' @references Jeffreys, H. (1961). Theory of Probability, Third Edition. Oxford University Press.
#' @references Kass, R., and A. Raftery, A. (1995). Bayes Factors. Journal of the American Statistical Association, 90, 773-795.
#' @references Spiegelhalter, D., Best, N., Carlin, B., Van Der Linde, A. (2002). Bayesian measures of model complexity and fit. Journal of the Royal Statistical Society: Series B (Statistical Methodology) 64(4), 583–639.
#' @references Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and widely applicable information criterion in singular learning theory. J. Mach. Learn. Res. 11, 3571–3594.
#'
#' @seealso \code{\link{plm0}} \code{\link{plm}}, \code{\link{gplm0}},\code{\link{gplm}} \code{\link{summary.tournament}} and \code{\link{plot.tournament}}
#' @examples
#' \donttest{
#' data(krokfors)
#' set.seed(1)
#' t_obj <- tournament(formula=Q~W,data=krokfors,num_cores=2)
#' t_obj
#' summary(t_obj)
#' }
#' @export
tournament <- function(formula=NULL,data=NULL,model_list=NULL,method='WAIC',winning_criteria=NULL,...) {
    default_win_crit <- c('WAIC'=2.2,'DIC'=2.2,'Posterior_probability'=0.75)
    error_msg <- "The method input must contain a string indicating the method to be used for comparing the models. The methods are 'WAIC' (default), 'DIC' and 'Posterior_probability'."
    if(is.null(method)){
        stop(error_msg)
    }else{
        if(length(method)>1){
            stop(error_msg)
        }
        if(!(method%in%c('WAIC','DIC','Posterior_probability'))){
            stop(error_msg)
        }
    }
    if(grepl("IC",method)){
        error_msg <- "The winning_criteria when the method is set to 'WAIC' or 'DIC' must be a numerical real value. This is the threshold value which a model comparison statistic of the more complex model must surpass to be declared the more appropriate model when compared with a less complex model. This value defaults to 2.2 when the models are compared by their 'WAIC' or 'DIC' values."
    }else{
        error_msg <- "The winning_criteria when the method is set to 'Posterior_probability' must be a numerical value between 0 and 1. This is the threshold for which the posterior probability of a more complex model, which is calculated using the Bayes factor, needs to surpass to be declared the appropriate model when compared with a less complex model. This value defaults to 0.75 when the method is set to 'Posterior_probability'."
    }
    if(!is.null(winning_criteria)){
        if(length(method)>1){
            stop(error_msg)
        }
        if(!inherits(winning_criteria,"numeric")){
            stop(error_msg)
        }else if( method=='Posterior_probability' & abs(winning_criteria-0.5)>0.5 ){
            stop(error_msg)
        }
    }
    if(is.null(winning_criteria)){
        winning_criteria <- default_win_crit[method]
    }
    error_msg <- 'Please provide; formula and data (name arguments explicitly) or model_list with four model objects of types gplm, gplm0, plm and plm0.'
    models <- list()
    if(!is.null(model_list) | (is.null(model_list) & inherits(formula,'list'))){
        if(is.null(model_list) & inherits(formula,'list')){
            model_list=formula
        }
        if(length(model_list)!=4){
            stop(error_msg)
        }
        models_class <- unlist(lapply(model_list,class))
        if(!all(sort(models_class)==c('gplm','gplm0','plm','plm0'))){
            stop(error_msg)
        }
        models <- model_list
        names(models) <- models_class
        models <- models[order(names(models))]
        #TODO: make sure data argument matches models in model_list
        if(length(unique(lapply(models,function(x) x$data)))!=1){
            stop('The four models added have to be fit on the same data set')
        }
        if(length(unique(lapply(models,function(x) x$c_param)))!=1){
            stop('The four models added have to be fit either all with the same stage of zero discharge (c), or all with unknown c')
        }
    }else{
        message('Running tournament:')
        models$gplm <- gplm(formula, data, ...)
        message('25% - gplm finished')
        models$gplm0 <- gplm0(formula, data, ...)
        message('50% - gplm0 finished')
        models$plm <- plm(formula, data, ...)
        message('75% - plm finished')
        models$plm0 <- plm0(formula, data, ...)
        message('100% - plm0 finished')
    }
    round1 <- list(list(models$gplm,models$gplm0),list(models$plm,models$plm0))
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
    out_obj$contestants <- models
    out_obj$winner <- round2[[which(round2_res$winner)]]
    out_obj$summary <- rbind(round1_res,round2_res)
    out_obj$info <- list("winner"=class(round2[[which(round2_res$winner)]]),"method"=method,"winning_criteria"=winning_criteria)
    return(out_obj)
}


