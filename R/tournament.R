#' Compare two models using Bayes factor
#'
#' evaluate_game uses the Bayes factor of two models to determine whether one model favors the other
#'
#' @param m a list of two model objects fit on the same dataset. Model objects allowed are "gplm", "gplm0", "plm" and "plm0"
#' @param winning_criteria a numerical value between 0 and 1 which sets the threshold for which the probability of the first model in the list, given the data and calculated using Bayes factor, must exceed for it to be declared the better model of the two. This value defaults to 0.75.
#' @return
#' A data.frame with the summary of the results of the game
#' @references B. Hrafnkelsson, H. Sigurdarson, S.M. Gardarsson, 2020, Generalization of the power-law rating curve using hydrodynamic theory and Bayesian hierarchical modeling. arXiv preprint 2010.04769.
#'
#' @seealso \code{\link{tournament}}
#' @keywords internal
evaluate_game <- function(m,winning_criteria=0.75){
    B_vec <- sapply(m,function(x) 1/mean(exp(0.5*x$Deviance_posterior)))
    BF <- B_vec[1]/B_vec[2]
    PR_m1 <- 1/(1+(1/BF))
    DIC_vec <- sapply(m,function(x) x$DIC)
    num_eff_param <- sapply(m,function(x) x$num_effective_param)
    winner <- if(PR_m1>=winning_criteria) 1 else 2
    data.frame(model=sapply(m,class),
               B=B_vec,
               DIC=DIC_vec,
               num_eff_param=num_eff_param,
               P=c(PR_m1,1-PR_m1),
               winner=1:2==winner)
}


#' Tournament - Model comparison
#'
#' tournament compares four rating curve models of different complexities and determines the model that provides the best fit of the data at hand.
#'
#' @param formula an object of class "formula", with discharge column name as response and stage column name as a covariate.
#' @param data data.frame containing the variables specified in formula.
#' @param ... optional arguments passed to the model functions. Also, if data and formula are NULL, one can add four model objects of types "gplm", "gplm0", "plm" and "plm0". This runs the tournament for the input models and prevents running all four models explicitly.
#' @param winning_criteria a numerical value between 0 and 1 which sets the threshold for which the probability of the more complex model given the data in each model comparison, must exceed for it to be declared the better model. This value defaults to 0.75 to favor the less complex models when the superiority of the more complex model is somewhat ambiguous. See the Details section.
#' @details Tournament is a comparison method that uses Bayes factor to compute the posterior probabilities of the models and select the most appropriate of the four models given the data. The first round of model comparisons sets up model types "gplm" vs. "gplm0" and "plm" vs. "plm0". If the posterior probability of the more complex model ("gplm" and "plm", respectively) exceeds the "winning_criteria" (default value = 0.75) then it is chosen as the more appropriate model and moves on to the second and final round, where the winners from the first round will be compared in the same way. In the second round, if the more complex model (now the generalized power-law model) exceeds the same "winning_criteria" then it is chosen as the overall tournament winner and deemed the most appropriate model given the data. In each of the three matches, the posterior probabilities of the models are computed using the Bayes factor and assuming a priori that the two models were equally likely [see Jeffreys (1961) and Kass and Raftery (1995)].
#' @return
#' An object of type "tournament" with the following elements
#' \describe{
#'  \item{\code{summary}}{a data frame with information on results of the different games in the tournament.}
#'  \item{\code{contestants}}{model objects of types "plm0","plm","gplm0" and "gplm" being compared.}
#'  \item{\code{winner}}{model object of the tournament winner.}
#' }
#'
#' @references B. Hrafnkelsson, H. Sigurdarson, S.M. Gardarsson, 2020, Generalization of the power-law rating curve using hydrodynamic theory and Bayesian hierarchical modeling. arXiv preprint 2010.04769.
#' @references Jeffreys, H. (1961). Theory of Probability, Third Edition. Oxford University Press.
#' @references Kass, R., and A. Raftery, A. (1995). Bayes Factors. Journal of the American Statistical Association, 90, 773-795.
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
tournament <- function(formula=NULL,data=NULL,...,winning_criteria=0.75) {
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
                    game_df <- evaluate_game(round1[[i]],winning_criteria)
                    round_df <- data.frame(round=1,game=i)
                    cbind(round_df,game_df)
                  })
    round1_res <- do.call(rbind,round1_res)
    round1_winners <- round1_res$model[round1_res$winner]

    round2 <- lapply(1:length(round1),function(i){
        round1[[i]][[which(round1_res$winner[round1_res$game==i])]]
    })
    round2_res <- cbind(data.frame(round=2,game=3),evaluate_game(round2,winning_criteria))
    round2_winner <- round2_res$model[round2_res$winner]
    out_obj <- list()
    attr(out_obj, "class") <- "tournament"
    out_obj$contestants <- args[!(names(args) %in% c('formula','data'))]
    out_obj$winner <- round2[[which(round2_res$winner)]]
    out_obj$summary <- rbind(round1_res,round2_res)
    return(out_obj)
}


