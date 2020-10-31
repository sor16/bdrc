print_fun <- function(x){
    cat("\nCall:\n",
        paste(deparse(x$formula), sep = "\n", collapse = "\n"), "\n\n", sep = "")
}

summary_fun <- function(x){
    names(x$param_summary)=paste0(names(x$param_summary),c('-2.5%','-50%','-97.5%'))
    cat("\nFormula: \n",
        paste(deparse(x$formula), sep = "\n", collapse = "\n"))
    cat("\nLatent parameters:\n")
    print(x$param_summary[1:2,],row.names = T,digits=3,right=F)
    cat("\nHyperparameters:\n")
    print(x$param_summary[3:nrow(x$param_summary),],row.names = T,digits=3,right=F)
    cat("\nDIC:",x$DIC_summary[,2])
}


plot_fun <- function(x,type){

}

predict_fun <- function(x,type){

}


print.bplm0 <- function(x,...){
    print_fun(x)
}

summary.bplm0 <- function(x,...){
    summary_fun(x)
}

plot.bplm0 <- function(x,type='rating_curve',...){
}

predict.bplm0 <- function(x,...){
}

print.bplm <- function(x,...){
    print_fun(x)
}

summary.bplm <- function(x,...){
    summary_fun(x)
}

plot.bplm <- function(x,type='rating_curve',...){
}

predict.bplm <- function(x,...){
}

print.bgplm0 <- function(x,...){
    print_fun(x)
}

summary.bgplm0 <- function(x,...){
    summary_fun(x)
}

plot.bgplm0 <- function(x,type='rating_curve',...){
}

predict.bgplm0 <- function(x,...){
}

print.bgplm <- function(x,...){
    print_fun(x)
}

summary.bgplm <- function(x,...){
    summary_fun()
}

plot.bgplm <- function(x,type='rating_curve',...){
}

predict.bgplm <- function(x,...){
}
