
print.bplm0 <- function(x,...){
    cat("\nCall:\n",
        paste(deparse(x$formula), sep = "\n", collapse = "\n"), "\n\n", sep = "")
}


summary.bplm0 <- function(x,...){
    names(x$param_summary)=paste0(names(x$param_summary),c('-2.5%','-50%','-97.5%'))
    cat("\nFormula: \n",
        paste(deparse(x$formula), sep = "\n", collapse = "\n"))
    cat("\nLatent parameters:\n")
    print(x$param_summary[1:2,],row.names = T,digits=3,right=F)
    cat("\nHyperparameters:\n")
    print(x$param_summary[3:nrow(x$param_summary),],row.names = T,digits=3,right=F)
    cat("\nDIC:",x$DIC[2])
}

plot.bplm0 <- function(x,type='rating_curve',xlab='Q',ylab='W',...){
    if(type=='rating_curve'){
        plot(x$rating_curve$median,x$rating_curve$W,type='l')
        points(x$data[,1,drop=T],x$data[,2,drop=T],...)
    }else if(type=='beta'){
        plot(sort(x$beta$W),type='l',...)
    }
}


# print.bplm <- function(x,...){
#     cat("\nCall:\n",
#         paste(deparse(x$formula), sep = "\n", collapse = "\n"), "\n\n", sep = "")
# }
#
#
# summary.bplm <- function(x,...){
#     names(x$param_summary)=paste0(names(x$param_summary),c('-2.5%','-50%','-97.5%'))
#     cat("\nFormula: \n",
#         paste(deparse(x$formula), sep = "\n", collapse = "\n"))
#     cat("\nLatent parameters:\n")
#     print(x$param_summary[1:2,],row.names = T,digits=3,right=F)
#     cat("\nHyperparameters:\n")
#     print(x$param_summary[3:nrow(x$param_summary),],row.names = T,digits=3,right=F)
#     cat("\nDIC:",x$DIC[2])
# }
#
# plot.bplm <- function(x,type='rating_curve',xlab='Q',ylab='W',...){
#     if(type=='rating_curve'){
#         plot(x$rating_curve$median,x$rating_curve$W,type='l')
#         points(x$data[,1,drop=T],x$data[,2,drop=T],...)
#     }else if(type=='beta'){
#         plot(sort(x$beta$W),type='l',...)
#     }
# }


print.bgplm <- function(x,...){
    cat("\nCall:\n",
        paste(deparse(x$formula), sep = "\n", collapse = "\n"), "\n\n", sep = "")
}


summary.bgplm <- function(x,...){
    names(x$param_summary)=paste0(names(x$param_summary),c('-2.5%','-50%','-97.5%'))
    cat("\nFormula: \n",
        paste(deparse(x$formula), sep = "\n", collapse = "\n"))
    cat("\nLatent parameters:\n")
    print(x$param_summary[1:2,],row.names = T,digits=3,right=F)
    cat("\nHyperparameters:\n")
    print(x$param_summary[3:nrow(x$param_summary),],row.names = T,digits=3,right=F)
    cat("\nDIC:",x$DIC[2])
}

plot.bgplm <- function(x,type='rating_curve',xlab='Q',ylab='W',...){
    if(type=='rating_curve'){
        plot(x$rating_curve$median,x$rating_curve$W,type='l')
        points(x$data[,1,drop=T],x$data[,2,drop=T],...)
    }else if(type=='beta'){
        plot(sort(x$beta$W),type='l',...)
    }
}


