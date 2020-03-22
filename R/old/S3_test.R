#S3 object test
#collect a list of results from gbplm
results_obj=list('theta'=1:9,'fit'=1:100,'lower'=2:101,'upper'=0:99,beta=rnorm(100))
#set class of the list to be gbplm
attr(results_obj, "class") <- "gbplm"

#implement generic functions for the class
print.gbplm <- function(x,...){
    return('þetta er svolítið sniðugt')
}

summary.gbplm <- function(x,...){
    return("Hér verður samantekt af niðurstöðum")
}

plot.gbplm <- function(x,...){
    plot(x[['beta']])
}
