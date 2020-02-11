#S3 object test
results_obj=list('theta'=1:9,'fit'=1:100,'lower'=2:101,'upper'=0:99,beta=rnorm(100))
attr(results_obj, "class") <- "gbplm"
print.gbplm <- function(x,...){
    return('þetta er svolítið sniðugt')
}

summary.gbplm <- function(x,...){
    return("Hér verður samantekt af niðurstöðum")
}

plot.gbplm <- function(x,...){
    plot(x[['beta']])
}
