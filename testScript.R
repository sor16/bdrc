library(devtools)
load_all()
library(RCmodels)
library(dplyr)
library(ggplot2)

#read data
data(V316_river)
rc_dat <- V316_river
rc_formula <- as.formula('Q~W')

## bplm0
bplm0.fit <- bplm0(formula = rc_formula,data = rc_dat)

summary(bplm0.fit)


ggplot(data=bplm0.fit$rating_curve) +
  geom_point(data=rc_dat,aes(Q,W)) +
  geom_line(aes(median,w)) +
  geom_line(aes(lower,w),linetype='dashed') +
  geom_line(aes(upper,w),linetype='dashed')


bplm0.fit_known_c <- bplm0(rc_formula,rc_dat,c_param=0.75)

summary(bplm0.fit_known_c)

ggplot(data=bplm0.fit_known_c$rating_curve) +
  geom_point(data=rc_dat,aes(Q,W)) +
  geom_line(aes(median,w)) +
  geom_line(aes(lower,w),linetype='dashed') +
  geom_line(aes(upper,w),linetype='dashed')

## bplm
bplm.fit <- bplm(rc_formula,rc_dat)

summary(bplm.fit)

ggplot(data=bplm.fit$rating_curve) +
    geom_point(data=rc_dat,aes(Q,W)) +
    geom_line(aes(median,w)) +
    geom_line(aes(lower,w),linetype='dashed') +
    geom_line(aes(upper,w),linetype='dashed')

bplm.fit_known_c <- bplm(rc_formula,rc_dat,c_param=0.75)

summary(bplm.fit_known_c)

ggplot(data=bplm.fit_known_c$rating_curve) +
    geom_point(data=rc_dat,aes(Q,W)) +
    geom_line(aes(median,w)) +
    geom_line(aes(lower,w),linetype='dashed') +
    geom_line(aes(upper,w),linetype='dashed')


#bgplm0
bgplm0.fit <- bgplm0(rc_formula,rc_dat)

summary(bgplm0.fit)

ggplot(data=bgplm0.fit$rating_curve) +
    geom_point(data=rc_dat,aes(Q,W)) +
    geom_line(aes(median,w)) +
    geom_line(aes(lower,w),linetype='dashed') +
    geom_line(aes(upper,w),linetype='dashed')

bgplm0.fit_known_c <- bgplm0(rc_formula,rc_dat,c_param = 0.75)

summary(bgplm0.fit_known_c)

ggplot(data=bgplm0.fit_known_c$rating_curve) +
    geom_point(data=rc_dat,aes(Q,W)) +
    geom_line(aes(median,w)) +
    geom_line(aes(lower,w),linetype='dashed') +
    geom_line(aes(upper,w),linetype='dashed')


####bgplm#####
bgplm.fit<- bgplm(rc_formula,rc_dat)

summary(bgplm.fit)

ggplot(data=bgplm.fit$rating_curve) +
  geom_point(data=rc_dat,aes(Q,W)) +
  geom_line(aes(median,w)) +
  geom_line(aes(lower,w),linetype='dashed') +
  geom_line(aes(upper,w),linetype='dashed')

bgplm.fit_known_c <- bgplm(rc_formula,rc_dat,c_param = 0.75)

summary(bgplm.fit_known_c)

ggplot(data=bgplm.fit_known_c$rating_curve) +
    geom_point(data=rc_dat,aes(Q,W)) +
    geom_line(aes(median,w)) +
    geom_line(aes(lower,w),linetype='dashed') +
    geom_line(aes(upper,w),linetype='dashed')

#for debugging
data=rc_dat
formula <- rc_formula
model_dat <- data[,all.vars(formula)]
model_dat <- model_dat[order(model_dat[,2,drop=T]),]
Q <- model_dat[,1,drop=T]
w <- model_dat[,2,drop=T]
y=log(Q)
c_param=0.75
#c_param=NULL
w_max=NULL
forcepoint=rep(FALSE,nrow(data))
num_chains=4
nr_iter=20000
burnin=2000
thin=5




