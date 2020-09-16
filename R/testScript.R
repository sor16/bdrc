library(readxl)
library(dplyr)
library(ggplot2)
source('R/B_splines.R')
source('R/W_unobserved.R')
source('R/help_functions.R')
source('R/plm_methods.R')
source('R/bplm0.R')
source('R/bplm.R')
source('R/bgplm0.R')
source('R/bgplm.R')


#read data
rc_dat <- read_excel('data/exceldata_RC.xlsx') %>% mutate(W=0.01*W)
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
    geom_point(data=rc_dat,aes(Q,w)) +
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





