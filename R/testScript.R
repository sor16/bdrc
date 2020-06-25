library(readxl)
library(dplyr)
library(ggplot2)
source('R/Adist.R')
source('R/B_splines.R')
source('R/W_unobserved.R')
source('R/priors.R')
source('R/bgplm.R')
source('R/bplm.R')
source('R/plm_methods.R')

#read data
rc_dat <- read_excel('data/exceldata_RC.xlsx') %>% mutate(W=0.01*W)
rc_formula <- as.formula('Q~W')

## Model 0 (bplm0)
bplm0.fit <- bplm0(rc_formula,rc_dat)

summary(bplm0.fit)


ggplot(data=bplm0.fit$rating_curve) +
  geom_point(data=rc_dat,aes(Q,W)) +
  geom_line(aes(median,W)) +
  geom_line(aes(lower,W),linetype='dashed') +
  geom_line(aes(upper,W),linetype='dashed')


bplm0.fit_known_c <- bplm0(rc_formula,rc_dat,c_param=0.75)

summary(bplm0.fit_known_c)

ggplot(data=bplm0.fit_known_c$rating_curve) +
  geom_point(data=rc_dat,aes(Q,W)) +
  geom_line(aes(median,W)) +
  geom_line(aes(lower,W),linetype='dashed') +
  geom_line(aes(upper,W),linetype='dashed')

#Model 2
bgplm.fit <- bgplm(rc_formula,rc_dat)

summary(bgplm.fit)

ggplot(data=bgplm.fit$rating_curve) +
    geom_point(data=rc_dat,aes(Q,W)) +
    geom_line(aes(median,W)) +
    geom_line(aes(lower,W),linetype='dashed') +
    geom_line(aes(upper,W),linetype='dashed')

bgplm.fit_known_c <- bgplm(rc_formula,rc_dat,c_param = 0.75)

summary(bgplm.fit_known_c)

ggplot(data=bgplm.fit_known_c$rating_curve) +
    geom_point(data=rc_dat,aes(Q,W)) +
    geom_line(aes(median,W)) +
    geom_line(aes(lower,W),linetype='dashed') +
    geom_line(aes(upper,W),linetype='dashed')

## Model 1 (bplm)
bplm.fit <- bplm(rc_formula,rc_dat)

summary(bplm.fit)

ggplot(data=bplm.fit$rating_curve) +
    geom_point(data=rc_dat,aes(Q,W)) +
    geom_line(aes(median,W)) +
    geom_line(aes(lower,W),linetype='dashed') +
    geom_line(aes(upper,W),linetype='dashed')

bplm.fit_known_c <- bplm(rc_formula,rc_dat,c_param=0.75)

summary(bplm.fit_known_c)

ggplot(data=bplm.fit_known_c$rating_curve) +
    geom_point(data=rc_dat,aes(Q,W)) +
    geom_line(aes(median,W)) +
    geom_line(aes(lower,W),linetype='dashed') +
    geom_line(aes(upper,W),linetype='dashed')


