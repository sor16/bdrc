library(readxl)
library(dplyr)
library(ggplot2)
source('R/Adist.R')
source('R/B_splines.R')
source('R/W_unobserved.R')
source('R/priors.R')
source('R/bgplm.R')
rc_dat <- read_excel('data/exceldata_RC.xlsx') %>% mutate(W=0.01*W)
rc_formula <- as.formula('Q~W')
bgplm.fit <- bgplm(rc_formula,rc_dat)

summary(bgplm.fit)

ggplot(data=bgplm.fit$rating_curve) +
    geom_point(data=rc_dat,aes(Q,W)) +
    geom_line(aes(median,W)) +
    geom_line(aes(lower,W),linetype='dashed') +
    geom_line(aes(upper,W),linetype='dashed')

ggplot(data=filter(bgplm.fit$rating_curve,W>=min(rc_dat$W) & W<=max(rc_dat$W))) +
    geom_point(data=rc_dat,aes(Q,W)) +
    geom_line(aes(median,W)) +
    geom_line(aes(lower,W),linetype='dashed') +
    geom_line(aes(upper,W),linetype='dashed')

bgplm.fit_known_c <- bgplm(rc_formula,rc_dat,c_param = 0.75)
