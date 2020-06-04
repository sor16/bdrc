library(readxl)
library(dplyr)
library(ggplot2)
source('R/Adist.R')
source('R/B_splines.R')
source('R/W_unobserved.R')
source('R/priors.R')
source('R/bgplm.R')
source('R/plm_methods.R')


rc_dat <- read_excel('data/exceldata_RC.xlsx') %>% mutate(W=0.01*W)
rc_formula <- as.formula('Q~W')
bgplm.fit <- bgplm(rc_formula,rc_dat)

summary(bgplm.fit)

ggplot(data=bgplm.fit$rating_curve) +
    geom_point(data=rc_dat,aes(Q,W)) +
    geom_line(aes(median,W)) +
    geom_line(aes(lower,W),linetype='dashed') +
    geom_line(aes(upper,W),linetype='dashed')

bgplm.fit_known_c <- bgplm(rc_formula,rc_dat,c_param = 1.199)

summary(bgplm.fit_known_c)

ggplot(data=bgplm.fit_known_c$rating_curve) +
    geom_point(data=rc_dat,aes(Q,W)) +
    geom_line(aes(median,W)) +
    geom_line(aes(lower,W),linetype='dashed') +
    geom_line(aes(upper,W),linetype='dashed')
