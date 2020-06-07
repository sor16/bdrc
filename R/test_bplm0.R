library(readxl)
library(dplyr)
source('R/priors.R')
source('R/bplm0.R')
rc_dat <- read_excel('data/exceldata_RC.xlsx') %>% mutate(W=0.01*W)
rc_formula <- as.formula('Q~W')
bplm0.fit <- bplm0(rc_formula,rc_dat)

summary(bplm0.fit)

ggplot(data=bplm0.fit$rating_curve) +
  geom_point(data=rc_dat,aes(Q,W)) +
  geom_line(aes(median,W)) +
  geom_line(aes(lower,W),linetype='dashed') +
  geom_line(aes(upper,W),linetype='dashed')

ggplot(data=filter(bplm0.fit$rating_curve,W>=min(rc_dat$W) & W<=max(rc_dat$W))) +
  geom_point(data=rc_dat,aes(Q,W)) +
  geom_line(aes(median,W)) +
  geom_line(aes(lower,W),linetype='dashed') +
  geom_line(aes(upper,W),linetype='dashed')

bgplm0.fit_known_c <- bplm0(rc_formula,rc_dat,c_param = 0.75)
