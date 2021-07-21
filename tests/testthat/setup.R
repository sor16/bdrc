data(krokfors)
krokfors <- krokfors[(1:nrow(krokfors))%%4==0,]
h_extrap <- 10
known_c <- 7.8
set.seed(1)
plm0.fit <- plm0(Q~W,krokfors,parallel=F)
set.seed(1)
plm.fit <- plm(Q~W,krokfors,parallel=F)
set.seed(1)
gplm0.fit <- gplm0(Q~W,krokfors,parallel=F)
set.seed(1)
gplm.fit <- gplm(Q~W,krokfors,parallel=F)


t_obj <- tournament(plm0.fit,plm.fit,gplm0.fit,gplm.fit)
