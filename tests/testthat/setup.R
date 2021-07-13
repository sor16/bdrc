data(krokfors)
krokfors <- krokfors[(1:nrow(krokfors))%%4==0,]
h_extrap <- 10
known_c <- 7.8
set.seed(1)
plm0.fit <- plm0(Q~W,krokfors,parallel=F)
set.seed(1)
plm0.fit_known_c <- plm0(Q~W,krokfors,c_param=known_c,h_max=h_extrap,parallel=F)
set.seed(1)
plm.fit <- plm(Q~W,krokfors,parallel=F)
set.seed(1)
plm.fit_known_c <- plm(Q~W,krokfors,c_param=known_c,h_max=h_extrap,parallel=F)
set.seed(1)
gplm0.fit <- gplm0(Q~W,krokfors,parallel=F)
set.seed(1)
gplm0.fit_known_c <- gplm0(Q~W,krokfors,c_param=known_c,h_max=h_extrap,parallel=F)
set.seed(1)
gplm.fit <- gplm(Q~W,krokfors,parallel=F)
set.seed(1)
gplm.fit_known_c <- gplm(Q~W,krokfors,c_param=known_c,h_max=h_extrap,parallel=F)

t_obj <- tournament(plm0.fit,plm.fit,gplm0.fit,gplm.fit)
