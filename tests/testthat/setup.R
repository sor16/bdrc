data(krokfors)
krokfors <- krokfors[(1:nrow(krokfors)) %% 4 == 0, ]
h <- krokfors$W
y <- log(krokfors$Q)
h_extrap <- 10
known_c <- 7.8

set.seed(1)
plm0.fit <- plm0(Q ~ W, krokfors, num_cores = 2)
set.seed(1)
plm.fit <- plm(Q ~ W, krokfors, num_cores = 2)
set.seed(1)
gplm0.fit <- gplm0(Q ~ W, krokfors, num_cores = 2)
set.seed(1)
gplm.fit <- gplm(Q ~ W, krokfors, num_cores = 2)

t_obj <- tournament(list(plm0.fit, plm.fit, gplm0.fit, gplm.fit))

# saveRDS(plm0.fit, file = "tests/cached_results/plm0.fit.rds")
# saveRDS(plm.fit, file = "tests/cached_results/plm.fit.rds")
# saveRDS(gplm0.fit, file = "tests/cached_results/gplm0.fit.rds")
# saveRDS(gplm.fit, file = "tests/cached_results/gplm.fit.rds")
# saveRDS(t_obj, file = "tests/cached_results/tournament.rds")
