# library(readr)
# GREA_dat=read_csv('~/Downloads/GREA.csv')
# MELBY_dat=read_csv('~/Downloads/MELBY.csv')
# LIFFEDARVE_dat <- read_csv('~/Downloads/LIFFEDARVE.csv')
# ## GREA was reported not to run with bplm. Fails to replicate. Check whether a fixed c was used
# bplm.GREA.fit <- bplm(Q~W,GREA_dat)
# bplm.GREA.fit.c <- bplm(Q~W,GREA_dat,c_param=7.419)
# ## error
# bplm.GREA.fit.c <- bplm(Q~W,GREA_dat,c_param=7.28)
#
# ## LIFFEDARVE was also reported not to run with bplm. Fails to replicate. Check whether a fixed c was used
# bplm.LIFFEDARVE.fit <- bplm(Q~W,LIFFEDARVE_dat)
# bplm.LIFFEDARVE.fit.c <- bplm(Q~W,LIFFEDARVE_dat,c_param=9.878)
# bplm.LIFFEDARVE.fit.c <- bplm(Q~W,LIFFEDARVE_dat,c_param=9.4)
#
#
# ## LIFFEDARVE was also reported not to run with bgplm. Fails to replicate. Check whether a fixed c was used
# bgplm.MELBY.fit <- bgplm(Q~W,MELBY_dat)
#
# bgplm.MELBY.fit.c <- bgplm(Q~W,MELBY_dat,c_param=7.617)
# ##error
# bgplm.MELBY.fit.c <- bgplm(Q~W,MELBY_dat,c_param=7.47)
#
# ##Debugging of fixed c model
# bplm0.fit <- bplm0(Q~W,V316_river,c_param=0.75)
# bplm.fit <- bplm(Q~W,V316_river,c_param=0.75)
# bgplm0.fit <- bgplm0(Q~W,V316_river,c_param=0.75)
# bgplm.fit <- bgplm(Q~W,V316_river,c_param=0.75)
# # #
