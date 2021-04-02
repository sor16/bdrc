# dat <- read.csv('~/Downloads/KK_7.8.csv')
# mod_list <- vector('list',10)
# for(i in 1:10){
#     mod_list[[i]] <- bgplm0(Q~W,dat)
# }
# th_dat <- lapply(1:length(mod_list),function(i){
#               p_s <- mod_list[[i]]$param_summary
#               tibble(mod=i,param=rownames(p_s)[3:nrow(p_s)],value=p_s[3:nrow(p_s),'median'])
#           }) %>% bind_rows()
# th_dat %>% ggplot(aes(factor(param),value)) + geom_boxplot() + facet_wrap(~factor(param),scales='free')
#
#
# dic_dat <- lapply(1:length(mod_list),function(i){
#                 tibble(mod=i,deviance_med=mod_list[[i]]$Deviance_summary[,'median'],dic=mod_list[[i]]$DIC)
#             }) %>% bind_rows()
# ggplot(dic_dat,aes(deviance_med,dic)) + geom_point()
