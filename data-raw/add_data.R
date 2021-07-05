
bunnerviken <- read.csv('data-raw/BUNNERVIKEN.csv',header=T)
flotemarken<- read.csv('data-raw/FLOTEMARKEN.csv',header=T)
lisjobacken <- read.csv('data-raw/LISJOBACKEN.csv',header=T)
halla <- read.csv('data-raw/HALLA.csv',header=T)

usethis::use_data(bunnerviken, overwrite = TRUE)
usethis::use_data(flotemarken, overwrite = TRUE)
usethis::use_data(lisjobacken, overwrite = TRUE)
usethis::use_data(halla, overwrite = TRUE)
