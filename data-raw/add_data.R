
norn <- read.csv('data-raw/norn.csv',header=T)
krokfors<- read.csv('data-raw/krokfors.csv',header=T)
spanga <- read.csv('data-raw/spanga.csv',header=T)
skogsliden <- read.csv('data-raw/skogsliden.csv',header=T)

usethis::use_data(norn, overwrite = TRUE)
usethis::use_data(krokfors, overwrite = TRUE)
usethis::use_data(spanga, overwrite = TRUE)
usethis::use_data(skogsliden, overwrite = TRUE)
