V316_river <- read.table('data-raw/V316_river.tsv',header=T)
V316_river$W <- V316_river$W*0.01

usethis::use_data(V316_river, overwrite = TRUE)
