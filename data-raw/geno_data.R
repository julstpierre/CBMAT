# load required packages ----
if (!require("pacman")) install.packages("pacman") 
pacman::p_load(usethis)

#write data
geno <- t(read.table("data-raw/G.012"))
usethis::use_data(geno, overwrite = TRUE)
