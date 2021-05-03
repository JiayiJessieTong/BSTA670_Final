## Create R package

# 1. Check if the name of the package is available
library(available)
available::available("extODAL")

library(usethis)
usethis::create_package("~/Documents/GitHub/extODAL")
usethis::use_readme_rmd()

usethis::use_git()
