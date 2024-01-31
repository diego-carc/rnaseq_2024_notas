library("sessioninfo")
library("here")
library("ggplot2")

## Hello world
print("Hi, I'm Diego!")

## Directories
dir_plots <- here::here("figures")
dir_rdata <- here::here("processed-data")

## Make directories for figures and files
dir.create(dir_plots, showWarnings = FALSE)
dir.create(dir_rdata, showWarnings = FALSE)

## Make an example image
pdf(file.path(dir_plots, "mtcars_gear_vs_mpg.pdf"),
    useDingbats = FALSE
)
ggplot(mtcars, aes(group = gear, y = mpg)) +
  geom_boxplot()
dev.off()

## Reproduce the code
options(width = 120)
sessioninfo::session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.2 (2023-10-31 ucrt)
# os       Windows 11 x64 (build 22621)
# system   x86_64, mingw32
# ui       RStudio
# language (EN)
# collate  Spanish_Mexico.utf8
# ctype    Spanish_Mexico.utf8
# tz       America/Mexico_City
# date     2024-01-30
# rstudio  2023.12.0+369 Ocean Storm (desktop)
# pandoc   NA
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package     * version date (UTC) lib source
# askpass       1.2.0   2023-09-03 [1] CRAN (R 4.3.2)
# cli           3.6.2   2023-12-11 [1] CRAN (R 4.3.2)
# colorspace    2.1-0   2023-01-23 [1] CRAN (R 4.3.1)
# crayon        1.5.2   2022-09-29 [1] CRAN (R 4.3.1)
# credentials   2.0.1   2023-09-06 [1] CRAN (R 4.3.2)
# curl          5.2.0   2023-12-08 [1] CRAN (R 4.3.2)
# dplyr         1.1.4   2023-11-17 [1] CRAN (R 4.3.2)
# fansi         1.0.6   2023-12-08 [1] CRAN (R 4.3.2)
# farver        2.1.1   2022-07-06 [1] CRAN (R 4.3.1)
# fs            1.6.3   2023-07-20 [1] CRAN (R 4.3.1)
# generics      0.1.3   2022-07-05 [1] CRAN (R 4.3.1)
# gert          2.0.1   2023-12-04 [1] CRAN (R 4.3.2)
# ggplot2     * 3.4.4   2023-10-12 [1] CRAN (R 4.3.2)
# gh            1.4.0   2023-02-22 [1] CRAN (R 4.3.2)
# gitcreds      0.1.2   2022-09-08 [1] CRAN (R 4.3.2)
# glue          1.7.0   2024-01-09 [1] CRAN (R 4.3.2)
# gtable        0.3.4   2023-08-21 [1] CRAN (R 4.3.2)
# here        * 1.0.1   2020-12-13 [1] CRAN (R 4.3.2)
# httr2         1.0.0   2023-11-14 [1] CRAN (R 4.3.2)
# jsonlite      1.8.8   2023-12-04 [1] CRAN (R 4.3.2)
# labeling      0.4.3   2023-08-29 [1] CRAN (R 4.3.1)
# lifecycle     1.0.4   2023-11-07 [1] CRAN (R 4.3.2)
# magrittr      2.0.3   2022-03-30 [1] CRAN (R 4.3.1)
# munsell       0.5.0   2018-06-12 [1] CRAN (R 4.3.1)
# openssl       2.1.1   2023-09-25 [1] CRAN (R 4.3.2)
# pillar        1.9.0   2023-03-22 [1] CRAN (R 4.3.1)
# pkgconfig     2.0.3   2019-09-22 [1] CRAN (R 4.3.1)
# purrr         1.0.2   2023-08-10 [1] CRAN (R 4.3.2)
# R6            2.5.1   2021-08-19 [1] CRAN (R 4.3.1)
# rappdirs      0.3.3   2021-01-31 [1] CRAN (R 4.3.1)
# rlang         1.1.3   2024-01-10 [1] CRAN (R 4.3.2)
# rprojroot     2.0.4   2023-11-05 [1] CRAN (R 4.3.2)
# rstudioapi    0.15.0  2023-07-07 [1] CRAN (R 4.3.1)
# scales        1.3.0   2023-11-28 [1] CRAN (R 4.3.2)
# sessioninfo * 1.2.2   2021-12-06 [1] CRAN (R 4.3.2)
# sys           3.4.2   2023-05-23 [1] CRAN (R 4.3.1)
# tibble        3.2.1   2023-03-20 [1] CRAN (R 4.3.1)
# tidyselect    1.2.0   2022-10-10 [1] CRAN (R 4.3.1)
# usethis       2.2.2   2023-07-06 [1] CRAN (R 4.3.2)
# utf8          1.2.4   2023-10-22 [1] CRAN (R 4.3.2)
# vctrs         0.6.5   2023-12-01 [1] CRAN (R 4.3.2)
# withr         3.0.0   2024-01-16 [1] CRAN (R 4.3.2)
#
# [1] C:/Users/Diego Carmona Campos/AppData/Local/R/win-library/4.3
# [2] C:/Program Files/R/R-4.3.2/library
