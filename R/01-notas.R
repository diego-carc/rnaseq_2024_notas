# Make this project
usethis::create_project("~/rnaseq_2024_notas")

# Create this file
usethis::use_r("01-notas.R")

# Second activity
usethis::use_r("02-visualizar-mtcars.R")

# Init git repo
usethis::use_git()
usethis::use_github()

# Bioconductor notes
usethis::use_r("03-Bioconductor.R")

# Gene expression Bioconductor data structures
usethis::use_r("04-Exp-Bioc-Struct.R")

# Recount
usethis::use_r("05_recount3.R")

# Make my app
usethis::use_r("app.R")

# Explore Model Matrix
usethis::use_r("06-ExploreModelMtrx.R")

# Smoking mouse
usethis::use_r("smokingMice.R")

# Differential Expression analysis
usethis::use_r("07_limma.R")

# Last activity
usethis::use_r("08_notas_repaso.R")
