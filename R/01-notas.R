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

# Gene expression Bioconductor data structures
usethis::use_r("04-Exp-Bioc-Struct.R")
