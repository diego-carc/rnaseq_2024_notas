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
