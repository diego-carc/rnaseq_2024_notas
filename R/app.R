## Para hacer
## https://libd.shinyapps.io/SRP009615/

## Primero necesitamos configurar RStudio con
## shinyapps.io. Para eso necesitaremos:
# install.packages("rsconnect")

## También necesitamos verificar que tengamos todos
## los paquetes en versiones nuevas. Eso
## lo podemos hacer con:
# BiocManager::valid()

## Después necesitamos copiar y pegar la información
## de nuestra cuenta (numéro y token de acceso)

## Ahora si ya podemos continuar
options(repos = BiocManager::repositories())

# Load package
library(recount3)

# Check available packages
human_projects <- available_projects()

# Check class
class(human_projects)
dim(human_projects)
head(human_projects)

## Encuentra tu proyecto de interés. Aquí usaremos
## SRP009615 de ejemplo

proj_info <- subset(
  human_projects,
  project == "SRP009615" & project_type == "data_sources"
)

## Crea un objeto de tipo RangedSummarizedExperiment (RSE)
## con la información a nivel de genes
rse_gene_SRP009615 <- create_rse(proj_info)

## Explora los proyectos disponibles de forma interactiva
#proj_info_interactive <- interactiveDisplayBase::display(human_projects)
## Selecciona un solo renglón en la tabla y da click en "send".

## Convirtamos las cuentas por nucleotido a cuentas por lectura
## usando compute_read_counts().
## Para otras transformaciones como RPKM y TPM, revisa transform_counts().
assay(rse_gene_SRP009615, "counts") <- compute_read_counts(rse_gene_SRP009615)


## Para este estudio en específico, hagamos más fácil de usar la
## información del experimento
rse_gene_SRP009615 <- expand_sra_attributes(rse_gene_SRP009615)
colData(rse_gene_SRP009615)[
  ,
  grepl("^sra_attribute", colnames(colData(rse_gene_SRP009615)))
]



iSEE::iSEE(rse_gene_SRP009615)
