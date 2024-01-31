# # Summarized experiment
# This is a package that provides a SummarizedExperiment
# structure for the representation of RNAseq results.
#
# A SummarizedExperiment includes:
# # colData - > Info about the samples
# rowRanges -> Info about the genes
# assays -> The counts or number of reads
#
# We can access to these tables by functions

# The Metadata section is an optional field
# in the SummarizedExperiment struct

# We can access to the Features and Samples
# using the notation obj[rows, cols].

# The rows(Features) has options to filter
# by region of interest.

## Lets build our first SummarizedExperiment object
library("SummarizedExperiment")
## ?SummarizedExperiment

## De los ejemplos en la ayuda oficial

## Creamos los datos para nuestro objeto de tipo SummarizedExperiment
## para 200 genes a lo largo de 6 muestras
nrows <- 200
ncols <- 6
## Números al azar de cuentas
set.seed(20210223)
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
## Información de nuestros genes
rowRanges <- GRanges(
  rep(c("chr1", "chr2"), c(50, 150)),
  IRanges(floor(runif(200, 1e5, 1e6)), width = 100),
  strand = sample(c("+", "-"), 200, TRUE),
  feature_id = sprintf("ID%03d", 1:200)
)
names(rowRanges) <- paste0("gene_", seq_len(length(rowRanges)))
## Información de nuestras muestras
colData <- DataFrame(Treatment = rep(c("ChIP", "Input"), 3),
                     row.names = LETTERS[1:6])
## Juntamos ahora toda la información en un solo objeto de R
rse <- SummarizedExperiment(
  assays = SimpleList(counts = counts),
  rowRanges = rowRanges,
  colData = colData
)

## Exploremos el objeto resultante
rse

# Dimension
dim(rse)

#
assayNames(rse)

#
assay(rse, 1)

#
colData(rse)

rse$Treatment

# Exercise
## Comando 1
rse[1:2, ]
# This command acces to the first two features across all samples

rse[, c("A", "D", "F")]
# This command acces to all the features across samples A D and F

# iSEE
# This is a package generates an interactive web-page for gene Expression
#analisys
#Load package
library("iSEE")
library("spatialLIBD")

# Download data
sce_layer <- spatialLIBD::fetch_data("sce_layer")

# Open iSEE

iSEE::iSEE(sce_layer)

# We generate the images:
# "figures/ReducedDimensionPlot1.pdf"
# "figures/ComplexHeatmapPlot1.pdf"
#



