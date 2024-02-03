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

# iSEE config
initial <- list()

################################################################################
# Settings for Reduced dimension plot 1
################################################################################

initial[["ReducedDimensionPlot1"]] <- new("ReducedDimensionPlot", Type = "PCA", XAxis = 1L, YAxis = 2L,
                                          FacetRowByColData = "sample_name", FacetColumnByColData = "sample_name",
                                          ColorByColumnData = "layer_guess", ColorByFeatureNameAssay = "logcounts",
                                          ColorBySampleNameColor = "#FF0000", ShapeByColumnData = "sample_name",
                                          SizeByColumnData = NA_character_, TooltipColumnData = character(0),
                                          FacetRowBy = "None", FacetColumnBy = "None", ColorBy = "Column data",
                                          ColorByDefaultColor = "#000000", ColorByFeatureName = "ENSG00000243485",
                                          ColorByFeatureSource = "---", ColorByFeatureDynamicSource = FALSE,
                                          ColorBySampleName = "151507_Layer1", ColorBySampleSource = "---",
                                          ColorBySampleDynamicSource = FALSE, ShapeBy = "None", SizeBy = "None",
                                          SelectionAlpha = 1, ZoomData = numeric(0), BrushData = list(),
                                          VisualBoxOpen = TRUE, VisualChoices = c("Color", "Size"),
                                          ContourAdd = FALSE, ContourColor = "#0000FF", PointSize = 5,
                                          PointAlpha = 1, Downsample = FALSE, DownsampleResolution = 200,
                                          CustomLabels = FALSE, CustomLabelsText = "151507_Layer1",
                                          FontSize = 1, LegendPointSize = 1, LegendPosition = "Bottom",
                                          HoverInfo = TRUE, LabelCenters = FALSE, LabelCentersBy = "sample_name",
                                          LabelCentersColor = "#000000", VersionInfo = list(iSEE = structure(list(
                                            c(2L, 14L, 0L)), class = c("package_version", "numeric_version"
                                            ))), PanelId = c(ReducedDimensionPlot = 1L), PanelHeight = 500L,
                                          PanelWidth = 4L, SelectionBoxOpen = FALSE, RowSelectionSource = "---",
                                          ColumnSelectionSource = "---", DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE,
                                          ColumnSelectionDynamicSource = FALSE, RowSelectionRestrict = FALSE,
                                          ColumnSelectionRestrict = FALSE, SelectionHistory = list())

################################################################################
# Settings for Row data table 1
################################################################################

initial[["RowDataTable1"]] <- new("RowDataTable", Selected = "ENSG00000197971", Search = "MBP;",
                                  SearchColumns = c("", "", "", "", "", "", "", "", "", ""),
                                  HiddenColumns = character(0), VersionInfo = list(iSEE = structure(list(
                                    c(2L, 14L, 0L)), class = c("package_version", "numeric_version"
                                    ))), PanelId = c(RowDataTable = 1L), PanelHeight = 500L,
                                  PanelWidth = 4L, SelectionBoxOpen = FALSE, RowSelectionSource = "---",
                                  ColumnSelectionSource = "---", DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE,
                                  ColumnSelectionDynamicSource = FALSE, RowSelectionRestrict = FALSE,
                                  ColumnSelectionRestrict = FALSE, SelectionHistory = list())

################################################################################
# Settings for Feature assay plot 1
################################################################################

initial[["FeatureAssayPlot1"]] <- new("FeatureAssayPlot", Assay = "logcounts", XAxis = "None",
                                      XAxisColumnData = "sample_name", XAxisFeatureName = "ENSG00000243485",
                                      XAxisFeatureSource = "---", XAxisFeatureDynamicSource = FALSE,
                                      YAxisFeatureName = "ENSG00000243485", YAxisFeatureSource = "---",
                                      YAxisFeatureDynamicSource = FALSE, FacetRowByColData = "sample_name",
                                      FacetColumnByColData = "sample_name", ColorByColumnData = "sample_name",
                                      ColorByFeatureNameAssay = "logcounts", ColorBySampleNameColor = "#FF0000",
                                      ShapeByColumnData = "sample_name", SizeByColumnData = NA_character_,
                                      TooltipColumnData = character(0), FacetRowBy = "None", FacetColumnBy = "None",
                                      ColorBy = "None", ColorByDefaultColor = "#000000", ColorByFeatureName = "ENSG00000243485",
                                      ColorByFeatureSource = "---", ColorByFeatureDynamicSource = FALSE,
                                      ColorBySampleName = "151507_Layer1", ColorBySampleSource = "---",
                                      ColorBySampleDynamicSource = FALSE, ShapeBy = "None", SizeBy = "None",
                                      SelectionAlpha = 0.1, ZoomData = numeric(0), BrushData = list(),
                                      VisualBoxOpen = FALSE, VisualChoices = "Color", ContourAdd = FALSE,
                                      ContourColor = "#0000FF", PointSize = 1, PointAlpha = 1,
                                      Downsample = FALSE, DownsampleResolution = 200, CustomLabels = FALSE,
                                      CustomLabelsText = "151507_Layer1", FontSize = 1, LegendPointSize = 1,
                                      LegendPosition = "Bottom", HoverInfo = TRUE, LabelCenters = FALSE,
                                      LabelCentersBy = "sample_name", LabelCentersColor = "#000000",
                                      VersionInfo = list(iSEE = structure(list(c(2L, 14L, 0L)), class = c("package_version",
                                                                                                          "numeric_version"))), PanelId = c(FeatureAssayPlot = 1L),
                                      PanelHeight = 500L, PanelWidth = 4L, SelectionBoxOpen = FALSE,
                                      RowSelectionSource = "---", ColumnSelectionSource = "---",
                                      DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE,
                                      RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE,
                                      SelectionHistory = list())

################################################################################
# Settings for Column data plot 1
################################################################################

initial[["ColumnDataPlot1"]] <- new("ColumnDataPlot", XAxis = "None", YAxis = "sample_name",
                                    XAxisColumnData = "sample_name", FacetRowByColData = "sample_name",
                                    FacetColumnByColData = "sample_name", ColorByColumnData = "sample_name",
                                    ColorByFeatureNameAssay = "logcounts", ColorBySampleNameColor = "#FF0000",
                                    ShapeByColumnData = "sample_name", SizeByColumnData = NA_character_,
                                    TooltipColumnData = character(0), FacetRowBy = "None", FacetColumnBy = "None",
                                    ColorBy = "None", ColorByDefaultColor = "#000000", ColorByFeatureName = "ENSG00000243485",
                                    ColorByFeatureSource = "---", ColorByFeatureDynamicSource = FALSE,
                                    ColorBySampleName = "151507_Layer1", ColorBySampleSource = "---",
                                    ColorBySampleDynamicSource = FALSE, ShapeBy = "None", SizeBy = "None",
                                    SelectionAlpha = 0.1, ZoomData = numeric(0), BrushData = list(),
                                    VisualBoxOpen = FALSE, VisualChoices = "Color", ContourAdd = FALSE,
                                    ContourColor = "#0000FF", PointSize = 1, PointAlpha = 1,
                                    Downsample = FALSE, DownsampleResolution = 200, CustomLabels = FALSE,
                                    CustomLabelsText = "151507_Layer1", FontSize = 1, LegendPointSize = 1,
                                    LegendPosition = "Bottom", HoverInfo = TRUE, LabelCenters = FALSE,
                                    LabelCentersBy = "sample_name", LabelCentersColor = "#000000",
                                    VersionInfo = list(iSEE = structure(list(c(2L, 14L, 0L)), class = c("package_version",
                                                                                                        "numeric_version"))), PanelId = c(ColumnDataPlot = 1L), PanelHeight = 500L,
                                    PanelWidth = 4L, SelectionBoxOpen = FALSE, RowSelectionSource = "---",
                                    ColumnSelectionSource = "---", DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE,
                                    ColumnSelectionDynamicSource = FALSE, RowSelectionRestrict = FALSE,
                                    ColumnSelectionRestrict = FALSE, SelectionHistory = list())

################################################################################
# Settings for Row data plot 1
################################################################################

initial[["RowDataPlot1"]] <- new("RowDataPlot", XAxis = "None", YAxis = "source", XAxisRowData = "source",
                                 FacetRowByRowData = "source", FacetColumnByRowData = "source",
                                 ColorByRowData = "source", ColorBySampleNameAssay = "logcounts",
                                 ColorByFeatureNameColor = "#FF0000", ShapeByRowData = "source",
                                 SizeByRowData = NA_character_, TooltipRowData = character(0),
                                 FacetRowBy = "None", FacetColumnBy = "None", ColorBy = "None",
                                 ColorByDefaultColor = "#000000", ColorByFeatureName = "ENSG00000243485",
                                 ColorByFeatureSource = "---", ColorByFeatureDynamicSource = FALSE,
                                 ColorBySampleName = "151507_Layer1", ColorBySampleSource = "---",
                                 ColorBySampleDynamicSource = FALSE, ShapeBy = "None", SizeBy = "None",
                                 SelectionAlpha = 0.1, ZoomData = numeric(0), BrushData = list(),
                                 VisualBoxOpen = FALSE, VisualChoices = "Color", ContourAdd = FALSE,
                                 ContourColor = "#0000FF", PointSize = 1, PointAlpha = 1,
                                 Downsample = FALSE, DownsampleResolution = 200, CustomLabels = FALSE,
                                 CustomLabelsText = "ENSG00000243485", FontSize = 1, LegendPointSize = 1,
                                 LegendPosition = "Bottom", HoverInfo = TRUE, LabelCenters = FALSE,
                                 LabelCentersBy = "source", LabelCentersColor = "#000000",
                                 VersionInfo = list(iSEE = structure(list(c(2L, 14L, 0L)), class = c("package_version",
                                                                                                     "numeric_version"))), PanelId = c(RowDataPlot = 1L), PanelHeight = 500L,
                                 PanelWidth = 4L, SelectionBoxOpen = FALSE, RowSelectionSource = "---",
                                 ColumnSelectionSource = "---", DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE,
                                 ColumnSelectionDynamicSource = FALSE, RowSelectionRestrict = FALSE,
                                 ColumnSelectionRestrict = FALSE, SelectionHistory = list())

################################################################################
# Settings for Sample assay plot 1
################################################################################

initial[["SampleAssayPlot1"]] <- new("SampleAssayPlot", Assay = "logcounts", XAxis = "None", XAxisRowData = "source",
                                     XAxisSampleName = "151507_Layer1", XAxisSampleSource = "---",
                                     XAxisSampleDynamicSource = FALSE, YAxisSampleName = "151507_Layer1",
                                     YAxisSampleSource = "---", YAxisSampleDynamicSource = FALSE,
                                     FacetRowByRowData = "source", FacetColumnByRowData = "source",
                                     ColorByRowData = "source", ColorBySampleNameAssay = "logcounts",
                                     ColorByFeatureNameColor = "#FF0000", ShapeByRowData = "source",
                                     SizeByRowData = NA_character_, TooltipRowData = character(0),
                                     FacetRowBy = "None", FacetColumnBy = "None", ColorBy = "None",
                                     ColorByDefaultColor = "#000000", ColorByFeatureName = "ENSG00000243485",
                                     ColorByFeatureSource = "---", ColorByFeatureDynamicSource = FALSE,
                                     ColorBySampleName = "151507_Layer1", ColorBySampleSource = "---",
                                     ColorBySampleDynamicSource = FALSE, ShapeBy = "None", SizeBy = "None",
                                     SelectionAlpha = 0.1, ZoomData = numeric(0), BrushData = list(),
                                     VisualBoxOpen = FALSE, VisualChoices = "Color", ContourAdd = FALSE,
                                     ContourColor = "#0000FF", PointSize = 1, PointAlpha = 1,
                                     Downsample = FALSE, DownsampleResolution = 200, CustomLabels = FALSE,
                                     CustomLabelsText = "ENSG00000243485", FontSize = 1, LegendPointSize = 1,
                                     LegendPosition = "Bottom", HoverInfo = TRUE, LabelCenters = FALSE,
                                     LabelCentersBy = "source", LabelCentersColor = "#000000",
                                     VersionInfo = list(iSEE = structure(list(c(2L, 14L, 0L)), class = c("package_version",
                                                                                                         "numeric_version"))), PanelId = c(SampleAssayPlot = 1L),
                                     PanelHeight = 500L, PanelWidth = 4L, SelectionBoxOpen = FALSE,
                                     RowSelectionSource = "---", ColumnSelectionSource = "---",
                                     DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE,
                                     RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE,
                                     SelectionHistory = list())

################################################################################
# Settings for Column data table 1
################################################################################

initial[["ColumnDataTable1"]] <- new("ColumnDataTable", Selected = "151507_Layer1", Search = "",
                                     SearchColumns = c("", "", "", "", "", "", "", "", "", "",
                                                       "", "", ""), HiddenColumns = character(0), VersionInfo = list(
                                                         iSEE = structure(list(c(2L, 14L, 0L)), class = c("package_version",
                                                                                                          "numeric_version"))), PanelId = c(ColumnDataTable = 1L),
                                     PanelHeight = 500L, PanelWidth = 4L, SelectionBoxOpen = FALSE,
                                     RowSelectionSource = "---", ColumnSelectionSource = "---",
                                     DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE,
                                     RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE,
                                     SelectionHistory = list())

################################################################################
# Settings for Complex heatmap 1
################################################################################

initial[["ComplexHeatmapPlot1"]] <- new("ComplexHeatmapPlot", Assay = "logcounts", CustomRows = TRUE,
                                        CustomRowsText = "ENSG00000183036\nENSG00000168314\nENSG00000197971",
                                        ClusterRows = TRUE, ClusterRowsDistance = "spearman", ClusterRowsMethod = "ward.D2",
                                        DataBoxOpen = FALSE, VisualChoices = "Annotations", ColumnData = "layer_guess",
                                        RowData = character(0), CustomBounds = FALSE, LowerBound = NA_real_,
                                        UpperBound = NA_real_, AssayCenterRows = FALSE, AssayScaleRows = FALSE,
                                        DivergentColormap = "purple < black < yellow", ShowDimNames = "Rows",
                                        LegendPosition = "Bottom", LegendDirection = "Horizontal",
                                        VisualBoxOpen = FALSE, NamesRowFontSize = 10, NamesColumnFontSize = 10,
                                        ShowColumnSelection = TRUE, OrderColumnSelection = TRUE,
                                        VersionInfo = list(iSEE = structure(list(c(2L, 14L, 0L)), class = c("package_version",
                                                                                                            "numeric_version"))), PanelId = c(ComplexHeatmapPlot = 1L),
                                        PanelHeight = 500L, PanelWidth = 4L, SelectionBoxOpen = FALSE,
                                        RowSelectionSource = "---", ColumnSelectionSource = "---",
                                        RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE,
                                        RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE,
                                        SelectionHistory = list())
# Open iSEE

iSEE::iSEE(sce_layer, initial = initial)

# We generate the images:
# "figures/ReducedDimensionPlot1.pdf"
# "figures/ComplexHeatmapPlot1.pdf"
#



