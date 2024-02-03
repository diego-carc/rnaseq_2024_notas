## ----download_SRP045638----------------
library("recount3")

human_projects <- available_projects()

rse_gene_SRP045638 <- create_rse(
  subset(
    human_projects,
    project == "SRP045638" & project_type == "data_sources"
  )
)
assay(rse_gene_SRP045638, "counts") <- compute_read_counts(rse_gene_SRP045638)

## ----describe_issue--------------------
rse_gene_SRP045638$sra.sample_attributes[1:3]

## ----solve_issue-----------------------
rse_gene_SRP045638$sra.sample_attributes <- gsub("dev_stage;;Fetal\\|", "", rse_gene_SRP045638$sra.sample_attributes)
rse_gene_SRP045638$sra.sample_attributes[1:3]

## ----describe_issue--------------------
rse_gene_SRP045638$sra.sample_attributes[1:3]

## ----attributes------------------------
rse_gene_SRP045638 <- expand_sra_attributes(rse_gene_SRP045638)

colData(rse_gene_SRP045638)[
  ,
  grepl("^sra_attribute", colnames(colData(rse_gene_SRP045638)))
]

## ----re_cast---------------------------
## Pasar de character a numeric o factor
rse_gene_SRP045638$sra_attribute.age <- as.numeric(rse_gene_SRP045638$sra_attribute.age)
rse_gene_SRP045638$sra_attribute.disease <- factor(tolower(rse_gene_SRP045638$sra_attribute.disease))
rse_gene_SRP045638$sra_attribute.RIN <- as.numeric(rse_gene_SRP045638$sra_attribute.RIN)
rse_gene_SRP045638$sra_attribute.sex <- factor(rse_gene_SRP045638$sra_attribute.sex)

## Resumen de las variables de interés
summary(as.data.frame(colData(rse_gene_SRP045638)[
  ,
  grepl("^sra_attribute.[age|disease|RIN|sex]", colnames(colData(rse_gene_SRP045638)))
]))

## Encontraremos diferencias entre muestra prenatalas vs postnatales
rse_gene_SRP045638$prenatal <- factor(ifelse(rse_gene_SRP045638$sra_attribute.age < 0, "prenatal", "postnatal"))
table(rse_gene_SRP045638$prenatal)

## http://rna.recount.bio/docs/quality-check-fields.html
rse_gene_SRP045638$assigned_gene_prop <- rse_gene_SRP045638$recount_qc.gene_fc_count_all.assigned / rse_gene_SRP045638$recount_qc.gene_fc_count_all.total
summary(rse_gene_SRP045638$assigned_gene_prop)

with(colData(rse_gene_SRP045638), plot(assigned_gene_prop, sra_attribute.RIN))

## Hm... veamos si hay una diferencia entre los grupos
with(colData(rse_gene_SRP045638), tapply(assigned_gene_prop, prenatal, summary))

## Guardemos nuestro objeto entero por si luego cambiamos de opinión
rse_gene_SRP045638_unfiltered <- rse_gene_SRP045638

## Eliminemos a muestras malas
hist(rse_gene_SRP045638$assigned_gene_prop)

table(rse_gene_SRP045638$assigned_gene_prop < 0.3)


rse_gene_SRP045638 <- rse_gene_SRP045638[, rse_gene_SRP045638$assigned_gene_prop > 0.3]

## Calculemos los niveles medios de expresión de los genes en nuestras
## muestras.
## Ojo: en un análisis real probablemente haríamos esto con los RPKMs o CPMs
## en vez de las cuentas.
## En realidad usariamos:
# edgeR::filterByExpr() https://bioconductor.org/packages/edgeR/ https://rdrr.io/bioc/edgeR/man/filterByExpr.html
# genefilter::genefilter() https://bioconductor.org/packages/genefilter/ https://rdrr.io/bioc/genefilter/man/genefilter.html
# jaffelab::expression_cutoff() http://research.libd.org/jaffelab/reference/expression_cutoff.html
#
gene_means <- rowMeans(assay(rse_gene_SRP045638, "counts"))
summary(gene_means)

## Eliminamos genes
rse_gene_SRP045638 <- rse_gene_SRP045638[gene_means > 0.1, ]

## Dimensiones finales
dim(rse_gene_SRP045638)

## Porcentaje de genes que retuvimos
round(nrow(rse_gene_SRP045638) / nrow(rse_gene_SRP045638_unfiltered) * 100, 2)


library("edgeR") # BiocManager::install("edgeR", update = FALSE)
dge <- DGEList(
  counts = assay(rse_gene_SRP045638, "counts"),
  genes = rowData(rse_gene_SRP045638)
)
dge <- calcNormFactors(dge)

dge

# Revisar si assigned fene prop tiene diferencias entre los grupos
library("ggplot2")
ggplot(as.data.frame(colData(rse_gene_SRP045638)), aes(y = assigned_gene_prop, x = prenatal)) +
  geom_boxplot() +
  theme_bw(base_size = 20) +
  ylab("Assigned Gene Prop") +
  xlab("Age Group")

# Creamos el modelo
mod <- model.matrix(~ prenatal + sra_attribute.RIN + sra_attribute.sex + assigned_gene_prop,
                    data = colData(rse_gene_SRP045638)
)
colnames(mod)

# Analisis de expresión deiferencial
library("limma")
vGene <- voom(dge, mod, plot = TRUE)


#
eb_results <- eBayes(lmFit(vGene))

de_results <- topTable(
  eb_results,
  coef = 2,
  number = nrow(rse_gene_SRP045638),
  sort.by = "none"
)
dim(de_results)
head(de_results)

## Visualicemos los resultados estadísticos
plotMA(eb_results, coef = 2)


volcanoplot(eb_results, coef = 2, highlight = 3, names = de_results$gene_name)

de_results[de_results$gene_name %in% c("ZSCAN2", "VASH2", "KIAA0922"), ]

## Extraer valores de los genes de interés
exprs_heatmap <- vGene$E[rank(de_results$adj.P.Val) <= 50, ]
geneNames <- rowRanges(rse_gene_SRP045638)$gene_name[match(row.names(exprs_heatmap),rowRanges(rse_gene_SRP045638)$gene_id)]
row.names(exprs_heatmap) <- geneNames
## Creemos una tabla con información de las muestras
## y con nombres de columnas más amigables
df <- as.data.frame(colData(rse_gene_SRP045638)[, c("prenatal", "sra_attribute.RIN", "sra_attribute.sex")])

colnames(df) <- c("AgeGroup", "RIN", "Sex")

## Hagamos un heatmap
library("pheatmap")
pheatmap(
  exprs_heatmap,http://127.0.0.1:28855/graphics/8d322f21-57b1-4789-9e67-af5006fa20ab.png
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_col = df
)

#

## Para colores
library("RColorBrewer")

## Conviertiendo los grupos de edad a colores
col.group <- df$AgeGroup
levels(col.group) <- brewer.pal(nlevels(col.group), "Set1")

col.group <- as.character(col.group)

## MDS por grupos de edad
#plotMDS(vGene$E, labels = df$AgeGroup, col = col.group)
