# Regresiones lineales

# Una regresión lineal esun modelo en el que una variable
# de resppuesta Y es el resultado de una variable independiente
# reescalada por un factor B1 + un valor B0. B1 es la esperanza
# de Y o cambio promedio de Y

## ?model.matrix
mat <- with(trees, model.matrix(log(Volume) ~ log(Height) + log(Girth)))
mat

# Intercept is B0: The value of Y when X == 0

summary(lm(log(Volume) ~ log(Height) + log(Girth), data = trees))

# We use logarithmic scale to make this linear!!!!

# Un análisis de expresión diferencial es, en esencia, un conjunto de regresiones lineales

# t = Estimado / sd(Estimado)

## Datos de ejemplo
(sampleData <- data.frame(
  genotype = rep(c("A", "B"), each = 4),
  treatment = rep(c("ctrl", "trt"), 4)
))

## Creemos las imágenes usando ExploreModelMatrix
vd <- ExploreModelMatrix::VisualizeDesign(
  sampleData = sampleData,
  designFormula = ~ genotype + treatment,
  textSizeFitted = 4
)

## Veamos las imágenes
cowplot::plot_grid(plotlist = vd$plotlist)

# Esta tabla nos ayuda a interpretar los valores de los coeficientes
# del modelo linear


# En el output de modelmatrix tenemos columnas indicator con valores
# binarios (0 = referencia, 1 = no referencia), para las variables categóricas


(sampleData <- data.frame(
  Response = rep(c("Resistant", "Sensitive"), c(12, 18)),
  Patient = factor(rep(c(1:6, 8, 11:18), each = 2)),
  Treatment = factor(rep(c("pre","post"), 15)),
  ind.n = factor(rep(c(1:6, 2, 5:12), each = 2))))

vd <- ExploreModelMatrix::VisualizeDesign(
  sampleData = sampleData,
  designFormula = ~ Response + Response:ind.n + Response:Treatment,
  textSizeFitted = 3
)
cowplot::plot_grid(plotlist = vd$plotlist, ncol = 1)

# o + ... es una sintaxis de model matrix para
# indicar que no queremos intersection

# Ejercicio
(sampleData <- data.frame(
  Response = rep(c("Resistant", "Sensitive"), c(12, 18)),
  Patient = factor(rep(c(1:6, 8, 11:18), each = 2)),
  Treatment = factor(rep(c("pre","post"), 15)),
  ind.n = factor(rep(c(1:6, 2, 5:12), each = 2))))
vd <- ExploreModelMatrix::VisualizeDesign(
  sampleData = sampleData,
  designFormula = ~ Response + Response:ind.n + Response:Treatment,
  textSizeFitted = 3
)
cowplot::plot_grid(plotlist = vd$plotlist, ncol = 1)
# Ver imagen :

# ¿Por qué es clave el 0 al inicio de la fórmula en el ejercicio 3?
# Porque es la sintaxis para especificar que se elimina la columna de intersection
