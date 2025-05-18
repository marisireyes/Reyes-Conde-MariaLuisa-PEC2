### PEC 2 - ANALISIS DE DATOS OMICOS - Msria Luisa Reyes Conde
## APARTADO 1
# librerias q nos serviran
library(GEOquery)
library(EnsDb.Hsapiens.v86)
library(SummarizedExperiment)
# Descarga la matriz de expresión correspondiente y los metadatos, y cárgalos en R. 
getGEOSuppFiles("GSE161731", makeDirectory = TRUE) #para extraer los archivos de GEO
# cargamos la matriz de expresiones
matriz_expr <- as.matrix(read.csv("GSE161731/GSE161731_counts.csv", row.names = 1, check.names = FALSE))
mode(matriz_expr) <- "numeric" 
# metadatos
metadatos <- read.csv("GSE161731/GSE161731_key.csv", stringsAsFactors = FALSE)
rownames(metadatos)<-metadatos$rna_id

# Construye un objeto SummarizedExperiment que contenga ambos. 
# Agrega también las coordenadas génicas como rowRanges. Ten en cuenta que necesitas tener muestras en común entre los metadatos 
# y la matriz de expresión antes de crear el objeto SummarizedExperiment, así como genes en común entre la matriz de expresión y 
# la anotación (por ejemplo, EnsDb.Hsapiens.v86).

# como tenemos que asegurarnos de que las muestras coincidan entre los metadatos y la expresion hacemos lo siguiente; buscamos las muestras
# que están tanto en la matriz de expresión como en los metadatos y "filtramos" estos dos conjuntos de datos (nos quedamos con las muestras
# comunes de la matriz y de los metadatos)
muestras_comunes <- intersect(colnames(matriz_expr), rownames(metadatos))
matriz_expr <- matriz_expr[, muestras_comunes]
metadatos <- metadatos[muestras_comunes, ]

# para obtener las coordenadas genicas de estas muestras usamos, como nos indica, la funcion "genes(EnsDb.Hsapiens.v86)"
coord_gen <- genes(EnsDb.Hsapiens.v86)
# para asegurarnos de que los genes de nuestra matriz de expresion estan incluidos en "coord_gen" vovlemos a intersecar como antes y asi
# nos quedamos con los genes que sean comunes:
genes_comunes <- intersect(rownames(matriz_expr), names(coord_gen))
matriz_expr <- matriz_expr[genes_comunes, ]
coord_gen <- coord_gen[genes_comunes]


# creamos ahora el summarizedExperiment
se <- SummarizedExperiment(
  assays = list(counts = matriz_expr),
  colData = metadatos,
  rowRanges = coord_gen
)
se

## APARTADO 2
# Limpia los metadatos y selecciona solamente tres cohortes: COVID19, Bacterial y healthy.
cohortes <- c("COVID-19", "Bacterial", "healthy")
metadatos <- metadatos[metadatos$cohort %in% cohortes, ]

#  Sugerencias: (i) elimina individuos duplicados (conserva solo la primera entrada cuando haya duplicados); 
metadatos <- metadatos[!duplicated(rownames(metadatos)), ]

# (ii) asegúrate de que las variables sean del tipo adecuado (por ejemplo, la edad debe ser numérica, etc.); 
metadatos$age <- as.numeric(metadatos$age)
metadatos$batch <- as.factor(metadatos$batch) #solo toma vslores 1 y 2


# Importante: una vez hecho esto, selecciona exclusivamente 75 muestras de manera aleatoria, utilizando la semilla que se
# proporciona a continuación (basada en tu primer nombre, primer y segundo apellidos sin tildes ni espacios y en minúscula).
myseed <- sum(utf8ToInt("marialuisareyesconde")) 
set.seed(myseed)
# seleccionamos las muestras ocn la semilla seleccionada y filtramos los metadatos
muestras_selec <- sample(rownames(metadatos), 75)
metadatos <- metadatos[muestras_selec, ]

# filtramos tambien la matriz y el se
matriz_expr <- matriz_expr[, muestras_selec]
se <- se[, muestras_selec]

# eliminamos guiones, barras y espacios

metadatos[] <- lapply(metadatos, function(x) {
  if (is.character(x) || is.factor(x)) {
    x <- as.character(x)           # convertir factor a carácter
    gsub("[- /]", "", x)            # eliminar guiones y espacios
  } else {
    x
  }
})

## APARTADO 3
library(DESeq2)
# creamos un objeto DESeqDataSet para poder hacer el analisis de expresion diferencial
dds <- DESeqDataSet(se, design = ~ cohort)
# Lleva a cabo el preprocesado inicial de los datos (eliminación de genes con baja expresión, etc.) 
genes_buenos <- rowSums(counts(dds) >= 10) >= (0.1 * ncol(dds))
dds <- dds[genes_buenos, ]
# y la transformación/normalización que consideres apropiada.
vst_datos <- varianceStabilizingTransformation(dds, blind = TRUE)
vst_matriz <- assay(vst_datos) #matriz normalizada

## APARTADO 4
library(ggplot2)
library(pheatmap)
library(limma)
# Realiza un análisis exploratorio sobre los datos transformados/normalizados.
metadatos <- as.data.frame(colData(se)) # metadatos transformados/normalizados
#Puedes usar PCA, MDS, clustering y/o heatmaps. 
# PCA:
pca <- prcomp(t(vst_matriz), scale. = TRUE)
# agregamos los pc a los metadatos para graficarlo
metadatos$PC1 <- pca$x[,1]
metadatos$PC2 <- pca$x[,2]
# grafico coloreado x cohorte
p_cohort<-ggplot(metadatos, aes(x = PC1, y = PC2, color = cohort)) +
  geom_point(size = 3) +
  theme_minimal() 
p_cohort
# MDS
mds <- plotMDS(vst_matriz, plot = FALSE)
metadatos$MDS1 <- mds$x
metadatos$MDS2 <- mds$y
# grafico coloreado x cohorte
ggplot(metadatos, aes(x = MDS1, y = MDS2, color = cohort)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "MDS - coloreado por cohorte")


# Identifica y elimina muestras atípicas (outliers). 
# Ejemplo: eliminar muestras fuera de +/- 3 SDs en PC1 o PC2
outlier_idx <- which(abs(metadatos$PC1) > 3 * sd(metadatos$PC1) |
                       abs(metadatos$PC2) > 3 * sd(metadatos$PC2))

# Nombres de muestras atípicas
outliers <- rownames(metadatos)[outlier_idx]

# Eliminar outliers
if (length(outliers) > 0) {
  se <- se[, !colnames(se) %in% outliers]
  metadatos <- as.data.frame(colData(se))
  vst_matriz <- assays(se)$vst
}

# De acuerdo con los resultados de estos análisis, identifica qué variables 
# en los metadatos pueden considerarse variables confusoras y por tanto deberían incluirse en la matriz de diseño. 

# grafico coloreado x cohorte + genero
p_cohort + 
  geom_point(aes(shape = gender), color = "black", size = 2.5) +
  scale_shape_manual(values = c("Female" = 1, "Male" = 4)) +
  labs(title = "PCA - coloreado por cohorte y genero")
  theme_minimal()
  # grafico coloreado x cohorte + raza
  p_cohort + 
    geom_point(aes(shape = race), color = "black", size = 2.5) +
    #scale_shape_manual(values = c("Female" = 1, "Male" = 4)) +
    labs(title = "PCA - coloreado por cohorte y raza")
  theme_minimal()
  # grafico coloreado x cohorte + hospitalizado
  p_cohort + 
    geom_point(aes(shape = hospitalized), color = "black", size = 2.5) +
    scale_shape_manual(values = c("Yes" = 1, "No" = 4)) +
    labs(title = "PCA - coloreado por cohorte y raza")
  theme_minimal()

# Añadiria genero y raza 
  
  
## APARTADO 5
# Construye la matriz de diseño y las matrices de contrastes adecuadas para evaluar la expresión génica diferencial 
# en las comparaciones Bacterial vs healthy y COVID19 vs healthy y realiza un análisis de expresión diferencial.
# semilla y determinar metodo
set.seed(myseed)
metodos <- sample(c("edgeR", "voom+limma", "DESeq2"), size = 1)
metodos # = voom+limma
library(edgeR)
dge <- DGEList(counts = assay(se, "counts"))
# Crear matriz de diseño
matriz_dis <- model.matrix(~ 0 + cohort + race + gender, data = metadatos)  # Añade aquí tus variables confusoras
# Aplicar voom
v <- voom(dge, matriz_dis, plot = TRUE)
metadatos[] <- lapply(metadatos, function(x) {
  if (is.character(x) || is.factor(x)) {
    x <- as.character(x)           # convertir factor a carácter
    gsub("[- /]", "", x)            # eliminar guiones y espacios
  } else {
    x
  }
})
colnames(matriz_dis) <- gsub("[- /]", "", colnames(matriz_dis))
# hacemos los contrastes para comparar:
matriz_contrs <- makeContrasts(
  Bacterial_vs_healthy = cohortBacterial - cohorthealthy,
  COVID19_vs_healthy = cohortCOVID19 - cohorthealthy,
  levels = matriz_dis
)

# por ultimo ajustamos con limma
ajuste <- lmFit(v, matriz_dis)
ajuste2 <- contrasts.fit(ajuste, matriz_contrs)
ajuste2 <- eBayes(ajuste2)

# resultados (logFC > 1.5 como umbral)
res_bacterial <- topTable(ajuste2, coef = "Bacterial_vs_healthy", number = Inf, adjust = "fdr", lfc = log2(1.5))
res_covid <- topTable(ajuste2, coef = "COVID19_vs_healthy", number = Inf, adjust = "fdr", lfc = log2(1.5))

# filtramos tb solo significativos
sig_bacterial <- subset(res_bacterial, adj.P.Val < 0.05)
sig_covid <- subset(res_covid, adj.P.Val < 0.05)



## APARTADO 6
# Compara los resultados de ambos contrastes (Bacterial vs healthy y COVID19 vs healthy).
# vamos a usar un diagrama de VEnn 
library(VennDiagram)

# Obtenemos los genes diferencisdos d cada uno
genes_bact <- rownames(sig_bacterial)
genes_covid <- rownames(sig_covid)

# Diagrama de Venn
grid.newpage()
diag_venn <- venn.diagram(
  list(Bacterial = genes_bact, COVID19 = genes_covid),
  filename = NULL,
  fill = c("red", "blue"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  main = "Genes diferencialmente expresados"
)

grid::grid.draw(diag_venn)

## APARTADO 7
# Realiza un análisis de sobrerepresentación para identificar las funciones enriquecidas entre los
# genes sobreexpresados en pacientes con COVID19 en comparación con los controles sanos. 
# Utiliza solo el dominio de Gene Ontology “Biological Process”. 
library(clusterProfiler) # realiza el analisis go
library(org.Hs.eg.db) # para mapear los id de genes a terminos biologicos 

# Convertir Ensembl a Entrez ID
genes_up_covid <- rownames(subset(sig_covid, logFC > 0))  # cgemos solo los genes sobreexpresados
entrez_ids <- mapIds(org.Hs.eg.db,
                     keys = genes_up_covid,
                     keytype = "ENSEMBL",
                     column = "ENTREZID",
                     multiVals = "first")
entrez_ids <- na.omit(entrez_ids)

# GO Biological Process:  test de hipergeométrica para evaluar si los genes de nuestra lista 
# están sobre-representados en ciertas categorías GO
ego_bp <- enrichGO(
  gene = entrez_ids,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  keyType = "ENTREZID",
  pAdjustMethod = "fdr",
  qvalueCutoff = 0.05,
  readable = TRUE
)

# visualizacion (sin revigo)
barplot(ego_bp, showCategory = 20, title = "GO: Biological Process (COVID19)")
dotplot(ego_bp, showCategory = 20)

