#script to perform standar workflow stepts to analyze single cell rnaseq
getwd()

library("Seurat")
library("tidyverse")

nsclc.sparse.m <- Read10X_h5(filename = "20k_NSCLC_DTC_3p_nextgem_Multiplex_count_raw_feature_bc_matrix.h5",
                     use.names = TRUE, #parametros predeterminados
                     unique.features = TRUE) #caracteristicas (genes) unicas)

#da un warning que tiene multiples modalidades (genome matrix has multiple modalities...)
#vemos las modalidades
str(nsclc.sparse.m)
#vemos List of 3  Gene Expression, Antibody Capture y Multiplexing Capture.
#estamos interesados en gene expression

cts <- nsclc.sparse.m$'Gene Expression'

#iniciamos el objeto seurat
#min cell queremos mantener todas las caracteristicas que tengan expresion en al menos 3 celulas, y mantenemos todas las celulas que tengan al menos 200 genes
nsclc.surat.obj = CreateSeuratObject(counts = cts, project = "NSCLC", min.cells = 3, min.features = 200)

#2. control de calidad -----------------------
# porcentaje de genes mitocondriales (contaminacion)
View(nsclc.surat.obj@meta.data)
nsclc.surat.obj[["percent.mt"]] <- PercentageFeatureSet(nsclc.surat.obj, pattern = "^MT-")
View(nsclc.surat.obj@meta.data)

#grafico de violin en seurat
VlnPlot(nsclc.surat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#inspeccion visual de la calidad

#nFeature_RNA celulas con numero diferente de genes
#nCount_RNA celulas con numeor de moleculas detectadas

FeatureScatter(nsclc.surat.obj, feature1="nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method="lm")

#filtrando celulas con baja calidad

nsclc.surat.obj <- subset(nsclc.surat.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
nsclc.surat.obj
#hay gente que tambien filtra x ribosomal

#3. Normalizando datos --------
#nsclc.surat.obj <- NormalizeData(nsclc.surat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
# esos son valores por defecto, puedo hacer solo
nsclc.surat.obj <- NormalizeData(nsclc.surat.obj)

str(nsclc.surat.obj)

#4. identificar genes/caracteristicas altamente variables ------
nsclc.surat.obj <- FindVariableFeatures(nsclc.surat.obj, selection.method = "vst", nfeatures= 2000)

#identificando los genes mas variables
top10 <- head(VariableFeatures(nsclc.surat.obj), 10)

#graficando caracteristicas variables con y sin etiquetas
#en rojo se veran las caracteristicas identificadas mas variables encontradas con
#variable features
plot1 <- VariableFeaturePlot(nsclc.surat.obj) # Genera el gráfico de genes variables
plot1 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) # Agrega etiquetas
plot1

# 5. escalamiento ------
#en single cell hay varias fuentes de variacion no asociado a evento biologico
#(ej ruido tecnico o fuente biologica como diferencias en ciclo celular) 
#ruido tecnico es batch effect, asi que queremos contar esas variaciones
#para los analisis posteriores y para eso usaremos scaledata, usando todos los genes
#como caracteristicas, porque despues usaremos tecnicas de reduccion de dimensionalidad
#y quiero que este todo escalado

all.genes <- rownames(nsclc.surat.obj)
nsclc.surat.obj <- ScaleData(nsclc.surat.obj, features = all.genes)

#si hago str(nsclc.surat.obj) en RNA formal class assay encuentro 3 slots 
# counts que contiene la raw data, data que tiene datos log normalizados y ahora scale.data los datos escalados

#6. reduccion de dimencionalidad ----
#por defecto considera las caracteristicas variables en el objeto seurat, pero si queremos podemos darle todos los gene
#que estamos interesados
nsclc.surat.obj <- RunPCA(nsclc.surat.obj, features = VariableFeatures(object = nsclc.surat.obj))

#el paso anterior nos muestra los componentes principales y los genes con score positivos y negativos del componente
#asi que ahora veremos los 5 primeros por dimensiones que nos dio 
print(nsclc.surat.obj[["pca"]], dims = 1:5, nfeatures = 5)

#heatmap para identificar score de los componentes (heatmap escalado por PC score)
# podemos ver si efectivamente captura heterogeneidad
DimHeatmap(nsclc.surat.obj, dims=1, cells=100, balanced = T) #dims = que pc quiero ver

#determinamos dimensionalidad de los datos, para usar solo los PC estadisticamente significativos que capturan dimensionalidad de la señal para analisis posteriores
#uno de las formas es un grafico de codo, asi cada PC se ranquea por desviacion
ElbowPlot(nsclc.surat.obj)

#entre 10 y 15 se explica mayor varianza, asi que el corte se hara en el PC 15

# 7. clusterizar -------
#agruparemos celulas similares (expresion similar), para eso encontraremos vecinos
#para eso le damos el objeto y los PC con mayor varianza en el dataset
nsclc.surat.obj <- FindNeighbors(nsclc.surat.obj, dims= 1:15)

#en el siguiente paso queremos asignar celulas a clusters y usamos un parametro de resolucion para determinar la granularidad, a mayor resolucion mayor numero de clusters
nsclc.surat.obj <- FindClusters(nsclc.surat.obj, resolution = c(0.1, 0.3, 0.5, 0.7, 1))
#vemos metadata, y nos genera columnas con el numero de grupos por resolucion (vemos pertenencia a un grupo)
View(nsclc.surat.obj@meta.data)
#podemos graficar cada grupo por resolucion (ej 0.1 nos da 8 grupos, 0.5 nos da 12)
DimPlot(nsclc.surat.obj, group.by = "RNA_snn_res.0.3", label=T)

#configurando identidad por clusters (que celula pertenece a cada cluster)
#por defecto nos muestra el numero de clusters a resulocion por defecto
Idents(nsclc.surat.obj)

#puedo cambiar la resolucion por defecto
Idents(nsclc.surat.obj) <- "RNA_snn_res.0.1"
Idents(nsclc.surat.obj)

# reduccion de dimensionalidad no lineal para agrupar celulas en un espacio de baja dimensionalidad
# y podemos explorar y visualizar los datos
nsclc.surat.obj <- RunUMAP(nsclc.surat.obj, dims = 1:15)
#notar que podemos usar label=TRUE o usar LabelClusters para ayudar a anotar clusters individuales
DimPlot(nsclc.surat.obj, reduction = "umap", label=T)
