#script to demostrate how to load single cell matrices in 
#several formats and converting to seurath object

getwd()

library("Seurat")
library(SeuratDisk)

# .RDS format
rds_obj <- readRDS("ependymal_cells.rds")
str(rds_obj)

# HD5 format
#son archivos 1 poco mas grande y tardan mas en cargar
hdf5_obj <- Read10X_h5(filename = "20k_PBMC_3p_HT_nextgem_Chromium_X_filtered_feature_bc_matrix.h5",
                     use.names = TRUE, #parametros predeterminados
                     unique.features = TRUE) #caracteristicas (genes) unicas)

#no es objeto seurat, asi que lo creamos
seurat_hd5 = CreateSeuratObject(counts = hdf5_obj)

# .mtx file
mtx_obj <- ReadMtx(mtx = "raw_feature_bc_matrix/matrix.mtx.gz",
        features = "raw_feature_bc_matrix/features.tsv.gz",
        cells = "raw_feature_bc_matrix/barcodes.tsv.gz")
seurat_mtx = CreateSeuratObject(counts = mtx_obj)

# .h5ad format (ann data format, es para scanpy)

#paso 1, convertir anndata a h5seurat
Convert("adata_SS2_for_download.h5ad", dest = "h5seurat", overwrite = TRUE)

#paso 2: cargar h5seurat como objeto seurat (paso anterior crea archivo nuevo con el mismo nombre pero extension h5seurat)
seurat_anndata <- LoadH5Seurat("adata_SS2_for_download.h5seurat")
