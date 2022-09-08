# SEURAT ANALYSIS USING SOLOTE OUTPUT

After SoloTE has finished processing the BAM file of interest, the following results should be available:

`
EXPERIMENTNAME_SoloTE_output
├── barcodes.tsv
├── features.tsv
└── matrix.mtx
`

Where "EXPERIMENTNAME" corresponds to the "OutputName" value used when running SoloTE. 


Then, processing with Seurat could be done as follows (in this example, the [Seurat - Guided Clustering Tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html) commands are used):

##Load libraries and files<br/>
`
library(Seurat)
library(ggplot2)
matrix <- "matrix.mtx"
barcodes <- "barcodes.tsv"
features <- "features.tsv"
solote_matrix <- ReadMtx(matrix,barcodes,features,feature.column=1)
solote_seuratobj <- CreateSeuratObject(count=solote_matrix,min.cells=3,project="SoloTE")
`

##Pre-processing workflow
`
solote_seuratobj$percent_mt <- PercentageFeatureSet(solote_seuratobj,pattern="^mt-")
seuratobj <- subset(solote_seuratobj,subset = nFeature_RNA>=200 & percent_mt<=5)
seuratobj <- NormalizeData(seuratobj)
seuratobj <- FindVariableFeatures(seuratobj, selection.method = "vst", nfeatures = 1000)
all.genes <- rownames(seuratobj)
seuratobj <- ScaleData(seuratobj, features = all.genes)
seuratobj <- RunPCA(seuratobj, features = VariableFeatures(object = seuratobj))
seuratobj <- FindNeighbors(seuratobj, dims = 1:10)
seuratobj <- FindClusters(seuratobj, resolution = 0.5)
seuratobj <- RunUMAP(seuratobj, dims = 1:10)
DimPlot(seuratobj,reduction="umap")
`



##Clustering
`
code
`

##Marker analysis
`
code
`






