# SEURAT ANALYSIS USING SOLOTE OUTPUT

After SoloTE has finished processing the BAM file of interest, the following results should be available:

```
EXPERIMENTNAME_SoloTE_output  
├── barcodes.tsv
├── features.tsv
└── matrix.mtx
```

Where "EXPERIMENTNAME" corresponds to the "OutputName" value used when running SoloTE. For the remaining of this tutorial, the files from the 2C like data available [here](https://github.com/bvaldebenitom/SoloTE/tree/main/Data_2Clike_SoloTE) will be used.


Then, processing with Seurat could be done as follows (in this example, the [Seurat - Guided Clustering Tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html) commands are used):

## Load libraries and files
```
library(Seurat)
solote_matrix <- ReadMtx("matrix.mtx.gz","barcodes.tsv.gz","features.tsv.gz",feature.column=1)
solote_seuratobj <- CreateSeuratObject(count=solote_matrix,min.cells=3,project="SoloTE")
```

## Pre-processing workflow
```
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
```



## Marker analysis
```
#Rename cells based on the expression of known marker genes for the 2C-like group ("Zscan4c", "Tcstv3")
seuratobj <- RenameIdents(seuratobj,'0' = "Non_2C",'1' = "Non_2C",'2' = "Non_2C",'3' = "Non_2C",'4' = "Non_2C",'5' = "Non_2C",'6' = "Non_2C",'7' = "Non_2C",'8' = "2C-like")
newmarkers <- FindAllMarkers(seuratobj_cnt,only.pos=TRUE)
newmarkers_signif <- newmarkers[which(newmarkers$p_val_adj<=0.05),]
```
Up to this point, we have a statistically significant set of markers (both genes and TEs). Then, to filter the TE results, we can take advantage of the SoloTE feature of adding the "SoloTE" keyword to TE identifiers:

```
newmarkers_signif_te <- newmarkers_signif[grep("SoloTE",newmarkers_signif$gene),]
head(newmarkers_signif_te)
```

For this example, the following results can be obtained:
```
                                                              p_val avg_log2FC pct.1 pct.2 p_val_adj cluster                                                          gene
SoloTE-chr18-80020461-80020681-Lx3C:L1:LINE-16.7-+-INTERGENIC     0   3.632743 0.798 0.005         0 2C-like SoloTE-chr18-80020461-80020681-Lx3C:L1:LINE-16.7-+-INTERGENIC
SoloTE-chr18-80017194-80017417-Lx3B:L1:LINE-16.7-+-INTERGENIC     0   2.685660 0.706 0.004         0 2C-like SoloTE-chr18-80017194-80017417-Lx3B:L1:LINE-16.7-+-INTERGENIC
SoloTE-chr18-80030263-80030464-Lx3B:L1:LINE-22.4-+-INTERGENIC     0   2.507013 0.697 0.006         0 2C-like SoloTE-chr18-80030263-80030464-Lx3B:L1:LINE-22.4-+-INTERGENIC
SoloTE-chr18-80013563-80013812-Lx3B:L1:LINE-18.4-+-INTERGENIC     0   2.474788 0.697 0.005         0 2C-like SoloTE-chr18-80013563-80013812-Lx3B:L1:LINE-18.4-+-INTERGENIC
SoloTE-chr16-9874989-9875412-L1-Mus1:L1:LINE-7.8-+-INTRONIC       0   2.337001 0.743 0.001         0 2C-like   SoloTE-chr16-9874989-9875412-L1-Mus1:L1:LINE-7.8-+-INTRONIC
SoloTE-chr18-80023718-80023941-Lx3B:L1:LINE-17.7-+-INTERGENIC     0   2.228682 0.670 0.002         0 2C-like SoloTE-chr18-80023718-80023941-Lx3B:L1:LINE-17.7-+-INTERGENIC
```







