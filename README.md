scRNASeq Analysis Pipeline
================

# Table of Contents

1.  [Introduction](#introduction)  
2.  [Recommended Readings/Resources](#recommended-readingsresources)  
3.  [Setup](#setup)
    - [CellRanger](#cellranger)  
    - [Conda Environment](#conda-environment)  
    - [R Packages](#r-packages)  
4.  [Ambient RNA Contamination
    Correction](#ambient-rna-contamination-correction)  
5.  [Importing Samples](#importing-samples)  
6.  [QC and Doublet Detection](#qc-and-doublet-detection)  
7.  [Integration](#integration)  
8.  [Clustering](#clustering)  
9.  [Annotation](#annotation)  
10. [Plots](#plots)

# Introduction

This is code for standard scRNAseq analysis. It includes pre-processing
of raw data, QC, Normalization, Integration, Clustering and Differential
expression analysis.

The code is heavily inspired by Seurat pipeline but includes additional
tweaks for Pancreas tissue analysis

# Recommended Readings/Resources

- [Orchestrating Single-Cell Analysis with
  Bioconductor](https://bioconductor.org/books/release/OSCA/) : One of
  the core reading material to understand single-cell analysis workflow
- [Seurat](https://satijalab.org/seurat/) : Seurat is a popular R
  package for single-cell RNA-seq data analysis. Their website have
  extensive documentation for the package functionalities
- [Analysis of single cell RNA-seq
  data](https://www.singlecellcourse.org/index.html)
- [scVerse](https://scverse.org/): scVerse encompases different tools
  and methods for scRNAseq data analysis in python
- [Benchmarking atlas-level data integration in single-cell
  genomics](https://www.nature.com/articles/s41592-021-01336-8)
- [Considerations for building and using integrated single-cell
  atlases](https://www.nature.com/articles/s41592-024-02532-y)

# Setup

For this workflow we will need to have certain softwares installed (R,
Python, CellRanger). We will need to install several R packges. We will
also create conda environemnt that contains essential python packages
used in the workflow

## CellRanger

CellRanger can be downloaded
[here](https://www.10xgenomics.com/support/software/cell-ranger/latest)
and also read the documentation for its usage. Many academic HPC servers
have it already installed so you might only need to load the module.

Essential considerations whenever using CellRanger are the chemistry,
the reference genome and the software version. 10X keep releasing new
chemistries that enhance the sequencing workflow, and they keep
releasing new versions of the software that might have parse the
sequencing data differently. For example, starting cellranger7.0 they
started making use of the detected introns in counting the UMIs for each
gene, which indeed changes the overall output for a given sample if
analyzed using cellranger6.x vs cellranger7.x. Lastly, as new versions
of the reference genomes gets released (different versions of hg19,
hg38), the genes are annotated differently.

Generally make sure to use the same genome reference version and the
software version (and preferably the same chemistry) if you are
analyzing cohort of samples to be integrated in the same analysis

You can review the code for running CellRanger on one or multiple
samples [here]()

## Conda Environment

We will need to use CellBender for ambient RNA correction. We might want
to use scVI for integration. We will also use Scrublet for doublet
detection. Here I will use these methods both as Slurm script and inside
R scripts using reticulate for better flow. But you can use them
independently in jupyter notebook or python script. You can create a
conda environment that has the following commands in bash

``` bash
conda create --name single_cell_env python=3.8 
## or you can use `conda create --prefix /path/to/env python=3.8` if you would like to install it in certain directory rather than your home directory
conda activate single_cell_env
## if installed in specific directory use `conda activate /path/to/env`
pip install cellbender
pip install scrublet
pip install scvi-tools
pip install scanpy
pip install scanpy
```

## R Packages

``` r
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratWrappers)
  library(tidyverse)
  library(scCustomize)
  library(CellChat)
  library(slurmR)
  library(CoGAPS)
  library(reticulate)
  library(qs)
})

use_condaenv("/path/to/conda/envs/single_cell/", required = TRUE)
sc <- import("scanpy", convert = F)
anndata <- import("anndata", convert = F)
scipy <- import('scipy', convert = F)
```

# Ambient RNA contamination correction

Many tools have been proposed for correcting ambient RNA contamination.
While [soupX](https://github.com/constantAmateur/SoupX) is one of the
best methods for ambient RNA correction, we might want to use
[CellBender](https://cellbender.readthedocs.io/en/latest/#) as it is
completely unsupervised and does not need to setup many parameters
before running it.

You can run cellbender using bash thru the following line

``` bash
cellbender remove-background \
--cuda \
--input /path/to/cellranger/output/raw_feature_bc_matrix.h5 \
--output /path/to/outputDir/sample_id.h5
```

or you can use this slurm script if you have access to HPC and would
like to parallelize the processing of the samples

# Importing samples

We often analyze multiple sample together. To make our life easier, we
should create an spreadsheet that at least has two columns: sample id
and cellranger output path. One can add as much metadata to each sample
as needed as separate columns (e.g., disease condition, run ID, patient
ID, …). Then we use this spreadsheet as our base for importing samples
and merging them.

Note: If you run CellBender, you will need to read using
`Read_CellBender_h5_Mat()` from
[scCustomize](https://samuel-marsh.github.io/scCustomize/) R package
instead of using `Read10X()` from Seurat. We will not need the
CellRanger output path here, instead we’ll need to path for CellBender
output.

``` r
# Importing samples manifest spreadsheet ----------------------------------------

samples_info <- readxl::read_xlsx("samples_manifest.xlsx")
samples <- samples_info$sample_id

# Importing samples --------------------------------------------------------------

## Reading h5 files and adding metadata to each one based on sample_info
scData <- lapply(samples, function(x){
  sample_info <- samples_info %>% filter(sample_id == x)
  seurat <- Read_CellBender_h5_Mat(paste0("outputs/cellbender/scRNAseq/",x,"/", x, "_filtered.h5")) %>% CreateSeuratObject()
  seurat$sample_id <- x
  seurat$patient_id <- sample_info$patient_id
  seurat$tissue <- sample_info$tissue
  seurat <- RenameCells(seurat, add.cell.id = x)
  return(seurat)
})

## merging the samples into one object
scData <- reduce(scData, merge)
scData <- JoinLayers(scData, assay = 'RNA')
scData[["RNA"]] <- split(scData[["RNA"]], f = scData$sample_id)
```

# QC and Doublet Detection

``` r
## Percent mitochondrial
scData[["percent.mt"]] <- PercentageFeatureSet(object = scData, pattern = "^MT-")

## Running Scrublet
adata <- anndata$AnnData(
  X = scipy$sparse$csr_matrix(
    Matrix::t(LayerData(scData, assay = 'RNA', layer = "counts"))
  ),
  obs = scData@meta.data
)
adata$var_names = rownames(scData)
adata$obs_names = colnames(scData)
sc$pp$scrublet(adata, batch_key = 'sample_id')
scrublet_res <- py_to_r(adata$obs) %>% select(doublet_score, predicted_doublet)
scData <- AddMetaData(scData, scrublet_res)
```

If using the batch_key argument doesn't work, we can manually apply this on sample-by-sample basis using the following snippet

```
scrublet_res <- lapply(unique(scData$sample_id), function(x){
  message(x)
  sample <- subset(scData, subset = sample_id == x)
  adata <- anndata$AnnData(
    X = scipy$sparse$csr_matrix(
      Matrix::t(LayerData(sample, assay = 'RNA', layer = "counts"))
    )
  )
  sc$pp$scrublet(adata)
  scrublet_res <- py_to_r(adata$obs) %>% 
    select(doublet_score, predicted_doublet)
  rownames(scrublet_res) <- colnames(sample)
  scrublet_res$predicted_doublet <- unlist(scrublet_res$predicted_doublet)
  return(scrublet_res)
}) %>% bind_rows()
scData <- AddMetaData(scData, scrublet_res)
```

## If the previous code doesn’t work you can save the adata object and run scrublet using the python code in `scrublet.py`

``` r
adata.write("outputs/scRNASeq_integration/GoL_ref.h5ad")
scrublet_res <- read.csv("outputs/scRNASeq_integration/scrublet_analysis.csv")
scData <- AddMetaData(scData, scrublet_res)
```

We will remove the low quality cell based percent mitrochrondrial and
number of genes. We will not filter based on doublet score now, we will
wait until we cluster and see if a cluster has high doublet score and
has known markers of two distinct cell types.

``` r
scData <- subset(scData, subset = percent.mt < 15 & nFeature_RNA > 200)
```

# Integration

Integration is an extensive area of research in computational biology to
remove the batch effect. Seurat provided multiple ways to perform
integration as discussed
[here](https://satijalab.org/seurat/articles/seurat5_integration). Based
on my experience CCA Integration works best when you have samples of the
same tissue/condition/cellTypes but it takes long time. rPCA Integration
works best if you have samples that have common cell types but still
certain samples have distinct cell types (e.g. integrating tumor and
normal samples) and it takes shorter time compared to CCA ingtegration.
scVI integration is based on scVI model, it can scales to million of
cells as it can utilize GPU.

Note if you are using

``` r
scData <- FindVariableFeatures(scData)
scData <- ScaleData(scData)
scData <- RunPCA(scData)
scData <- IntegrateLayers(
  object = scData, method = rPCAIntegration,
  orig.reduction = 'pca',
  new.reduction = "integrated.rPCA", verbose = TRUE
)
# scData <- IntegrateLayers(
#   object = scData, method = scVIIntegration,
#   conda_env = "/path/to/conda/envs/single_cell_env/",
#   new.reduction = "integrated.scvi", verbose = TRUE
# )
scData <- RunUMAP(scData, reduction = 'integrated.rPCA', dims = 1:30)
scData <- FindNeighbors(scData, reduction = 'integrated.rPCA', dims = 1:30)
```

# Clustering

Here we will run unsupervised clustering based on the integrated
embedding (in our case rPCA.integrated) to define different clusters. We
will run DE to define markers of each cluster and start annotating those
clusters.

There are several ways to judge how many cluster would you like to
define (or what is the clustering resolution). One way is to visualize
the relation between clusters at different resolutions using
[clustree](https://github.com/lazappi/clustree) and one can judge what
is the optimal resolution that should be used and what are the stable
clusters across different resolutions

``` r
scData <- FindNeighbors(scData, reduction = 'integrated.rPCA', dims = 1:30)
res <- seq(from = 0.05, to = 1, by = 0.05)
for(x in res){
  scData <- FindClusters(scData, resolution = x, cluster.name = paste0("res.",x))
}
clustree(scData, prefix = 'res.', layout = "sugiyama")

res_of_choice <- 0.4
scData <- FindClusters(scData, resolution = res_of_choice)
scData <- JoinLayers(scData)
markers <- FindAllMarkers(scData)
scData@misc$markers <- markers
scData@misc$markers$score <- scData@misc$markers$avg_log2FC * scData@misc$markers$pct.1
```

We can visualize the UMAP embedding of the integrated data. Also we can
visualize the expression of known markers

``` r
features <- c("SPINK1", "CTRB1", "PRSS1", "AMY2A", "REG3A",
              "EPCAM","KRT19", "SOX9",
              "PDGFRA","DCN","LUM",
              "PDGFRB","RGS5",
              "PTPRC", "CD14", "APOE", "C1QA","CD68",
              "CD3E","CD8A", "CD4",
              "CD19","CD79A",'MS4A1',"JCHAIN",
              "FCGR3B", "S100A8",
              "TPSAB1","CPA3",
              "CDH5", "VWF",
              "INS", "CHGA", "GCG", "SST", "TAC1", 
              "NTS", 'PKHD1L1',
              'NCAM1', 'NEGR1', 'NRN1')
DimPlot(scData, group.by = 'sample_id') + DimPlot(scData)
DotPlot(scData, features = features) +
  scale_color_gradient2(low = "blue", mid = "yellow", high = "red", midpoint = 0) +
  theme_minimal() +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust =1, face = 'bold'),
        axis.text.y = element_text(face = 'bold', size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
```

# Annotation

Annotation here refers to the process of changing the names of the
clusters (1,2,3,…) to meaningful names (e.g. Acinar cells, Ductal cells,
Fibroblasts,…). This can be done manually by looking at the markers of
each cluster and comparing them to known markers of cell types. Although
one can create a named vector in R that has the correspondce between the
cluster number and the new name, I prefer to create an spreadsheet with
the cluster number in first column and annotaion in the second column.

``` r
annotations <- readxl::read_xlsx("outputs/scRNASeq_integration/scData_cluster_annotation.xlsx")
annotations <- setNames(annotations$annotation, annotations$cluster)
scData <- RenameIdents(scData, annotations)
scData$cellType <- Idents(scData)
scData@misc$markers <- merge(scData@misc$markers, data.frame(annotations) %>% rownames_to_column('cluster'), by = 'cluster')
```

At this point we save the integrated data as a R object instead of
re-running the whole analysis every time. We can save it as `RDS` object
using `saveRDS()`. However, it might be prefered to save it as `qs`
object using `qsave` from the [`qs`](https://github.com/qsbase/qs) R
package for faster read/write speed

``` r
qsave(scData, "outputs/scRNASeq_integration/scData.qs")
```

# Plots

After the annotation, you can visualize the same umap plots as above
using the new cluster names. You can also visualize the known markers
for each cluster. Lastly, you can visualize the cell type proportion in
each sample or each condition (tumor/normal)

``` r
cluster_cols <- CellChat::scPalette(length(levels(scData)))
names(cluster_cols) <- levels(scData)
DimPlot(scData, cols = cluster_cols, raster = F) 
```

``` r
Idents(scData) <- factor(Idents(scData),
                          levels = c("Acinar", "Acinar/Ductal", "Ductal", "Fibroblasts", "Pericytes",
                                     "Macrophages", "T_Cells", "B_Cells", "Granulocytes",
                                     "Mast", "Endothelial", "Endocrine", 
                                     "NTS_pos", "Neural"))
features <- c("CHRM3","CUZD1", "MECOM",
              "SPINK1", "CTRB1", "PRSS1", "AMY2A", "REG3A",
              "EPCAM","KRT19", "SOX9",
              "PDGFRA","DCN","LUM",
              "PDGFRB","RGS5",
              "PTPRC", "CD14", "APOE", "C1QA","CD68",
              "CD3E","CD8A", "CD4",
              "CD19","CD79A",'MS4A1',"JCHAIN",
              "FCGR3B", "S100A8",
              "TPSAB1","CPA3",
              "CDH5", "VWF",
              "INS", "CHGA", "GCG", "SST", "TAC1", 
              "NTS", 'PKHD1L1',
              'NCAM1', 'NEGR1', 'NRN1')

DotPlot(scData, features = features) +
  scale_color_gradient2(low = "blue", mid = "yellow", high = "red", midpoint = 0) +
  theme_minimal() +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust =1, face = 'bold'),
        axis.text.y = element_text(face = 'bold', size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
```

``` r
prop_df <- table(scData$cellType, scData$sample_id) %>% 
  prop.table(margin = 2) %>%
  reshape2::melt() %>%
  `colnames<-`(c("CellType", "Sample", "Fraction"))

ggplot(prop_df, aes(x = Sample, y = Fraction, fill = CellType)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = cluster_cols) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, face = 'bold', hjust = 1)) 
```
