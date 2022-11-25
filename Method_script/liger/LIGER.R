rm(list=ls())
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(SingleCellExperiment)
  library(scater)
  library(ggplot2)
  library(getopt)
  library(scran)
  library(batchelor)
  library(Seurat)
  library(SeuratWrappers)
  library(tidyverse)
  library(patchwork)
  library(cowplot)
  library(SeuratDisk)
  library(harmony)
  library(rliger)
})
start.time <- Sys.time()
spec <- matrix(
  c("dataset",  "d", 1, "character","dataset name", 
    "filepath", "f",1,"character","file folder of data", 
    'savedir',"s",1,"character","where to save result",
    "verbose", "v", 0, "integer","verbose information",
    "Save",  "w", 0, "integer","whethre to save preprocessed result",
    "help",   "h", 0, "logical","help information"
  ),
  byrow=TRUE, ncol=5)
opt <- getopt(spec=spec,debug=FALSE)
#print(opt)


# args <- commandArgs()
# print(args)

if ( is.null(opt$verbose)) {
  opt$verbose = 1
}
if ( is.null(opt$Save)) {
  opt$Save = 1
}

method="liger"
print(paste0("method=",method))
dataset=opt$dataset
filepath=opt$filepath
savedir=opt$savedir
verbose=opt$verbose
save= opt$Save


### create file
parent_dir=savedir
output_dir <- file.path(parent_dir, dataset)
# print(output_dir)
if (!dir.exists(output_dir)){
  dir.create(output_dir)
} else {
  print("Dir already exists!")
}

output_dir=file.path(parent_dir,dataset,method)
if (!dir.exists(output_dir)){
  dir.create(output_dir)
} else {
  print("Dir already exists!")
}

dataset_path=paste0(filepath,"/",dataset,"_raw.rds")
#print(dataset_path)
data=readRDS(dataset_path)

print("read data cost time:")
print(Sys.time()-start.time)
# 
if(verbose){
  print("==================data information================")
  print(data)
  print("==================crosstable information==========")
  print(t(table(colData(data)$BATCH , colData(data)$celltype)))
  print("==================================================")
}





print("===================Creating SeuraObject==========")
data_seurat=CreateSeuratObject(counts = counts(data),meta.data = as.data.frame(colData(data)))
print("===================Normalize SeuratObject=======")
data_seurat <- NormalizeData(data_seurat, verbose = FALSE)
print("===================Find HVG=====================")
data_seurat <- FindVariableFeatures(data_seurat, verbose = FALSE)
print("===================scale Data===================")
data_seurat <- ScaleData(data_seurat, split.by = "BATCH", do.center = FALSE)
print("===================Running PCA===================")
#data_seurat <- RunPCA(data_seurat, npcs = 30, verbose = F)
print("preprecessing done")
print("total cost time:")
print(Sys.time()-start.time)


data_seurat <- RunOptimizeALS(data_seurat, k = 20, lambda = 5, split.by = "BATCH")
data_seurat <- RunQuantileNorm(data_seurat, split.by = "BATCH")

#data_seurat <- FindNeighbors(data_seurat, reduction = "iNMF", dims = 1:20)
#data_seurat <- FindClusters(data_seurat, resolution = 0.3)
# Dimensional reduction and plotting


if(save){
  print("========================Running UMAP==================")
  data_seurat <- RunUMAP(data_seurat, dims = 1:ncol(data_seurat[["iNMF"]]), reduction = "iNMF")
  print("========================Visulize UMAP=================")
  p1=DimPlot(data_seurat, reduction = "umap", group.by = "BATCH", label.size = 10)+ggtitle("Integrated Batch")
  p2=DimPlot(data_seurat, reduction = "umap", group.by = "celltype",label.size = 10)+ggtitle("Integrated Celltype")
  p= p1 + p2 
  print("save figure")
  fig_save_path=paste0(savedir,dataset,"/",method,"/",dataset,"_",method,".png")
  ggsave(fig_save_path, p, width=12, height=5)
  print("save corrected data")
  data_save_path=paste0(savedir,dataset,"/",method,"/",dataset,"_",method,"_corrected")
  data_seurat$BATCH=as.character(data_seurat$BATCH)
  data_seurat$celltype=as.character(data_seurat$celltype)
  #########################################################################
  SaveH5Seurat(data_seurat, filename = paste0(data_save_path,".h5Seurat"),verbose = F,overwrite = T)
  Convert(paste0(data_save_path,".h5Seurat"), dest = "h5ad",overwrite = T,verbose=F)
  print("done")
}


