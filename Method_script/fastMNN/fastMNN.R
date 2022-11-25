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

method="fastMNN"
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
data_seurat=CreateSeuratObject(counts = counts(data),meta.data = as.data.frame(colData(data)),verbose=F)
print("===================Split SeuratOject============")
scRNAlist <- SplitObject(data_seurat, split.by = "BATCH")
print("===================Normalize SeuratObject=======")
scRNAlist <- lapply(scRNAlist, FUN = function(x) NormalizeData(x,verbose=F))
print("===================Find HVG=====================")
scRNAlist <- lapply(scRNAlist, FUN = function(x) FindVariableFeatures(x,verbose=F))
print("preprecessing done")
print("total cost time:")
print(Sys.time()-start.time)

print("save preprocessed dataset done")
print("run fastMNN correction....") 
scRNA <- RunFastMNN(object.list = scRNAlist,verbose=FALSE)
print("fastMNN correction cost time")
print(Sys.time()-start.time)

scRNA <- RunUMAP(scRNA, reduction = "mnn", dims = 1:30,verbose=FALSE)
scRNA <- FindNeighbors(scRNA, reduction = "mnn", dims = 1:30,verbose=FALSE)
scRNA <- FindClusters(scRNA,verbose=FALSE)
p1 <- DimPlot(scRNA, group.by = "BATCH",label.size = 10) + 
ggtitle("Integrated Batch")
p2 <- DimPlot(scRNA, group.by="celltype",label.size = 10) + 
ggtitle("Integrated Celltype")
p = p1 + p2
#print(p)
print("save corrected figure")
fig_save_path=paste0(savedir,dataset,"/",method,"/",dataset,"_",method,".png")
ggsave(fig_save_path, p, width=12, height=5)
print("save corrected data")
data_save_path=paste0(savedir,dataset,"/",method,"/",dataset,"_",method,"_corrected")
#print(data_save_path)
SaveH5Seurat(scRNA, filename = paste0(data_save_path,".h5Seurat"),verbose = F,overwrite = T)
Convert(paste0(data_save_path,".h5Seurat"), dest = "h5ad",overwrite = T,verbose=F)
print("done")



