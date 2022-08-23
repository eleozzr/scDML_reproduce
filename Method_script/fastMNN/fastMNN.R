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
  c("dataset",  "d", 1, "character","dataset name", ## 第二个参数只能一个字符好像
    "filepath", "f",1,"character","file folder of data", ## 传入目标文件所在的路径
    'savedir',"s",1,"character","where to save result",## 结果存在哪里
    "verbose", "v", 0, "integer","verbose information",
    "Save",  "w", 0, "integer","whethre to save preprocessed result",
     "help",   "h", 0, "logical","help information"
    ),
  byrow=TRUE, ncol=5)
opt <- getopt(spec=spec,debug=FALSE)
#print(opt)
###################### 设置默认值 ################

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

################################## DEBUG #############################
#Rscript fastMNN/fastMNN.R -d "4batch_4celltype_multi" -f "/Users/xiaokangyu/Desktop/单细胞学习/单细胞数据集/splatter_sim/"
#-sd "/Users/xiaokangyu/Desktop/tDCA_project/evaluation/"
# method="fastMNN"
# print(paste0("method=",method))
# dataset="4batch_4celltype_multi"
# filepath="/Users/xiaokangyu/Desktop/单细胞学习/单细胞数据集/splatter_sim/"
# savedir="/Users/xiaokangyu/Desktop/tDCA_project/evaluation/" # 我应该把结果存在evluation的结果里，这样就能统一读了
# verbose=1
# save=1

######################################################################
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

print("读取数据用时:")
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
print("加上预处理数据总用时:")
print(Sys.time()-start.time)
# if(save){
#   print("saving data to file")
#   write_path=paste0(savedir,dataset,"/",method,"/",dataset,"_",method,"_preprocessed.rds")
#   saveRDS(scRNAlist,file=write_path)
#   print("加上存储文件总用时:")
#   print(Sys.time()-start.time)
#   }
print("save preprocessed dataset done")
print("run fastMNN correction....") 
scRNA <- RunFastMNN(object.list = scRNAlist,verbose=FALSE)
print("加上fastMNN矫正用时")
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



