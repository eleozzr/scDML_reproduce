rm(list=ls())
suppressPackageStartupMessages({
    library(iSMNN)
    library(Seurat)
    library(cowplot)
    library(ggplot2)
    library(getopt)
})

### read dataset ###

start.time <- Sys.time()
spec <- matrix(
  c("dataset",  "d", 1, "character","dataset name", ## 
    "filepath", "f",1,"character","file folder of data", ## 
    'savedir',"s",1,"character","where to save result",##
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

method="iSMNN"
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

print("read data cost:")
print(Sys.time()-start.time)
# 
if(verbose){
  print("==================data information================")
  print(data)
  print("==================crosstable information==========")
  print(t(table(colData(data)$BATCH , colData(data)$celltype)))
  print("==================================================")
}
print(data)

matched.clusters=names(table(data$celltype))
#########################
batch_list=list()
batch.cluster.labels=list()
id=1
count_matrix=NULL
batch_id=NULL
#########################
for(batch_names in names(table(data$BATCH))){
  batch_list[[id]]=data[,data$BATCH==batch_names]
  batch.cluster.labels[[id]]=colData(batch_list[[id]])$celltype      
  names(batch.cluster.labels[[id]])=colnames(batch_list[[id]])
  count_matrix=cbind(count_matrix,counts(batch_list[[id]]))
  batch_id=c(batch_id,rep(batch_names,ncol(counts(batch_list[[id]]))))
  id=id+1
}

merge <- CreateSeuratObject(counts = count_matrix, project = "merge", min.cells = 0, min.features = 0)
names(batch_id) <- colnames(merge)
merge <- AddMetaData(object = merge, metadata = batch_id, col.name = "batch_id")
merge.list <- SplitObject(merge, split.by = "batch_id")


merge.list <- lapply(X = merge.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

corrected.results <- iSMNN(object.list = merge.list, batch.cluster.labels = batch.cluster.labels, matched.clusters = matched.clusters,
                           strategy = "Short.run", iterations = 5, dims = 1:20, npcs = 30, k.filter = 30)
# change k.filter to 15,otherwise it throw error 

p1 <- DimPlot(corrected.results, reduction = "umap", group.by = "batch_id") 
p2 <- DimPlot(corrected.results, reduction = "umap", group.by = "cell.anno")
p=p1+p2
#print(p)

#save figure and corrected.results
print("save corrected figure")
fig_save_path=paste0(savedir,dataset,"/",method,"/",dataset,"_",method,".png")
ggsave(fig_save_path, p, width=12, height=5)

print("save corrected data")
data_save_path=paste0(savedir,dataset,"/",method,"/",dataset,"_",method,"_corrected")
#print(data_save_path)
saveRDS(corrected.results,paste0(data_save_path,".rds"))



