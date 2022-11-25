rm(list=ls())

suppressPackageStartupMessages({
    library(iSMNN)
    library(Seurat)
    library(cowplot)
    library(ggplot2)
    library(getopt)
})

### read dataset ###
dataset="bct"
read_dir=paste0("/DATA2/zhangjingxiao/yxk/dataset/",dataset,"/",dataset,"_raw.rds")
data=readRDS(read_dir)
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

p1 <- DimPlot(corrected.results, reduction = "umap", group.by = "batch_id") 
p2 <- DimPlot(corrected.results, reduction = "umap", group.by = "cell.anno")
p=p1+p2
print(p)


#save figure and corrected.results
print("save corrected figure")
fig_save_path=paste0("/DATA2/zhangjingxiao/yxk/scDML_project/iSMNN/",dataset,"_iSMNN_run.png")
ggsave(fig_save_path, p, width=12, height=5)

# save corrected.results




