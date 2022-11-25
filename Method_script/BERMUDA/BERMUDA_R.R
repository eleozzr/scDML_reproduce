setwd(system("pwd", intern = T) )
#setwd("/Users/xiaokangyu/Desktop/scDML/scDML_project/BERMUDA/")
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(SingleCellExperiment)
  library(ggplot2)
  library(gdata)
  library(Seurat)
  library(getopt)
  library(stringr)
  source("./BERMUDA_raw_code/func_data.R")
  source("./BERMUDA_raw_code/2017-08-28-runMN-US.R")  
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

method="BERMUDA"
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
save_folder=paste0(savedir,"/",dataset,"/",method,"/")

print(Sys.time()-start.time)
# 
if(verbose){
  print("==================data information================")
  print(data)
  print("==================crosstable information==========")
  print(t(table(colData(data)$BATCH , colData(data)$celltype)))
  print("==================================================")
}

dat = counts(data)
# count to TMP
dat = apply(dat,2,function(x) (x*10^6)/sum(x))
# labels convert to factor,guarantee as.numeric to int
batch = as.factor(data@colData$BATCH)
cell = as.factor(data@colData$celltype)

##################write map_rule###################
map <- mapLevels(x=data@colData$BATCH)
map_df=t((do.call(rbind.data.frame,list(map))))
map_rule=data.frame(int=map_df[,1],batch=rownames(map_df))
print(map_rule)
write.csv(map_rule,paste0(save_folder,"bermuda_map_rule_batch.csv"),row.names=F)#
###############
map <- mapLevels(x=data@colData$celltype)
map_df=t((do.call(rbind.data.frame,list(map))))
map_rule=data.frame(int=map_df[,1],celltype=rownames(map_df))
print(map_rule)
write.csv(map_rule,paste0(save_folder,"bermuda_map_rule_celltype.csv"),row.names=F)## 
#############################################################

batch=as.numeric(batch)# convert factor to numeric
cell=as.numeric(cell)
#batch = unlist(lapply(batch,function(x) strtoi(substr(x, 6, 100))))
#cell = unlist(lapply(cell,function(x) strtoi(substr(x, 6, 100))))

# save simulated data
idx=list()
cnt=1
for(id in unique(batch)){
  idx[[cnt]]=which(batch==id)
  cnt=cnt+1
}
for (id in seq_len(length(unique(batch)))){
  idx_temp=idx[[id]]
  write_dataset(paste0(save_folder,dataset,"_batch_",id,".csv"),dat[,idx_temp],rep(id, ncol(dat[,idx_temp])), cell[idx_temp])
}

#memory.limit(size = 100000)
folder_name = save_folder
dataset_copy=dataset
dataset_names = list.files(folder_name,pattern = ".*batch_[0-9]\\.csv$")
dataset_names = str_sub(dataset_names,1,nchar(dataset_names)-4)
print(dataset_names)
dataset_list = list()
var_genes = list() #  a subset of highly variable genes
num_cells = 0


# Detect clusters in each dataset using Seurat
for (i in 1:length(dataset_names)) {
  filename = paste(folder_name, paste0(dataset_names[i], ".csv"), sep="/")
  print(paste0("Dataset: ", filename))
  # Seurat
  dataset = read_dataset(filename)
  dataset_list[[i]] = seurat_preprocessing(dataset, dataset_names[[i]])
  var_genes[[i]] = dataset_list[[i]]@var.genes
  num_cells = num_cells + dim(dataset_list[[i]]@data)[2]
}
var_genes = unique(unlist(var_genes))
for (i in 1:length(dataset_list)) {
  var_genes = intersect(var_genes, rownames(dataset_list[[i]]@data))
}

# combine datasets for metaneighbor
# log transformed TPM by Seurat, for metaneighbor
data = matrix(0, nrow = length(var_genes), ncol = num_cells) 
# labels, starting from 1
cluster_label_list = list()
dataset_label_list = list()
cell_idx = 1
cluster_idx = 1
for (i in 1:length(dataset_list)) {
  cell = dim(dataset_list[[i]]@data)[2]
  data[,cell_idx:(cell_idx+cell-1)] = as.matrix(dataset_list[[i]]@data[var_genes,])
  cluster_label_list[[i]] = as.integer(dataset_list[[i]]@meta.data$res.0.6) + cluster_idx
  cluster_idx = max(cluster_label_list[[i]]) + 1
  dataset_label_list[[i]] = rep(i, cell)
  cell_idx = cell_idx + cell
}

# write dataset with shifted cluster labels
# no cluster labels overlap between clusters 
for (i in 1:length(dataset_list)) {
  seurat_csv = paste(folder_name, paste0(dataset_names[i], "_seurat.csv"), sep="/")
  write_dataset_cluster(seurat_csv, as.matrix(dataset_list[[i]]@raw.data[var_genes,]), 
                        dataset_list[[i]]@meta.data$sample_labels,
                        dataset_list[[i]]@meta.data$cell_labels,
                        cluster_label_list[[i]])
}

# Metaneighbor
cluster_labels = unique(unlist(cluster_label_list))
rownames(data) = var_genes
pheno = as.data.frame(list(Celltype = as.character(unlist(cluster_label_list)),
                           Study_ID = as.character(unlist(dataset_label_list))),
                      stringsAsFactors=FALSE)
# run metaneighbor
cluster_similarity =run_MetaNeighbor_US(var_genes, data, cluster_labels, pheno)

# set cluster pairs from the same dataset to 0
for (i in 1:length(dataset_list)) {
  cluster_idx_tmp = unique(cluster_label_list[[i]])
  cluster_similarity[cluster_idx_tmp, cluster_idx_tmp] = 0
}

# order rows and columns
cluster_similarity = cluster_similarity[order(as.numeric(rownames(cluster_similarity))),]
cluster_similarity = cluster_similarity[,order(as.numeric(colnames(cluster_similarity)))]

# write out metaneighbor file
metaneighbor_file = paste(folder_name, paste0(dataset_copy, "_metaneighbor.csv"), sep="/")
write.table(cluster_similarity, metaneighbor_file, sep = ",", quote = F, col.names = T, row.names = F)
print("BERMUDA_R script running done....")



