#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 23:42:14 2021

@author: xiaokangyu
"""

import argparse
from time import time
import scanpy as sc
import numpy as np
from scanorama import correct, visualize, process_data
from scanorama import dimensionality_reduce
import scanorama
from datetime import timedelta
import scipy
import pickle
import os
import matplotlib
matplotlib.use('Agg')

method="Scanorama"

x0=time()         
parser = argparse.ArgumentParser()
parser.add_argument('--dataset',type=str,required=True, help='dataset name')
parser.add_argument("--filepath",type=str,required=True,help="folder path of stored data")
parser.add_argument("--verbose",type=bool,default=True,help="print additional information")
parser.add_argument("--savedir",type=str,required=True,help="where to save data")
parser.add_argument("--save",type=bool,default=True,help="whether to save data")
print("method=",method)
args=parser.parse_args()  
print("dataset=",args.dataset)
dataset=args.dataset
filepath=args.filepath
verbose=args.verbose
savedir=args.savedir
save=args.save
save_figdir=savedir

sc.settings.figdir=save_figdir+dataset+"/"+method+"/"
dataset_path=filepath+"/"+dataset+"_raw.h5ad"# 目前不用加入raw
adata=sc.read(dataset_path)
print("read data cost",time()-x0,"s")

## 还需要创建文件夹
if not os.path.exists(sc.settings.figdir):
    os.makedirs(sc.settings.figdir)

###############format converting########################
if issubclass(type(adata.X), scipy.sparse.csc.csc_matrix):
    adata.X=adata.X.tocsr()

####log1p
print("====normalize data and log1p===================")
sc.pp.normalize_per_cell(adata,counts_per_cell_after=1e4)# 
sc.pp.filter_genes_dispersion(adata,n_top_genes =1000) 
sc.pp.log1p(adata)

adata_ls=[]
for batch in np.unique(adata.obs["BATCH"]):
    sep_batch=adata[adata.obs["BATCH"]==batch,:].copy()
    adata_ls.append(sep_batch)

# if(save):
#     print("saving data to file")
#     write_path=savedir+dataset+"/"+method+"/"+dataset+"_"+method+"_preprocessed.pkl" 
#     #############################判断存储结果的文件夹是否存在，如果不存在，则创建##########
#     #pdb.set_trace()
# #     temp_dir = os.path.join(os.getcwd(), "../dataset/preprocessed_dataset/"+dataset+"/"+method+"/")
# #     if not os.path.exists(temp_dir): 
# #         os.makedirs(temp_dir)
# #     temp_dir = os.path.join(os.getcwd(), "../dataset/corrected_dataset/"+dataset+"/"+method+"/")
# #     if not os.path.exists(temp_dir): 
# #         os.makedirs(temp_dir)
# #     temp_dir = os.path.join(os.getcwd(), "../result/"+dataset+"/"+method+"/")
# #     if not os.path.exists(temp_dir): 
# #         os.makedirs(temp_dir)
#     ################################################################################

#     with open(write_path,"wb") as f:  # Python 3: open(..., 'rb')
#         pickle.dump(adata_ls, f)

print("preprocess dataset total cost",time()-x0,"s")
print("read data cost",time()-x0,"s")    
corrected = scanorama.correct_scanpy(adata_ls, batch_size=50, return_dense=True, knn=20)
adata_corrected = sc.AnnData(np.concatenate([corrected[i].X for i in range(len(corrected))]))
print("correct time cost",time()-x0)

if(save):
    #pdb.set_trace()
    adata=sc.AnnData.concatenate(*adata_ls)
    adata_corrected = sc.AnnData(np.concatenate([corrected[i].X for i in range(len(corrected))]))
    adata_corrected.obs['BATCH'] =adata.obs["BATCH"].values
    adata_corrected.obs['celltype'] = adata.obs["celltype"].values
    print("======================Visulization after Batch Effect Corecttion using Scanorama=============================")
    sc.tl.pca(adata_corrected)
    sc.pp.neighbors(adata_corrected)
    sc.tl.umap(adata_corrected)
    sc.pl.umap(adata_corrected,color=["BATCH","celltype"],save="_"+dataset+"_"+method+"_corrected.png")
    adata_corrected.write(savedir+dataset+"/"+method+"/"+dataset+"_"+method+"_corrected.h5ad")
    
print("done")




