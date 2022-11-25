#!/usr/bin/env python
# coding: utf-8

# In[1]:


import scvi
import argparse
import scanpy as sc
from time import time
import matplotlib.pyplot as plt
from datetime import timedelta
import scipy
import pickle
import os
sc.set_figure_params(figsize=(4, 4))
import matplotlib
matplotlib.use('Agg')

method="scVI"
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
dataset_path=filepath+"/"+dataset+"_raw.h5ad"
adata=sc.read(dataset_path)
print("read data cost",time()-x0,"s")

if not os.path.exists(sc.settings.figdir):
    os.makedirs(sc.settings.figdir)

sc.pp.filter_genes(adata, min_counts=3)
adata.layers["counts"] = adata.X.copy() # preserve counts
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata # freeze the state in `.raw`

if(dataset=="heart_140"): # sc.pp.hvg==bug
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=2000,
        subset=True,
        layer="counts",
        flavor="seurat_v3",
    )
else:
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=2000,
        subset=True,
        layer="counts",
        flavor="seurat_v3",
        batch_key="BATCH"
    )
    
# scvi.model.SCVI.setup_anndata(
#     adata,
#     layer="counts",
#     batch_key="BATCH",
# )
scvi.data.setup_anndata(adata, layer="counts", batch_key="BATCH")

model = scvi.model.SCVI(adata)
print(model)
model.train()
latent = model.get_latent_representation()
print(latent.shape)
adata_corrected=sc.AnnData(latent)
adata_corrected.obs['BATCH'] =adata.obs["BATCH"].values
adata_corrected.obs['celltype'] = adata.obs["celltype"].values
#adata_corrected.obsm["X_emb"]=adata.X.copy() ##
print("======================Visulization after Batch Effect Corecttion using scVI=============================")
sc.tl.pca(adata_corrected)
sc.pp.neighbors(adata_corrected)
sc.tl.umap(adata_corrected)
sc.pl.umap(adata_corrected,color=["BATCH","celltype"],save="_"+dataset+"_"+method+"_corrected.png")
adata_corrected.write(savedir+dataset+"/"+method+"/"+dataset+"_"+method+"_corrected.h5ad")
print("running scVI done.......")

