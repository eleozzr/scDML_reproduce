#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 23:42:14 2021

@author: xiaokangyu
"""
import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
import argparse
from time import time
import scanpy as sc
import numpy as np
from datetime import timedelta
import scipy
import os
from tnn.tnn import *
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
import matplotlib
matplotlib.use('Agg')

method="INSCT"
        
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
adata_insct=sc.read(dataset_path)
#print("read data cost",time()-x0,"s")

if not os.path.exists(sc.settings.figdir):
    os.makedirs(sc.settings.figdir)

###############format converting########################
# if issubclass(type(adata_insct.X), scipy.sparse.csc.csc_matrix):
#     adata_insct.X=adata_insct.X.tocsr()

####log1p
print("====normalize data and log1p===================")
sc.pp.normalize_total(adata_insct, target_sum=1e4)
sc.pp.log1p(adata_insct)
sc.pp.highly_variable_genes(adata_insct, n_top_genes=1000, subset = True)
sc.tl.pca(adata_insct)
#print("preprocess dataset total cost",time()-x0,"s")


model=TNN(k=20)
model.fit(X = adata_insct, batch_name = "BATCH", shuffle_mode = True)
adata_insct.obsm['X_insct'] = model.transform(adata_insct)

adata_corrected=sc.AnnData(adata_insct.obsm['X_insct'])
adata_corrected.obs['BATCH'] =adata_insct.obs["BATCH"].values
adata_corrected.obs['celltype'] = adata_insct.obs["celltype"].values
print(adata_corrected)
print("======================Visulization after Batch Effect Corecttion using INSCT=============================")
adata_corrected.obsm["X_emb"]=adata_corrected.X
sc.pl.embedding(adata_corrected,basis="emb",color=["BATCH","celltype"],save="_"+dataset+"_"+method+"_corrected.png")
adata_corrected.write(savedir+dataset+"/"+method+"/"+dataset+"_"+method+"_corrected.h5ad")

#################################################################################

print("done")





