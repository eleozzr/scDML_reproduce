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

method="Scanorama_raw"

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
dataset_path=filepath+"/"+dataset+"_raw.h5ad"#
adata0=sc.read(dataset_path)
print("read data cost",time()-x0,"s")


if not os.path.exists(sc.settings.figdir):
    os.makedirs(sc.settings.figdir)

###############format converting########################
if issubclass(type(adata0.X), scipy.sparse.csc.csc_matrix):
    adata0.X=adata0.X.tocsr()

adata=adata0.copy()
sc.pp.normalize_per_cell(adata,counts_per_cell_after=1e4)# 
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata,n_top_genes=2000,subset=True)

adata_ls=[]
for batch in np.unique(adata.obs["BATCH"]):
    sep_batch=adata0[adata0.obs["BATCH"]==batch,adata.var_names].copy()
    adata_ls.append(sep_batch)
    
corrected = scanorama.correct_scanpy(adata_ls, return_dimred=True)


if(save):
    #pdb.set_trace()
    adata_corrected=sc.concat(corrected)
    sc.pp.neighbors(adata_corrected,use_rep="X_scanorama")
    sc.tl.umap(adata_corrected)
    sc.pl.umap(adata_corrected,color=["BATCH","celltype"],save="_"+dataset+"_"+method+"_corrected.png")
    adata_corrected.write(savedir+dataset+"/"+method+"/"+dataset+"_"+method+"_corrected.h5ad")
    
print("done")




