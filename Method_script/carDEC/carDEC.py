import pandas as pd
import os
import numpy as np
import pickle
from copy import deepcopy
from shutil import move
import warnings
import argparse
from time import time
import sklearn.metrics as metrics
from sklearn.metrics import adjusted_rand_score as ari, normalized_mutual_info_score as nmi
import scanpy as sc
from anndata import AnnData
from CarDEC import CarDEC_API

import os
import matplotlib
matplotlib.use('Agg')

method="carDEC"

x0=time()         
parser = argparse.ArgumentParser()
parser.add_argument('--dataset',type=str,required=True, help='dataset name')
parser.add_argument("--filepath",type=str,required=True,help="folder path of stored data")
parser.add_argument("--ncelltype",type=int,required=True,help="number of celltype in dataset")
parser.add_argument("--verbose",type=bool,default=True,help="print additional information")
parser.add_argument("--savedir",type=str,required=True,help="where to save data")
parser.add_argument("--save",type=bool,default=True,help="whether to save data")
print("method=",method)
args=parser.parse_args()  
print("dataset=",args.dataset)
dataset=args.dataset
filepath=args.filepath
ncelltype=args.ncelltype
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

if(adata.shape[1]>2000):
    num_hvg=2000
else:
    num_hvg=int(adata.shape[1]*0.9)

print(adata)
CarDEC = CarDEC_API(adata, weights_dir = "evaluation/"+dataset+"/"+method+"/weight", batch_key = "BATCH", n_high_var = num_hvg, LVG = True)
print("load carDEC model")

CarDEC.build_model(n_clusters = ncelltype)
print("build carDEC model")
CarDEC.make_inference()
print("make inference")
embedded = deepcopy(CarDEC.dataset.obsm['embedding']) #The latent embedding numpy array

q = deepcopy(CarDEC.dataset.obsm['cluster memberships']) #The cluster membership numpy array
labels = np.argmax(q, axis=1)
labels = [str(x) for x in labels]

adata_corrected = AnnData(embedded)
adata_corrected.obs["celltype"] = list(CarDEC.dataset.obs['celltype'])
adata_corrected.obs["reassign_label"] = list(labels)
adata_corrected.obs["BATCH"] = list(CarDEC.dataset.obs['BATCH'])

sc.tl.pca(adata_corrected)
sc.pp.neighbors(adata_corrected)
sc.tl.umap(adata_corrected)
sc.pl.umap(adata_corrected,color=["BATCH","celltype","reassign_label"],save="_"+dataset+"_"+method+"_corrected.png")
print("save figure...")
adata_corrected.write(savedir+dataset+"/"+method+"/"+dataset+"_"+method+"_corrected.h5ad")
print("save adata_corrected...")
print("done...")
