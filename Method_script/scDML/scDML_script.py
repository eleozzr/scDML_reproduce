import scanpy as sc 
import matplotlib
import argparse
matplotlib.use("Agg")
from scDML import scDMLModel
# import os 
# os.system("clear")

method="scDML"
parser = argparse.ArgumentParser()
parser.add_argument('--dataset',type=str,required=True, help='dataset name')
parser.add_argument("--filepath",type=str,required=True,help="folder path of stored data")
parser.add_argument("--ncelltype",type=int,required=True,help="number of celltype in dataset")
parser.add_argument("--K_in",type=int,default=5,help="K value to calculate KNN pair")
parser.add_argument("--K_bw",type=int,default=10,help="K value to calculate MNN pair")
parser.add_argument("--cluster_method",type=str,default="louvain",help="clustering algorithm to initize cluster label")
parser.add_argument("--resolution",type=float,default=3.0,help="resolution of clustering algorithm")
parser.add_argument("--n_hvg",type=int,default=1000,help="number of highly variable genes to be selected")
parser.add_argument("--verbose",type=bool,default=True,help="print additional information")
parser.add_argument("--savedir",type=str,required=True,help="where to save data")
parser.add_argument("--save",type=bool,default=True,help="whether to save data")
print("method=",method)
args=parser.parse_args()  
print("dataset=",args.dataset)
dataset=args.dataset
filepath=args.filepath
ncelltype=args.ncelltype
K_in=args.K_in
K_bw=args.K_bw
cluster_method=args.cluster_method
resolution=args.resolution
n_hvg=args.n_hvg
verbose=args.verbose
savedir=args.savedir
save=args.save
save_figdir=savedir

dataset_path=filepath+"/"+dataset+"_raw.h5ad"# 
adata_raw=sc.read(dataset_path)
sc.settings.figdir=save_figdir+dataset+"/"+method+"/"

method="scDML"
scdml=scDMLModel(save_dir=save_figdir+dataset+"/"+method+"/")

adata=scdml.preprocess(adata_raw,cluster_method=cluster_method,resolution=resolution,n_high_var=n_hvg)
#print(adata)
scdml.integrate(adata,batch_key="BATCH",ncluster_list=[ncelltype],K_in=K_in,K_bw=K_bw,
               expect_num_cluster=ncelltype,merge_rule="rule2")
adata.obs["cluster_celltype"]=adata.obs["reassign_cluster"].copy()
# visulization 
#####################################################
sc.pp.neighbors(adata,random_state=0,use_rep="X_emb")
sc.tl.umap(adata)
#####################################################

sc.pl.umap(adata,color=["BATCH","celltype"],save="_"+dataset+"_"+method+"_corrected.png")
adata.write(savedir+dataset+"/"+method+"/"+dataset+"_"+method+"_corrected.h5ad")




