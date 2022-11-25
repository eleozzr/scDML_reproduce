import scanpy as sc 
import matplotlib
import argparse
matplotlib.use("Agg")
# import os 
# os.system("clear")

method="Raw"
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

dataset_path=filepath+"/"+dataset+"_raw.h5ad"# 
adata=sc.read(dataset_path)
sc.settings.figdir=save_figdir+dataset+"/"+method+"/"

sc.pp.normalize_total(adata,target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata,n_top_genes=2000,subset=True) # use top 2000 genes
sc.pp.scale(adata)
sc.tl.pca(adata)
sc.pp.neighbors(adata,random_state=0)
sc.tl.umap(adata)
sc.pl.umap(adata,color=["BATCH","celltype"],save="_"+dataset+"_"+method+"_raw.png")
del adata.raw ## if dont delte adata.raw, adata.write() will throw errors
adata.write(savedir+dataset+"/"+method+"/"+dataset+"_"+method+"_raw.h5ad")
print("done")