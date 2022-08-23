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


dataset_path=filepath+"/"+dataset+"_raw.h5ad"# 
adata=sc.read(dataset_path)

# method="scDML1.0"
# scdml1=scDMLModel(adata,mode="unsupervised",batch_key="BATCH",celltype_key="celltype",save_dir=save_figdir+dataset+"/"+method+"/")
# scdml1.full_run(resolution=1.0,mode="unsuprevised",fixed_ncluster=ncelltype,expect_num_cluster=ncelltype,do_plot=False)   

# method="scDML2.0"
# scdml2=scDMLModel(adata,mode="unsupervised",batch_key="BATCH",celltype_key="celltype",save_dir=save_figdir+dataset+"/"+method+"/")
# scdml2.full_run(resolution=2.0,mode="unsuprevised",fixed_ncluster=ncelltype,expect_num_cluster=ncelltype,do_plot=False)

method="scDML"
scdml3=scDMLModel(adata,mode="unsupervised",batch_key="BATCH",celltype_key="celltype",save_dir=save_figdir+dataset+"/"+method+"/")
scdml3.full_run(resolution=3.0,mode="unsuprevised",fixed_ncluster=ncelltype,expect_num_cluster=ncelltype,do_plot=False,flag=dataset)




