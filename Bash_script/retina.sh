echo "scDML_porject running,add ncelltype parmater"

dataset="retina"
n_cluster=14

log_dir="/DATA2/zhangjingxiao/yxk/scDML_project/Log"
project_dir="/DATA2/zhangjingxiao/yxk/scDML_project/"
filepath="/DATA2/zhangjingxiao/yxk/dataset/"${dataset}"/"
scANpath="/home/zhangjingxiao/.conda/envs/DESC/bin/python"
INSCTpath="/home/zhangjingxiao/.conda/envs/DESC/bin/python"
bbknnpath="/home/zhangjingxiao/.conda/envs/DESC/bin/python"
scDMLpath="/home/zhangjingxiao/.conda/envs/scDML/bin/python"
Seurat4path="/home/zhangjingxiao/.conda/envs/R4/bin/"
scVIpath="/home/zhangjingxiao/.conda/envs/scvi-env/bin/python"
Seurat2path="/home/zhangjingxiao/.conda/envs/Seurat2/bin/"
BERMUDApath="/home/zhangjingxiao/.conda/envs/BERMUDA/bin/python"

echo $dataset

if [[ ! -e ${log_dir}/${dataset} ]]; then
    mkdir ${log_dir}/${dataset}
elif [[ ! -d ${log_dir}/${dataset} ]]; then
    echo "${log_dir}/${dataset} already exists but is not a directory"
fi

method="fastMNN"

echo $method
cd ${project_dir}

file="./evaluation/"${dataset}"/"${method}"/"${dataset}"_"${method}"_corrected.h5ad"
if [ -f "$file" ]
then
	echo "$file found.$method had been already runned"
else
	echo "$file not found.rerun this method for this method"

    Log_out=${dataset}"_"${method}"_log.txt"
    nohup ${Seurat4path}/Rscript ${method}/fastMNN.R -d $dataset -f ${filepath} -s ${project_dir}/evaluation/ > ${log_dir}/${dataset}/${Log_out} 2>&1&
    echo $!
fi



method="harmony"
echo $method
cd ${project_dir}

file="./evaluation/"${dataset}"/"${method}"/"${dataset}"_"${method}"_corrected.h5ad"
if [ -f "$file" ]
then
	echo "$file found.$method had been already runned"
else
	echo "$file not found.rerun this method for this method"

    Log_out=${dataset}"_"${method}"_log.txt"
    nohup ${Seurat4path}/Rscript ${method}/harmony.R -d $dataset -f ${filepath} -s ${project_dir}/evaluation/ > ${log_dir}/${dataset}/${Log_out} 2>&1&
    echo $!
fi


method="Seurat3"
echo $method
cd ${project_dir}


file="./evaluation/"${dataset}"/"${method}"/"${dataset}"_"${method}"_corrected.h5ad"
if [ -f "$file" ]
then
	echo "$file found.$method had been already runned"
else
	echo "$file not found.rerun this method for this method"
    Log_out=${dataset}"_"${method}"_log.txt"
    nohup ${Seurat4path}/Rscript ${method}/Seurat3.R -d $dataset -f ${filepath} -s ${project_dir}/evaluation/ > ${log_dir}/${dataset}/${Log_out} 2>&1&
    echo $!
fi


method="liger"
echo $method
cd ${project_dir}

file="./evaluation/"${dataset}"/"${method}"/"${dataset}"_"${method}"_corrected.h5ad"
if [ -f "$file" ]
then
	echo "$file found.$method had been already runned"
else
	echo "$file not found.rerun this method for this method"
    Log_out=${dataset}"_"${method}"_log.txt"
    nohup ${Seurat4path}/Rscript ${method}/LIGER.R -d $dataset -f ${filepath} -s ${project_dir}/evaluation/ > ${log_dir}/${dataset}/${Log_out} 2>&1&
    echo $!
fi


method="Scanorama"
echo $method
cd ${project_dir}

file="./evaluation/"${dataset}"/"${method}"/"${dataset}"_"${method}"_corrected.h5ad"
if [ -f "$file" ]
then
	echo "$file found.$method had been already runned"
else
	echo "$file not found.rerun this method for this method"

    Log_out=${dataset}"_"${method}"_log.txt"
    nohup ${scANpath} ${method}/Scanorama.py --dataset $dataset --filepath $filepath --savedir ${project_dir}/evaluation/ > ${log_dir}/${dataset}/${Log_out} 2>&1&
    echo $!
fi


method="INSCT"
echo $method
cd ${project_dir}

file="./evaluation/"${dataset}"/"${method}"/"${dataset}"_"${method}"_corrected.h5ad"
if [ -f "$file" ]
then
	echo "$file found.$method had been already runned"
else
	echo "$file not found.rerun this method for this method"
    
    Log_out=${dataset}"_"${method}"_log.txt"
    nohup ${INSCTpath} ${method}/INSCT.py --dataset $dataset --filepath $filepath --savedir ${project_dir}/evaluation/ > ${log_dir}/${dataset}/${Log_out} 2>&1&
    echo $!
fi


method="BBKNN"
echo $method
cd ${project_dir}

file="./evaluation/"${dataset}"/"${method}"/"${dataset}"_"${method}"_corrected.h5ad"
if [ -f "$file" ]
then
	echo "$file found.$method had been already runned"
else
	echo "$file not found.rerun this method for this method"
    Log_out=${dataset}"_"${method}"_log.txt"
    nohup ${bbknnpath} ${method}/BBKNN.py --dataset $dataset --filepath $filepath --savedir ${project_dir}/evaluation/ > ${log_dir}/${dataset}/${Log_out} 2>&1&
    echo $!
fi


method="scDML"
echo $method
cd ${project_dir}

file="./evaluation/"${dataset}"/"${method}"/"${dataset}"_"${method}"_corrected.h5ad"
if [ -f "$file" ]
then
	echo "$file found.$method had been already runned"
else
	echo "$file not found.rerun this method for this method"
    Log_out=${dataset}"_"${method}"_log.txt"
    nohup ${scDMLpath} scDML_script.py --dataset $dataset --filepath $filepath --savedir ${project_dir}/evaluation/ --ncelltype ${n_cluster} > ${log_dir}/${dataset}/${Log_out} 2>&1&
    echo $!
fi



method="scVI"
echo $method
cd ${project_dir}

file="./evaluation/"${dataset}"/"${method}"/"${dataset}"_"${method}"_corrected.h5ad"
if [ -f "$file" ]
then
	echo "$file found.$method had been already runned"
else
	echo "$file not found.rerun this method for this method"
    Log_out=${dataset}"_"${method}"_log.txt"
    nohup ${scVIpath} ${method}/scVI.py --dataset $dataset --filepath $filepath --savedir ${project_dir}/evaluation/ > ${log_dir}/${dataset}/${Log_out} 2>&1&
    echo $!
fi


method="BERMUDA"

echo $method

cd ${project_dir}

file="./evaluation/"${dataset}"/"${method}"/"${dataset}"_"${method}"_corrected.h5ad"
if [ -f "$file" ]
then
	echo "$file found.$method had been already runned"
else
	echo "$file not found.rerun this method for this method"
    echo "BERMUDA_R_preprocessing"
    cd ${project_dir}/${method}
    echo $(pwd)

    Log_out=${dataset}"_"${method}"_log_R.txt"

    nohup ${Seurat2path}/Rscript BERMUDA_R.R -d $dataset -f ${filepath} -s ${project_dir}/evaluation/ > ${log_dir}/${dataset}/${Log_out} 2>&1&
    echo $!

    while [ -e  /proc/$!/statm ]; do
        sleep 1
    done

    if [[ ! -e ${log_dir}/${dataset} ]]; then
        mkdir ${log_dir}/${dataset}
    elif [[ ! -d ${log_dir}/${dataset} ]]; then
        echo "${log_dir}/${dataset} already exists but is not a directory"
    fi

    method="BERMUDA"
    echo "BERMUDA_python_batch_correction"
    echo $method
    cd ${project_dir}/${method}
    echo $(pwd)
    Log_out=${dataset}"_"${method}"_log_python.txt"

    nohup ${BERMUDApath} run_BERMUDA.py --dataset $dataset --filepath $filepath --savedir ${project_dir}/evaluation/ > ${log_dir}/${dataset}/${Log_out} 2>&1&
    echo $!

    echo "bermuda running done,all done!!!!"
fi



