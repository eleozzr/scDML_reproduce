echo "BERMUDA python script !!!!!!"

log_dir="/DATA2/zhangjingxiao/yxk/scDML_project/Log"
project_dir="/DATA2/zhangjingxiao/yxk/scDML_project/"
dataset="bct_multi"
filepath="/DATA2/zhangjingxiao/yxk/dataset/bct/"

BERMUDApath="/home/hanfang_yang/anaconda3/bin/python"


echo $dataset

if [[ ! -e ${log_dir}/${dataset} ]]; then
    mkdir ${log_dir}/${dataset}
elif [[ ! -d ${log_dir}/${dataset} ]]; then
    echo "${log_dir}/${dataset} already exists but is not a directory"
fi

method="BERMUDA"
echo $method
cd ${project_dir}/${method}
echo $(pwd)
Log_out=${dataset}"_"${method}"_python_log.txt"

nohup ${BERMUDApath} run_BERMUDA.py --dataset $dataset --filepath $filepath --savedir ${project_dir}/evaluation/ > ${log_dir}/${dataset}/${Log_out} 2>&1&
echo $!






echo "done"



