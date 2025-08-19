#!/usr/bin/bash

#------------------
set -e 
set -o pipefail
shopt -s expand_aliases

#----------------------------
#           configures
conda_env='/share/data1/yuhail/software/semibin2'
alias ls='/usr/bin/ls'
alias exit='builtin exit'

if [ $# -ne 7 ];then
    echo -e "\n\t$0 <environment> <depth> <fasta> <path_output_prefix> <seed:2024> <minLengthContig:2000> <minLengthBin:200000>"
    echo -e "\n\tenvironment:"
    echo -e "\n\t\thuman_gut/dog_gut/ocean/soil/cat_gut/human_oral/mouse_gut/pig_gut/built_environment/wastewater/chicken_caecum/global"
    echo -e "\n"
    exit 127
fi

#---------------------------
#      parse parameters
env=$1
dep=`realpath -s $2`
infa=`realpath -s $3`
pref=`realpath -s $4`
seed=$5
mc=$6
mb=$7

pref_name=${pref##*/} # 输出文件的名字
mymb=`echo ${mb}/1000 | bc`


#------------------
clean_tmp(){
    echo "Accidental termination"
    rm ${pref}.bins.running
    exit 127
}

trap 'clean_tmp' SIGINT SIGTERM

## 判断是否已经有了结果，或者正在运行
if [ -f ${pref}.bins.running ];then
    echo "It's running." && exit 127;
elif [ -f ${pref}.bins.ok ];then
    echo "exists ${pref}.bins" && exit 0;
fi

touch ${pref}.bins.running

## 创建输出目录
( [ ! -d ${pref} ] ) && { mkdir -p "${pref}" || exit 127 ; }
cd ${pref}

source activate ${conda_env} || exit 127;


SemiBin2 single_easy_bin  \
    -i ${infa} \
    -m ${mc} \
    -o ${pref} \
    --minfasta-kbs ${mymb} \
    --environment ${env} \
    --random-seed ${seed} \
    --depth-metabat2 ${dep} > ${pref}.log 2>&1 && touch ${pref}.bins.ok || exit 127;

## 整理分箱的结果
ls output_bins/*.gz | perl -ne 'chomp;$_=~/SemiBin_(\d+).fa/;print "mv $_ $1.fa.gz\n"' | sh
chmod 444 *.fa.gz

## 清理中间文件
rm -r  ${pref}.bins.running contig_bins.tsv data.csv recluster_bins_info.tsv SemiBinRun.log output_bins/

