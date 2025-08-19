#!/usr/bin/bash
##########################################################
# Creater       :  夜下凝月
# Created  date :  2023-09-08, 16:42:31
# Modiffed date :  2023-09-08, 16:42:31
##########################################################


( [ $# -ne 3 ] ) && \
    echo -e "\n\tUsage:\t$0 <fasta|.gz> <depth> <path_output_prefix>\n" && exit 127

in_fa=`realpath -s $1` # file
depth=`realpath -s $2` # file
pref=`realpath -s $3` # prefix of output. (relative/absolute path)

pref_name=${pref##*/} # 输出文件的名字

#------------------
set -e 
set -o pipefail
shopt -s expand_aliases

alias pigz='/usr/local/bin/pigz'
alias bzip2='/usr/bin/bzip2'
conda_env='/home/zhangy2/.conda/envs/metabinner'

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

metabinner_dir="/share/data1/software/MetaBinner-master/"

( [ ! -d ${pref} ] ) && { mkdir -p "${pref}" || exit 127 ; }
cd ${pref}

source activate $conda_env  || exit 127;

## link fasta
if [ ! -L ${pref_name}.tmp.fa ] || [ ! -f ${pref_name}.tmp.fa ];then
    if [[ $in_fa =~ ".gz" ]];then
        pigz -dc ${in_fa} > ${pref_name}.tmp.fa || exit 127;
    elif [[ $in_fa =~ ".bz2" ]];then
        bzip2 -dc ${in_fa} > ${pref_name}.tmp.fa || exit 127;
    else
        ln -s ${in_fa} ${pref_name}.tmp.fa || exit 127;
    fi
fi

## format depth
( [ ! -f ${pref_name}.tmp.depth ] ) && { cut -f 1,4 ${depth} > ${pref_name}.tmp.depth || exit 127; }
## generate kmer
( [ ! -f ${pref_name}.kmer.ok ] ) && { python ${metabinner_dir}/scripts/gen_kmer.py.zy ${pref_name}.tmp.fa 0 4 ${pref_name}.kmer.csv && touch ${pref_name}.kmer.ok || ! rm ${pref_name}.kmer.csv || exit 127 ;}
## get cluster
( [ ! -f ${pref_name}.logtrans.tsv ] ) && { python ${metabinner_dir}/scripts/component_binning.py.zy \
    --software_path ${metabinner_dir} \
    --contig_file ${pref_name}.tmp.fa \
    --coverage_profiles ${pref_name}.tmp.depth \
    --composition_profiles ${pref_name}.kmer.csv \
    --output ./ --log ${pref_name}.out.log \
    --threads 80 \
    --dataset_scale large \
    --contig_length_threshold 500 \
    && touch ${pref_name}.logtrans.ok \
    || exit 127 ; }

## get bins
( [ ! -f ${pref_name}.bin.ok ] ) && { 
    python ${metabinner_dir}/scripts/gen_bins_from_tsv.py \
        -f ${pref_name}.tmp.fa \
        -r ${pref}/intermediate_result/kmeans_length_weight_X_t_logtrans_result.tsv \
        -o ${pref_name}.bins \
        && mv ${pref_name}.bins ${pref}.bins \
        && touch ${pref_name}.bin.ok \
        && cd .. \
        && touch ${pref}.bins.ok \
        && rm ${pref}.bins.running \
        || exit 127 ; }

## remove tmp
( [ -f ${pref_name}.bins.ok ] ) && { rm -rf ${pref} || exit 127; }
