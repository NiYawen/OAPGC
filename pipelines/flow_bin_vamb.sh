#!/usr/bin/bash
##########################################################
# Creater       :  夜下凝月
# Created  date :  2024-01-15, 14:43:35
# Modiffed date :  2024-01-15, 14:43:35
##########################################################

( [ $# -ne 3 ] ) && \
        echo -e "\n\tUsage:\t$0 <fasta> <depth> <path_output_prefix>\n" && exit 127

in_fa=`realpath -s $1` # file
depth=`realpath -s $2` # file
pref=`realpath -s $3` # prefix of output. (relative/absolute path)

pref_name=${pref##*/} # 文件名前缀
BASH_LOG=${pref}.run.log
ERROR_LOG=${pref}.error
LOG=${pref}.log

######################
#       Config
shopt -s expand_aliases
set -o pipefail # 如果管道符第一个出错，就报错

 
######################
#       alias
alias seqkit="/usr/local/bin/seqkit"
alias parallel="/usr/local/bin/parallel"
alias ls="/usr/bin/ls"
get_fasta='/share/data1/zhangy2/scripts/find_fasta_from_list.py'
minBin=200000
CONDA_ENV=/share/data1/zhangy2/conda/envs/vamb3


function check_file(){
    # 检查文件的有无
    inf=$1
    log=$2
    if [ ! -f $inf ];then
        echo -e "ERROR:\tNo file $inf" > $log && exit 127
    fi
}

function check_size(){
    # 检查fasta文件大小是否满足条件
    fa=$1
    minBin=$2
    log=$3
    size=`seqkit stats $fa | tail -n 1 | perl -ne 'chomp;$_=~s/,//g;@l=split/\s+/;print "$l[4]"'`
    if [ $size -lt $minBin ];then
        echo -e "ERROR:\tinput file too small. ( $size\t$fa )" > $log
        exit 127
    fi
}

#######################
#       check
check_size $in_fa $minBin ${ERROR_LOG}
check_file $in_fa ${ERROR_LOG}
check_file $depth ${ERROR_LOG}


( [ -f ${pref}.ok ] ) && echo "exists ${pref}" && exit 127
( [ ! -d ${pref} ] ) && { mkdir -p ${pref} || exit 127 ;}
( [ ${pref_name} == "" ] )  && { echo "输出必须是前缀，而不是目录!!!!!!" && exit 127;}

source activate ${CONDA_ENV} || exit 127;

cd ${pref};

cut -f 1 ${depth} | seqkit grep -f - ${in_fa} -o ${pref_name}.tmp.fa
seqkit fx2tab --id-regexp "(\S+)" -i -n -l ${pref_name}.tmp.fa > ${pref_name}.tmp.fa.len # 获取序列长度

batch_size=256
while [[ 2 -gt 1 ]] ;
do
    # run vamb
    if [ -f vamb.o/clusters.tsv ];then
        break
    else
        cmd="vamb -e 100 -q 20 40 60 80 --outdir vamb.o --fasta ${pref_name}.tmp.fa --jgi ${depth} -t ${batch_size}"
        echo cd ${pref} > ${BASH_LOG}
        echo $cmd >> ${BASH_LOG}
        { $cmd  2> ${ERROR_LOG}.tmp && rm ${ERROR_LOG}.tmp && touch vamb.ok && break ;} || rm -r vamb.o
    fi

    if [ ${batch_size} -lt 2 ]; then
        echo -e "ERROR:\tbatch size ERROR(too small):${batch_size}" >> ${ERROR_LOG}.tmp && mv ${ERROR_LOG}.tmp ${ERROR_LOG}  && exit 127
    fi
    batch_size=$(( batch_size/2 ))
done

# extract bins
if [ -f vamb.o/clusters.tsv ];then
    ( [ ! -d ${pref}/bins ] ) && { mkdir -p ${pref}/bins || exit 127 ;}

    # 如果使用了seqkit grep命令，有可能出现找不到序列的情况，而且不会报错
    perl -e '%x;%y;%len;open I, "vamb.o/clusters.tsv";open L,"'"${pref_name}.tmp.fa.len"'";while(<L>){chomp;@l=split/\t/;$len{$l[0]}=$l[1]};while(<I>){chomp;$_=~/^(\d+)\s+(\S+)/;push @{$x{$1}},$2 ; $y{$1}+=$len{$2}}; foreach $k(keys %y){if ($y{$k} > '"$minBin"'){open O,">bins/$k.fa.list";@i=@{$x{$k}}; $j=join("\n",@i); print O $j} }' \
        && { ls bins/*.list | parallel -j 10 $get_fasta {} ${pref_name}.tmp.fa -e \> {.} 2\> ${pref_name}.name.err && touch extra.ok; }

    if [ -f extra.ok ];then
        if [[ `stat -c "%s" ${pref_name}.name.err` -ne 0 ]];then
            echo -e "\n\n\nextra err" >> ${pref}.log && exit 127
        else 
            rm -r vamb.ok extra.ok ${pref_name}.* \
                && mv vamb.o/log.txt ${pref}.log \
                && mv bins/*.fa ./ \
                && rm -r bins vamb.o \
                && touch ${pref}.ok
        fi
    fi
fi
