#!/usr/bin/python3

import tempfile, functools
import pandas as pd
import argparse
import pickle
import gzip, bz2, re
from multiprocessing import Pool
import itertools
import os
import logging
import subprocess
import asyncio
import sqlite3

__doc__ = f'''
    genomes list:
        GUT_GENOME000001        /share/data1/lvqb/immune/01.data/UHGG.fasta.data/MGYG-HGUT-00001.fna.gz
        GUT_GENOME000004        /share/data1/lvqb/immune/01.data/UHGG.fasta.data/MGYG-HGUT-00002.fna.gz
        GUT_GENOME000008        /share/data1/lvqb/immune/01.data/UHGG.fasta.data/MGYG-HGUT-00003.fna.gz
        GUT_GENOME000010        /share/data1/lvqb/immune/01.data/UHGG.fasta.data/MGYG-HGUT-00004.fna.gz

    genomes taxo:
        ID                      superkingdom    phylum          class                   order   family  genus   species             subspecies
        GUT_GENOME118144        Bacteria        Proteobacteria  Alphaproteobacteria     RF32    CAG-239 51-20   GUT_GENOME118144    xxxx
        GUT_GENOME244811        Bacteria        Proteobacteria  Alphaproteobacteria     RF32    CAG-239 51-20   GUT_GENOME118144    asdf
        GUT_GENOME219382        Bacteria        Proteobacteria  Alphaproteobacteria     RF32    CAG-239 51-20   GUT_GENOME219382    
        GUT_GENOME245766        Bacteria        Proteobacteria  Alphaproteobacteria     RF32    CAG-239 51-20   GUT_GENOME245766    sub.x
'''

# ANSI 转义序列定义颜色
class ColoredFormatter(logging.Formatter):
    # 定义颜色
    BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE = range(30, 38)

    RESET_SEQ = "\033[0m"
    COLOR_SEQ = "\033[1;%dm"
    # 定义日志级别与颜色的映射
    LEVEL_COLORS = {
        logging.DEBUG: CYAN,
        logging.INFO: GREEN,
        logging.WARNING: YELLOW,
        logging.ERROR: RED,
        logging.CRITICAL: MAGENTA,
    }
    
    def __init__(self, fmt=None, datefmt=None):
        super().__init__(fmt, datefmt)

    def format(self, record):
        # 获取日志级别对应的颜色
        color = self.COLOR_SEQ % (self.LEVEL_COLORS.get(record.levelno, self.WHITE))
        # 格式化信息
        levelname = record.levelname
        colored_levelname = f"[{color}{levelname}{self.RESET_SEQ}]"
        record.levelname = colored_levelname
        message = super().format(record)
        record.levelname = levelname
        return message




def get_args():
    readsLen = "150"
    parser = argparse.ArgumentParser(description = __doc__, formatter_class = argparse.RawTextHelpFormatter)
    parser.add_argument("-l", help="genomes list")
    parser.add_argument("-t", help="genomes taxo")
    parser.add_argument("-d", "--database", help="database path.")
    parser.add_argument("-p", "--threads", type=int, default=None, help="threads. [defaul: all]")
    parser.add_argument("-a", "--add", action = "store_true", help="add new genomes to the database.")
    parser.add_argument("-m", "--readsLen", default=readsLen, help=f"length of reads.[default:{readsLen}]")
    args = parser.parse_args()
    return(args)


def read_file(filepath):
    logging.debug(f"read file: filepath")
    with open(filepath, 'rb') as f:
        header = f.read(4)
        if header[:2] == b'\x1f\x8b':  # Gzip
            return gzip.open(filepath, 'rt')
        elif header[:2] == b'BZ':  # Bzip2
            return bz2.open(filepath, 'rt')
        else:        
            return open(filepath, 'r')


def load_data(filepath):
    logging.info(f"open file:\t{filepath}\n")
    try:
        with open(filepath, 'rb') as f:
            mydata = pickle.load(f)
            species_count = mydata['species_count']
            genome2taxid = mydata['genome2taxid']
            taxids = mydata['taxids']
            id_map = mydata['id_map']
            names_str = mydata['names_str']
            nodes_str = mydata['nodes_str']
        return (species_count, genome2taxid, taxids, id_map, names_str, nodes_str)
    except Exception as e:
        logger.error(f"{e}")
        exit(127)


def taxo_init():
    species_count = dict() # 计数每个species有多少个基因组， {"species1":count, "species2": count, ...}
    genome2taxid = dict() # {"genome_id": ("speiec taxid", order) } # 基因组属于哪个物种，以及这个基因组在物种当中是第几个基因组(这个顺序之和出现的顺序有关)
    levels_name = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species' , 'supspecies']
    levels_num = [1, 2, 3, 4, 5, 6, 7, 8]
    taxids = {x:1 for x in levels_num} # 用于计数每个层级有多少个分类
    id_map = {1: 'no rank'} # 每个taxo_name对应的taxid
    parent_id = 1 # 当前taxid对应的上级taxid
    taxid = 11 # 当前taxid，默认是no rank = 11
    names_str = "1\t|\tno rank\t|\t\t|\tscientific name\t|\n"
    nodes_str = "1\t|\t1\t|\tno rank\t|\n"
    return (species_count, genome2taxid, taxids, id_map, names_str, nodes_str, parent_id, taxid, levels_name, levels_num)

 
def parser_taxo(filename, title=True, db_path="", old_map=None):
    '''
        input:
            genome_file.tsv
        output:
            genome2taxid, data
            genome2taxid = dict() # {"genome_id": ("speiec taxid", order, 0) } # 最后一个值代表之前是否已经添加过了
    '''
    ## 初始化变量
    species_count, genome2taxid, taxids, id_map, names_str, nodes_str, parent_id, taxid, levels_name, levels_num = taxo_init()
    
    ## 如果是向已有的文件中添加新的基因组，就加载旧的文件,这一步会替换上面提到的变量，所以必须在初始化之后
    if old_map:
        species_count, genome2taxid, taxids, id_map, names_str, nodes_str = load_data(f"{db_path}/taxonomy/my.pkl")
    
    f = read_file(filename)
    f.readline() # 跳过tittle
    
    for line in f:
        # linef = line.strip().split() 
        linef = re.split("\t", line.strip()) # ['genome_id', 'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'supspecies']
        genome_id = linef[0]
        for index, level_name in enumerate(levels_name):
            try:
                level_num = levels_num[index]
                taxo_name = linef[level_num] # 当前基因组的分类名称
            
                # 如果为空，则跳过
                if taxo_name == "":
                    continue
                
                # 如果没有被记录在内，则计入names.dmp, nodes.dmp
                if not id_map.get(taxo_name):
                    taxid = f"{level_num}{taxids[level_num]}" # 相当于第几个分类层级的第几个物种
                    id_map[taxo_name] = f"{taxid}" # taxoname: taxo_id
                    nodes_str += f"{taxid}\t|\t{parent_id}\t|\t{level_name}\t|\n"
                    names_str += f"{taxid}\t|\t{taxo_name}\t|\t\t|\tscientific name\t|\n"
                    taxids[level_num] += 1
                parent_id = id_map[taxo_name] # 将当前节点的ID作为下一层级的父节点
            except Exception as e:
                pass
        # 如果
        if taxid == 11:
            continue
        ## 因为每次循环完以后，都是最小的分类层级，所以直接计数就可以
        if species_count.get(taxid):
            species_count[taxid] += 1
        else:
            species_count[taxid] = 1

        genome2taxid[genome_id] = [taxid, species_count[taxid], 0]
    f.close()
    
    data = {
        'taxids': taxids, 'id_map': id_map,
        'names_str': names_str, 'nodes_str': nodes_str,
        "species_count":species_count,
        "genome2taxid": genome2taxid
    }
    if not os.path.exists(f"{db_path}/taxonomy"):
        run_cmd(f"mkdir -p {db_path}/taxonomy")
    
    with open(f"{db_path}/taxonomy/names.dmp", 'w') as namefile,\
        open(f"{db_path}/taxonomy/nodes.dmp", 'w') as nodefile:
        namefile.write(names_str)
        nodefile.write(nodes_str)
        
    return genome2taxid, data


async def async_write_to_db(df, db_name):
    loop = asyncio.get_running_loop()
    await loop.run_in_executor(None, lambda: df.to_sql('results', sqlite3.connect(db_name), if_exists='append', index=False))
       
def process_genomes(args, db_path=None, temp_dir = None):
    filepath, genome_id = args
    species_id, _ = genome_id.split("C")
    db_sql = f"{db_path}/taxonomy/db.accession2taxid.sql"
    # gid = f"{species_id}{species_order}"
    f = read_file(filepath)
    res = ""
    ctgs = []
    index = 1
    for  line in f:
        if RE_SEQ_NAME.match(line):
            res += f">G{genome_id}CTG{index}\n"
            ctgs.append(f"G{genome_id}CTG{index}")
            index += 1
            continue
        res += line
    f.close()
    df = pd.DataFrame({'accession':ctgs, 'accession_version':ctgs, "taxid": species_id })
    asyncio.run(async_write_to_db(df, db_sql))
    fa_path = os.path.join(db_path, temp_dir, f"{genome_id}.fna")
   
    with open(fa_path, 'w') as f:
        f.write(res)
        
    cmd = f"kraken2-build --add-to-library {fa_path} --db {db_path} && rm {fa_path} "
    run_cmd(cmd)


def file_generator(filename, genome2taxid, data):
    '''
        renamed = taxid + count
        input:
            /path/to/input/genome.fa, {'genome1': (taxid, order_of_taxid)}
        return:
            /path/to/input/genome.fa, /path/to/output, renamed
    '''
    with read_file(filename) as f:
        for line in f:
            name, filepath = re.split("\t", line.strip())
            if data['genome2taxid'][name][2]:
                continue
            else:
                data['genome2taxid'][name][2] = 1
            taxid, order, _ = genome2taxid[name]
            yield (filepath, f"{taxid}C{order}")


def format_genomes(genomelist, genome2taxid, db_path, data, chunk_size=1000, threads=5):
    '''
        filename: 基因组列表
            genome1 /path/to/genome1.fa
            genome2 /path/to/genome2.fa.gz
            ...     ...
    '''
    
    with sqlite3.connect(f"{db_path}/taxonomy/db.accession2taxid.sql") as conn:
        conn.execute('CREATE TABLE IF NOT EXISTS results (accession TEXT, accession_version TEXT, taxid TEXT, gi INTEGER PRIMARY KEY AUTOINCREMENT)')
    new_genome = False
    with tempfile.TemporaryDirectory(dir=db_path) as temp_dir:
        # 创建进程池
        with Pool(processes=threads) as pool:
            # 使用生成器逐步处理文件
            tasks = file_generator(genomelist, genome2taxid, data)
            while True:
                # 从生成器中获取一批任务
                task_batch = list(itertools.islice(tasks, chunk_size))
                if not task_batch:
                    break
                # pool.map(rename_seq, task_batch)
                new_genome = True
                pool.map( functools.partial(process_genomes, db_path=db_path, temp_dir=temp_dir), task_batch)
    # 将db.accession2taxid.sql, 转为文件db.accession2taxid
    conn = sqlite3.connect(f"{db_path}/taxonomy/db.accession2taxid.sql")
    final_df = pd.read_sql_query("SELECT * FROM results", conn)
    conn.close()
    final_df.columns = ['accession','accession.version','taxid','gi']
    final_df.to_csv(f"{db_path}/taxonomy/db.accession2taxid", sep="\t", index=False)
    return new_genome


def run_cmd(cmd):
    logger.info(f"RUN: {cmd}")
    res = subprocess.Popen(cmd, shell=True,
           #  stdout=None,
            stderr=subprocess.PIPE,
            stdout=subprocess.PIPE,
            text = True
            )
    # 如果stderr/stdout被捕获，则需要执行下面的命令，否则就会造成管道阻塞
    _, stderr = res.communicate()
    if res.returncode != 0:
        logger.error(stderr)
        raise


def main(args):
    db_path = os.path.abspath(args.database)
    taxofile = args.t
    genome_list_file = args.l
    readsLen = args.readsLen
    threads = os.cpu_count() if args.threads is None else args.threads
    old_map = args.add
    # 1、解析分类文件
    genome2taxid, data = parser_taxo(taxofile, db_path=db_path, old_map=old_map)

    # 2.1、生成文件db.accession2taxid， names.dmp, nodes.dmp
    # 2.2、重命名序列名称 && kraken2-build --add-to-library
    new_genome = format_genomes(genome_list_file, genome2taxid, db_path, data, threads=threads)
    if not new_genome:
        logger.warning("No new genome add.")
        exit(0)
    with open(f"{db_path}/taxonomy/my.pkl", 'wb') as f:
        pickle.dump(data, f)
    # 3、 kraken2-build --build --db $DBNAME
    cmd = f"kraken2-build --build --db {db_path} --threads {threads}"
    run_cmd(cmd)
    
    # 4、建库bracken
    cmd = f"bracken-build -d {db_path} -t {threads} -k 35 -l {readsLen}"
    run_cmd(cmd)

if __name__ == "__main__":
    #logger = ColoredFormatter("%(asctime)s - %(levelname)s - %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
   # 创建日志记录器
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    # 禁用默认的日志处理器传播
    logger.propagate = False

    # 创建控制台处理器
    ch = logging.StreamHandler()
    # 设置日志输出格式
    formatter = ColoredFormatter("%(asctime)s - %(levelname)s - %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    
    RE_SEQ_NAME = re.compile("^>")
    args = get_args()
    main(args)
