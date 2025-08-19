## OAPGC -  Oral and Airway Prokaryotic Genome Catalogue

OAPGC was constructed by integrating 12,616 public metagenomic samples, 933 newly sequenced airway metagenomes, and 64,169 isolate genomes from relevant habitats. 99,215 non-redundant strain-level genomes were reconstructed, which were subsequently clustered into 2,474 species-level clusters. In addition, by incorporating Global gut Microbiome Reference (GMR) collection, we further extended our analysis to explore the distribution and potential roles of oral microbes within the gut environment.</font>

<br>
The customized Kraken2 and Bracken database of OAPGC , integrated with GMR, is available at https://zenodo.org/records/16880410
<br>

------

### Single-coverage and Mash-based multiple-coverage binning

To maximize the recovery of microbial genomes and improve genome quality, metagenome-assembled genomes in the OAPGC were obtained through four single-coverage binning methods (MetaBAT, MetaBinner, Semibin, and VAMB) along with a Mash-based multiple-coverage binning  approach based on MetaBAT 2 (namely MetaBAT2-multi). Details of depth files generation and MetaBAT2-multi can be found in our previous study,  VMGC ([VMGC/README.md at main · RChGO/VMGC · GitHub](https://github.com/RChGO/VMGC/blob/main/README.md)). The specific operational steps about four single-coverage binning methods are outlined below.

- Required dependencies

  > GNU parallel 20201122<br>
  > MetaBAT2 v2.15<br>
  > MetaBinner v1.4.4<br>
  > SemiBin2 v2.0.2<br>
  > VAMB v3.0.2<br>
  > [flow_bin_metabinner.sh](https://github.com/NiYawen/OAPGC/blob/main/pipelines/flow_bin_metabinner.sh)<br>
  > [flow_bin_semibin2.sh](https://github.com/NiYawen/OAPGC/blob/main/pipelines/flow_bin_semibin2.sh)<br>
  > [flow_bin_vamb.sh](https://github.com/NiYawen/OAPGC/blob/main/pipelines/flow_bin_vamb.sh)<br>

- Step 1: MetaBAT 2

  > ls 01.depth/*.depth| parallel  -j 10 --plus metabat2 -i {.}.filter.fa.gz -a {} -o 01.metabat2_single/{/.}_Mbat2s -m 2000 -s 200000 --saveCls --seed 2020

- Step 2: MetaBinner

  > ls 01.depth/*.depth| parallel -j 10 --plus flow_bin_metabinner.sh {.}.filter.fa.gz {} 02.metabinner/{/.}

- Step 3:  Semibin

  > ls 01.depth/*.depth | parallel -j 10  --plus flow_bin_semibin2.sh human_oral {} {.}.filter.fa.gz 03.SemiBin/{/.} 2024 2000 200000 

- Step 4: VAMB

  > ls 01.depth/*.depth | parallel -j 10 --plus flow_bin_vamb.sh  {.}.filter.fa.gz  {}  04.vamb/{/.} 
  >
  > <br>

------

### Taxonomic profiling

Based on the 2474 species-level genome bins (SGBs) in the OAPGC, we reconstructed the prokaryotic composition of the oral and airway samples using Kraken2 and Bracken tools.

- Required dependencies

  > Kraken 2.1.3<br>
  > Bracken 2.8<br>
  > Python 3.8.16<br>
  > GNU parallel 20201122<br>
  > [kraken_build.py](https://github.com/NiYawen/OAPGC/blob/main/pipelines/kraken_build.py)<br>

- Step 1: Create customized Kraken2 and Bracken databases, including  prokaryotic species from the OAPGC, all viral operational taxonomic units (vOTUs) from the OAVGC (https://github.com/RChGO/AVGC), and the human reference genome.

  > cat ../01.SGB/genomes.v1.map  > merged.genomes.path
  > cat ../03.vOTU/all.virus.path >> merged.genomes.path 
  > cat human.path >> merged.genomes.path 
  >
  > cat ../01.SGB/genomes.v1.taxonomy  > merged.genomes.taxonomy
  > cat ../03.vOTU/all.virus.taxonomy >> merged.genomes.taxonomy 
  > cat human.taxonomy >> merged.genomes.taxonomy 
  >
  > ./kraken_build.py -l merged.genomes.path -t merged.genomes.taxonomy -d kraken2_db -p 80 > oral.build.log 2> oral.build.err<br>

- Step 2: The clean reads from each sample are mapped to the database, generating compositions at various taxonomic levels.

  > find clean_reads/*.1.fq.gz | sed 's/.1.fq.gz//'| parallel -j 5 kraken2 --threads 10 --confidence 0.1 --db kraken2_db --report prof/{/}.report --report-minimizer-data --output prof/{/}.output {}.1.fq.gz {}.2.fq.gz

  > find prof/*.report | parallel -j 5 bracken -d kraken2_db -i {} -o {}.bracken -r 150 -l S -t 1

<br>

*Note that the associated scripts nested within the shell scripts are also available on [pipelines](https://github.com/NiYawen/OAPGC/tree/main/pipelines).

<br>

------

<br>

<font size=2> Correspondence and requests for materials should be addressed to nywdoctor@163.com </font>


