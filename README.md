# REST

REST: a RNA 3'End Sequencing data analysis Toolkit. 

![REST](https://github.com/christear/REST/blob/main/images/REST.jpg)

#### 0. Prerequisites

- Python >= 3.6
- BEDTOOLS
- DNABERT (https://github.com/jerryji1993/DNABERT)
- featureCounts (only required for count)

Download the REST with GitHub or manually \
`git clone https://github.com/christear/REST.git`

To run REST, plese install DNABERT firstly with the below codes, indeally under the directory of REST \
`conda create -n rest python=3.6` \
`conda activate rest` \
`conda install pytorch torchvision cudatoolkit=10.0 -c pytorch` \
if the conda install takes too long time, please try \
`pip install torch` it works with torch=1.10.2=pypi_0

to run with GPU, please use this for pytorch \
`conda install pytorch torchvision torchaudio pytorch-cuda=11.8 -c pytorch -c nvidia` 

Download DNABERT with GitHub or manually \
`git clone https://github.com/jerryji1993/DNABERT.git` \
`cd DNABERT` \
`python3 -m pip install --editable . ` or `pip install --editable . --ignore-intsalled certifi` if certifi got errors\
please use the dnabert.requirements.txt rather than the one under DNABERT/examples \
`python3 -m pip install -r ../dnabert.requirements.txt` or `pip install -r ../dnabert.requirements.txt` \
`cd ..`

After sucessfuly install of DNABERT, install other required packages for REST \
`python3 -m pip install -r requirements.txt` or `pip install -r requirements.txt`

To retrain the model (PASBERT) for polyadenylation site (PAS) analysis , please download the initial pretrained DNABERT model as mentioned in (https://github.com/jerryji1993/DNABERT)

To directly predict/filter true PAS, please download our finetuned model from Google Driver  (https://drive.google.com/drive/folders/1Ns28TP24QpTPqTX-JeYSd92iNG-ZTS2q?usp=share_link)

To run count functions, please install featureCounts \
`conda install -c bioconda subread`

#### 1. Introduction
The REST is a toolkit for RNA 3'End Sequencing data analysis, which could be applied on both bulk RNA-seq data from mutiple 3'End Sequencing methods and common Single cell sequencing (SC-seq) data, including 10X genomics Chromium System.

The tools has been tested on:

1, 3'End Sequencing data from Lexogen QuanSeq 3' mRNA-Seq Library Prep Kit REV (https://www.lexogen.com/quantseq-3mrna-sequencing-rev/)

2, SC-seq data from 10X genomics Chomium System (https://www.10xgenomics.com/instruments/chromium-x-series)

#### 2. Getting Started
The REST includes several key functions (some are still under development/debugging). 

1. call_cluster: call cluster based on 3' end-seq data using the method from previous studies[1,2]. 

2. filter_cluster: filter the cluster based on sequence features fine-tuned by the DNABERT model. It has two options, including 'retrain' or 'pred'. 
- The 'retrain' option will fine-tune a pre-trained DNABERT model based on the identified cluster and annotation to define true/false PAS as well as an active-learning-based strategy to re-label true/false for each round of prediction. Finally, A fine-tuned model (PASBERT) will be returned. 
- The ‘pred’ option will predict the true/false PAS for each cluster based on a fine-tuned model (PASBERT). 

3. merge_cluster: merge the filtered cluster from multiple samples and refine the peak based on the most frequent position across multiple samples. 

4. annotate_PAS: annotate PAS with a gene annotation GTF file as described in a previous study[3].

5. count_PAS: count the number of supporting reads for each PAS from bam file. The bam file could be bulk 3\'-end sequencing data or scRNA-seq data. 

6. count_gene: count the read number for each gene by aggregating counts of PAS from the corresponding gene 
 
6. calculate PAS usage (PAU) or weighted 3' UTR length index (WULI) as mentioned in the previous study[3].
#### 3. Tutorial
The details of each function are shown below. 

1. call cluster based on RNA 3'End Sequencing data, required input should be in bam or bed format. The method was derived from the previous publications [1-2]. 
	- `python rest.py call_cluster --input_file [input.file] --input_format [bam/bed] --strand 2 --output [output.cluster] --output_dis [cluster.distance]`. The --strand option stands for the strand information of the sequencing data. 1: forward strand; 2: reverse strand.

2. filter reliable clusters as putative polyadenylation sites (PASs)
	
	a. retrain/fune-tune a pre-trained DNABERT model to the PASBERT model 
	- `python rest.py filter_cluster --run train --DNABERT_path [DNABERT.path] --input_file [input.cluster] --out_dir [output.dir] --reference [reference.genome] --annotation [pas.annotation] --model [pre-trained.DNABERT.model] --round 5 --kmer 5 --motif_file human.pas.motif`
	
	b. predict based on the pre-trained PASBERT model  
	- `python rest.py filter_cluster --run pred --DNABERT_path [DNABERT.path] --input_file [input.cluster] --out_dir [output.dir]  --model [model] --reference [reference.genome]`

3. merge the cluster passed PASBERT filtering from multiple samples
	- `python rest.py merge_cluster --file_list [file1.txt,file2.txt,file3.txt...] --read 5 --sam_num 2 --distance 25 --output [merged.pas.txt]`

4. annotate the putative PAS based on gene annotation in GTF format[3]. 
	- `python rest.py annotate_PAS --input_pas pas_file --gtf annotation_gtf --output output --header`
	
5. count the sequencing reads for each putative PAS

	a. bulk RNA-seq data 
	- `python rest.py count_PAS --bam [bam.file] --data_type bulk --pas [pas.file] --output [output.count] --strand [0/1/2] --win [24] --end [5/3] --thread 20`
	
	b. SC-seq data (under development/debug)
	- `python rest.py count_PAS --bam [bam.file] --data_type sc --pas [pas.file] --output [output.count] --strand [0/1/2] --win [24] --end [5/3] --thread 20`

6. count the reads for each gene by aggregating count of PAS
	- `python rest.py count_gene --pas_count [pas.count.file] --pas_anno [pas.anno.file] --output [output.gene.count] --pas_type [last_exon,intron_anno,intron_unanno,ups_exon,ncRNA]`

7. APA analyis (under development)
	
	a. calculate PAS usage (PAU) for each PAS 
	- `python rest.py cal_pau -i [input.count.tab] -c [condition.file] -o [output.apa.events] -a [condition.a] -b [condition.b]`
	
	a. calculate weighted 3' UTR length index (WULI) as described in previous study[3].  
	- `python rest.py cal_wuli -i [input.count.tab] -c [condition.file] -o [output.apa.events] -a [condition.a] -b [condition.b]`

#### 4. The output
Output of each functions usally in a bed-like format, except that the 'retrain' option in 'filter_cluster' will output a directory contains the fine-tuned DNABERT model (PASBERT), inlcuding several files named as:\
config.json, pytorch_model.bin, special_tokens_map.json, tokenizer_config.json, training_args.bin, vocab.txt.  

#### 5. Reference
1. Hoque, M., Ji, Z., Zheng, D., Luo, W., Li, W., You, B., ... & Tian, B. (2013). Analysis of alternative cleavage and polyadenylation by 3′ region extraction and deep sequencing. Nature methods, 10(2), 133-139. 
2. Xiao, M. S., Zhang, B., Li, Y. S., Gao, Q., Sun, W., & Chen, W. (2016). Global analysis of regulatory divergence in the evolution of mouse alternative polyadenylation. Molecular systems biology, 12(12), 890. 
3. Long, Y., Zhang, B., Tian, S., Chan, J. J., Zhou, J., Li, Z., ... & Gao, X. (2023). Accurate transcriptome-wide identification and quantification of alternative polyadenylation from RNA-seq data with APAIQ. Genome Research, 33(4), 644-657.




