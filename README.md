# REST

REST: a RNA 3'End Sequencing data analysis Toolkit. 

#### 0. Prerequisites

- Python 3.6
- bedtools
- R
- bash
- featureCounts
- DNABERT (https://github.com/jerryji1993/DNABERT)

To install DNABERT, please use the dnabert.requirements.txt rather than the one under DNABERT/examples
`conda create -n rest python=3.6` \
`conda activate rest` \
`conda install pytorch torchvision cudatoolkit=10.0 -c pytorch` \
to run with GPU, please use this for pytorch
`conda install pytorch torchvision torchaudio pytorch-cuda=11.8 -c pytorch -c nvidia`
`git clone https://github.com/jerryji1993/DNABERT.git` \
`cd DNABERT` \
`python3 -m pip install --editable . ` \
`python3 -m pip install -r dnabert.requirements.txt`
#### 1. Introduction

The REST is a toolkit for RNA 3'End Sequencing data analysis, which could be applied on both bulk RNA-seq data from mutiple 3'End Sequencing methods and common Single cell sequencing (SC-seq) data, including 10X genomics Chromium System.

The tools has been tested on:

1, 3'End Sequencing data from Lexogen QuanSeq 3' mRNA-Seq Library Prep Kit REV (https://www.lexogen.com/quantseq-3mrna-sequencing-rev/)

2, SC-seq data from 10X genomics Chomium System (https://www.10xgenomics.com/instruments/chromium-x-series)

#### 2. Getting Started

To run REST, plase install DNABERT firstly. Ideally, DNABERT should be located under the main directory of REST 

The REST includes several key steps (under developing/debuging). 

1. call cluster based on RNA 3'End Sequencing data, required input should be in bam or bed format. The method was derived from the previous publications [1-2]. 
	- `python rest.py callcluster --input_file [input.file] --input_format [bam/bed] --strand 2 --output [output.cluster] --output_dis [cluster.distance]`
2. filter reliable clusters as putative polyadenylation sites (PASs)
	
	a. retrain/fune-tuen a pre-trained DNABERT model to the PASBERT model 
	- `python rest.py filtercluster --run train --DNABERT_path [DNABERT.path] --input_file [input.cluster] --out_dir [output.dir] --reference [reference.genome] --annotation [pas.annotation] --model [pre-trained.DNABERT.model] --round 5 --kmer 5 --motif_file human.pas.motif`
	
	b. predict based on the pre-tained PASBERT model  
	- `python rest.py filtercluster --run pred --DNABERT_path [DNABERT.path] --input_file [input.cluster] --out_dir [output.dir]  --model [model] --reference [reference.genome]`
		
3. count the sequencing reads for each putative PAS 
	
	a. bulk RNA-seq data 
	- `REST countbulk -bam [bam.file] -i [input.cluster] -o [output.count] -s [strandness] --win [width]`
	
	b. SC-seq data 
	- `REST countsc -bam [bam.file] -i [input.cluster] -o [output.count] -s [strandness] --win [width] -f [output.type/count or usage]`
4. APA analyis

	- `REST apa -i [input.count.tab] -c [condition.file] -o [output.apa.events] -a [condition.a] -b [condition.b]`

#### 3. Synopsis
#### 4. Tutorial
#### 5. The output




