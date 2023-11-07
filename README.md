# REST

REST: a RNA 3'End Sequencing data analysis Toolkit. 

#### 0. Prerequisites

- Python 3.6
- DNABERT (https://github.com/jerryji1993/DNABERT)
- scikit-learn
- pandas
- numpy
- pysam
- portion
- cigar
- progressbar
- bedtools
- Perl 5
- R 4.1
- bash
- featureCounts

#### 1. Introduction

The REST is a toolkit for RNA 3'End Sequencing data analysis, which could be applied on both bulk RNA-seq data from mutiple 3'End Sequencing methods and common Single cell sequencing (SC-seq) data, including 10X genomics Chromium System.

The tools has been tested on:

1, 3'End Sequencing data from Lexogen QuanSeq 3' mRNA-Seq Library Prep Kit REV (https://www.lexogen.com/quantseq-3mrna-sequencing-rev/)

2, SC-seq data from 10X genomics Chomium System (https://www.10xgenomics.com/instruments/chromium-x-series)

#### 2. Getting Started

To run REST, plase run the data.preparation.sh to download prerequisites at first.

The REST includes several key steps. 

1. call cluster based on RNA 3'End Sequencing data, required input should be in bam or bed format. The method was derived from the previous publications [1-2]. 
	- REST callcluster -i [input.file] -f [bam/bed] -o [output.file] 
2. filter reliable clusters as putative polyadenylation sites (PASs)
	
	a. using the pre-trained deep-learning model 
	- REST filtercluster -i [input.file] -o [output.file] --model [pre-trained model]
	
	b. training the customized model using annotation 
	- REST filtercluster -i [input.file] -o [output.file] --outmodel [output.model]
		
3. count the sequencing reads for each putative PAS 
	
	a. bulk RNA-seq data 
	- REST countbulk -bam [bam.file] -i [input.cluster] -o [output.count] -s [strandness] --win [width]
	
	b. SC-seq data 
	- REST countsc -bam [bam.file] -i [input.cluster] -o [output.count] -s [strandness] --win [width] -f [output.type/count or usage]
4. APA analyis

	- REST apa -i [input.count.tab] -c [condition.file] -o [output.apa.events] -a [condition.a] -b [condition.b]

#### 3. Synopsis
#### 4. Tutorial
#### 5. The output




