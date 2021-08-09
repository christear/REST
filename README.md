# R3ESeq

A toolkit for RNA 3'End Sequencing data analysis.

#### 0. Prerequisites

- Python 3.6
- DNABERT (https://github.com/jerryji1993/DNABERT)
- scikit-learn
- pandas
- numpy
- bedtools
- Perl 5
- R 4.1
- bash
- umi_tools (only required for single cell sequening data analysis)

#### 1. Introduction

The R3ESeq is a toolkit for RNA 3'End Sequencing data analysis, which could be applied on both bulk RNA-seq data from mutiple 3'End Sequencing methods and common Single cell sequencing (SC-seq) data from 10X genomics Chromium System.

The tools has been tested on:

1, 3'End Sequencing data from Lexogen QuanSeq 3' mRNA-Seq Library Prep Kit REV (https://www.lexogen.com/quantseq-3mrna-sequencing-rev/)

2, SC-seq data from 10X genomics Chomium System (https://www.10xgenomics.com/instruments/chromium-x-series)

#### 2. Getting Started

The R3ESeq includes several key steps. 

1. call cluster based on RNA 3'End Sequencing data, required input should be in bam or bed format. The method was derived from the previous publications [1-2]. 
	- R3ESeq callcluster -i [input.file] -f [bam/bed] -o [output.file] 
2. filter reliable clusters as putative polyadenylation sites (PASs)
	a. using the pre-built deep-learning model 
	b. training the owning model using annotation 
3. count the sequencing reads for each putative PAS 
	a. bulk RNA-seq data 
	b. SC-seq data 
4. APA analyis
 
#### 3. Synopsis
#### 4. Tutorial
#### 5. The output




