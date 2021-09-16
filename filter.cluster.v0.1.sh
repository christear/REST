#!/usr/bin/env bash
:<<USE
filter the cluster by DNA-BERT:
a. using the pre-trained deep-learning model
b. training the customized model using annotation
## input
1, option [training/prediction] 
2, identified cluster from 3'end sequencing data 
3, genome reference 
4, model / could be the pre-built model or DNA-BERT initial model for training/fine-tune
5, output        
6, PAS annotation, only required for training 
7, round of training, only required fot training 
8, output model, only required for training 
... 
USE

option=${1}
cluster=${2}
genome=${3}
model=${4}
output=${5}

# 
if [ $option =~ 'train' ]; then
	anno=${6}
	n=${7}
	outmodel=${8}
	# dataset preparation
	echo dataset preparation ...
	grep "chr[1-9XY]" $cluster | awk '{OFS="\t"}{print $1,$8 - 24,$8 + 24,$4,$5,$6}' > $cluster.win48
	intersectBed -a $cluster.win48 -b $anno -s -u | awk '{OFS="\t"}{print $1,$2 - 76,$3 + 76,$4,$5,$6}'| fastaFromBed -bed -fi $genome -s -name -tab -fo $cluster.anno.txt
	intersectBed -a $cluster.win48 -b $anno -s -v | awk '{OFS="\t"}{print $1,$2 - 76,$3 + 76,$4,$5,$6}'| fastaFromBed -bed -fi $genome -s -name -tab -fo $cluster.unanno.txt
	awk '{print $0,1}' $cluster.anno.txt > $cluster.data.txt
	awk '{print $0,0}' $cluster.unanno.txt >> $cluster.data.txt
	# 
	echo training the model ...
	xx.py/xx.sh $cluster.data.txt $model $output $n $outmodel   
	# 
elif [ $option =~ 'pred']; then
	echo runing the prediction ...
	grep "chr[1-9XY]" $cluster | awk '{OFS="\t"}{print $1,$2 - 76,$3 + 76,$4,$5,$6}'| fastaFromBed -bed -fi $genome -s -name -tab -fo $cluster.data.txt 
	#
	echo running the prediction
	xx.py/xx.sh $cluster.data.txt $model $out  
	#
else
	echo running option: $option was not defined
	exit 
fi

