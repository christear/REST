#!/usr/bin/env bash
:<<USE
download prerequisites to run R3ESeq toolkit, including 
1, deep-learning models 
2, Python modules
3, bedtools  
4, featureCounts
5, annotations  
...
USE

function usage () {
        echo
        echo
        echo Usage: species annotation work_directory 
        echo data.prepapration.sh [species] [annotation]
        echo supported species: human or mouse
        echo supported annotation resource: gencode, polyA database 
        echo
        echo
}

if [ -z $1 ]; then
        usage
        exit
fi

species=${1}
annotation=${2}
cd ${3}
#echo $species $annotation ${3}
##

# deep-learning models 
echo checking models ...
if [ ! -f DNABERT.model ]; then
	echo downloading DNA-BERT model ...
	#wget ... 
else
	echo DNABERT.model ... OK ...
fi

echo setup DNABERT model enviroment 
python ./DNABERT/setup.py install
conda create -n R3ESeq python=3.6
conda activate R3ESeq
conda install pandas
conda install tqdm
conda install --file requirements.txt 
conda install pytorch torchvision torchaudio -c pytorch


if [[ $species =~ "human" ]]; then
	echo downloading pre-built deep-learning model 'for' human ...
	# wget ...
elif [[ $species =~ "mouse" ]]; then
	echo downloading pre-built deep-learning model 'for' mouse ...
	# wget ...
else
	echo $species are not supported with current version 
fi

# python modules 
echo checking python modules ...
echo checking scikit-learn ... 
echo checking numpy ...
echo checking pandas ...
echo checking pysam ...
echo checking portion ... 
echo checking cigar ...
echo checking progressbar ...

# pip install pysam 
# pip install portion 
# pip install cigar
# pip install progressbar 

# required tools 
echo checking bedtools ... 
echo checking featureCounts ... 

if [[ $annotation =~ 'gencode' ]]; then
	echo downloading annotation from Gencode ...
	#wget ...
elif [[ $annotation =~ 'PolyA' ]]; then
	echo downloading annotation from polyA database ...
	#wget 
elif [[ $annotation =~ 'both' ]]; then
	echo downloading annotation from Gencode and polyA database ...
	# wget 
	# wget 
else
	echo no annotation will be downloaded ...
fi

 


