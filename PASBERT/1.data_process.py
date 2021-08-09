# Data Process

import random
from sklearn.model_selection import train_test_split
from DNABERT.motif.motif_utils import seq2kmer
import pandas as pd
import os
import argparse

parser = argparse.ArgumentParser(description='data_process')
parser.add_argument('--root',help='eg: /sample_L1')
parser.add_argument('--round',type=int,help='eg: 0')
parser.add_argument('--mer',type=int,help='eg: 5')
args = parser.parse_args()
root=args.root
round=args.round
mer=args.mer

if True:
    data_path='{}'.format(root)
    output_path='{}/active_r{}_m{}-0'.format(data_path,round,mer)
    os.system('mkdir -p {}'.format(output_path))
    if round==0:
        os.system('cp {}/positive.events.txt {}'.format(data_path,output_path))
        os.system('cp {}/negative.events.txt {}'.format(data_path,output_path))
    _positive=open('{}/positive.events.txt'.format(output_path),'r')
    _negative=open('{}/negative.events.txt'.format(output_path),'r')
    _positives=[]
    _negatives=[]
    for line in _positive.readlines():
        line=data_path+'/data/'+line.rstrip('\n')+'.data.txt'
        _positives.append(line)
    for line in _negative.readlines():
        line=data_path+'/data/'+line.rstrip('\n')+'.data.txt'
        _negatives.append(line)
    X=_positives+_negatives
    Y=[1 for i in range(len(_positives))]+[0 for i in range(len(_negatives))]
    #X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0)
    
    print('finished split train and test')

    def read_sequence(path):
        line=open(path).readlines()[0].split('\t')[0]
        return line
    
    print('start process tsv')

    # train.tsv
    with open('{}/train.tsv'.format(output_path),'w') as w:
        w.write('sequence\tlabel\tpath\n')
        for i in range(len(X)):
            w.write('{}\t{}\t{}\n'.format(seq2kmer(read_sequence(X[i]),5),Y[i],X[i]))
    print('train.tsv saved')

