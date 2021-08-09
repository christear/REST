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

data_path='{}/active_r{}_m{}-0'.format(root,round,mer)
relabel_data=pd.read_csv('{}/dev_active_r{}_m{}.relabeled.tsv'.format(data_path,round-1,mer),sep='\t')

with open('{}/positive.events.txt'.format(data_path),'w') as p:
    with open('{}/negative.events.txt'.format(data_path),'w') as n:
        for i in range(len(relabel_data)):
            idx=relabel_data.loc[i,'path'].split('/')[-1].rstrip('.data.txt')
            label=relabel_data.loc[i,'relabel']
            if label==1:
                p.write('{}\n'.format(idx))
            elif label==0:
                n.write('{}\n'.format(idx))
