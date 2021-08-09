# Data Process

import random
import numpy as np
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
data=pd.read_csv('{}/pred/dev.tsv'.format(data_path),sep='\t')
_pred=np.load('{}/pred/pred_results.npy'.format(data_path))
threshold=0.5
_pred01=np.int8(_pred>=threshold)
data['pred_value']=_pred
data['pred_01']=_pred01
data.to_csv('{}/pred/dev_active_r{}_m{}.tsv'.format(data_path,round,mer),index=False,sep='\t')
