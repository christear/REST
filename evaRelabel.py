# countGene
# python=3.6-3.10
# Bin Zhang
# Data: Dec 10, 2023
# version: v1
# ENV: 
#

import re
import pandas as pd
import numpy as np
from utils.PAS_utils import kmer2seq, determineIP, determineMotif

# the input should be the relablled tsv and the pred tsv from the next round of prediction 
def evaRelabel(relabel_tsv,pred_tsv):
    df1 = pd.read_csv(relabel_tsv,sep = '\t')
    df2 = pd.read_csv(pred_tsv,sep = '\t')
    cseq_list = df1["sequence"].apply(kmer2seq)
    ## get index
    # PAS predicted as 1 
    indices1 = df1[df1['pred_01'] == 1].index 
    # PAS that have been relabed from 1 to 0 
    indices1_0 = df1[(df1['pred_01'] == 1) & (df1['relabel'] == 0)].index
    # PAS that have been relabed from 0 to 1
    indices0_1 = df1[(df1['pred_01'] == 0) & (df1['relabel'] == 1)].index
    # PAS predicted as 1 after relabelling-predicting
    indices2 = df2[df2['pred_01'] == 1].index
    indices_list = [indices1,indices1_0,indices0_1,indices2]
    # count the PAS motif and internap primming frequency
    motif_list = ['AATAAA','ATTAAA']
    pasfreq_list = []
    ipfreq_list = []
    for each_i in indices_list:
        pasfreq_list.append(cseq_list[each_i].apply(determineMotif,motif_list = motif_list,k = 2).sum())
        ipfreq_list.append(cseq_list[each_i].apply(determineIP,start = 85,end = 115).sum())
    # count internal primming frequency 
    out_df = pd.DataFrame({'pas_freq':pasfreq_list,'ip_freq':ipfreq_list,'num':list(map(len,indices_list))})
    return out_df

def evaPred(pred_tsv):
    df = pd.read_csv(relabel_tsv,sep = '\t')
    cseq_list = df["sequence"].apply(kmer2seq)
    indices1 = df1[df1['pred_01'] == 1].index 
    indices0 = df1[df1['pred_01'] == 0].index 
    indices_list = [indices1,indices0]
    #
    motif_list = ['AATAAA','ATTAAA']
    pasfreq_list = []
    ipfreq_list = []
    for each_i in indices_list:
        pasfreq_list.append(cseq_list[each_i].apply(determineMotif,motif_list = motif_list,k = 2).sum())
        ipfreq_list.append(cseq_list[each_i].apply(determineIP,start = 85,end = 115).sum())
    # count internal primming frequency 
    out_df = pd.DataFrame({'pas_freq':pasfreq_list,'ip_freq':ipfreq_list,'num':list(map(len,indices_list))})
    return out_df
      
    
    
    
    




'''
    relabel_tsv="~/c2092/myc_apa/04_singleCell_data/scpas_retrain1/ini_pred/pred.relabeled.tsv"
    pred_tsv="~/c2092/myc_apa/04_singleCell_data/scpas_retrain2/ini_pred/pred.tsv"
'''

