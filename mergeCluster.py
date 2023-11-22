# merge cluster output from filtering step with DNABERT prediction
# python=3.6-3.10
# Bin Zhang
# Data: Nov 21, 2023
# version: 1.0
# ENV:

import os
import sys
import numpy as np
import pandas as pd
import pybedtools
from pybedtools import BedTool

# calculate number of sample for bed row based on the name filed 
def calSamNum (name_str):
    saminf_list = name_str.split(',')
    samdict = {}
    for _sam in saminf_list:
        _sam_id,_ = _sam.split(':')
        samdict[_sam_id] = 1
    sam_num = len(samdict.keys())
    return sam_num
    
# merge cluster from multiple samples, the cluster file from each samples should be seperated with "," 
def mergeCluster (cluster_tsv_files,output_cluster,read_cut,sam_num_cut,dis_cut):
    file_list = cluster_tsv_files.split(',')
    print(f'### Merging cluster from {len(file_list)} files')
    i = 1
    df_list = []
    for _file in file_list:
        c_df = pd.read_csv(_file,sep = '\t')
        # check if the file have head 
        if c_df.columns[3] != 'name':
            c_df.iloc[:,3] = 'S' + str(i) + ':' + c_df.iloc[:,3]
            if read_cut != None:
                c_df = c_df.loc[c_df.iloc[:,4] > read_cut].copy()
        else:
            c_df['name'] = 'S'+ str(i) + ':' + c_df['name']
            if read_cut != None:
                c_df = c_df.loc[c_df['score'] > read_cut].copy()
        df_list.append(c_df)
    union_df = pd.concat(df_list)
    union_bed = BedTool.from_dataframe(union_df)
    union_bed = union_bed.sort()
    merge_bed = union_bed.merge(c = '4,5,6',d = dis_cut, o = 'collapse,sum,distinct')
    merge_df = merge_bed.to_dataframe()
    sam_nums = merge_df['name'].apply(calSamNum)
    filtered_df = merge_df[sam_nums > sam_num_cut].copy()
    if output_cluster != None:
        filtered_df.to_csv(output_cluster, index = False, sep = '\t')
    return filtered_df

if len(sys.argv) < 4:
    print("Usage:python mergeCluster.py [input_file1,input_file2,input_file3...] [output_file] [read_cut] [sample_number_cut] [distance_cut]")
    sys.exit(1)
else:
    mergeCluster(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5])
    
    
