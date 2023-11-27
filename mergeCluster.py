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
#
from collections import Counter


# calculate number of sample for bed row based on the name filed 
def calSamNum (name_str):
    saminf_list = name_str.split(',')
    sam_dict = {}
    for _sam in saminf_list:
        _sam_id,_ = _sam.split(':')
        sam_dict[_sam_id] = 1
    sam_num = len(sam_dict.keys())
    return sam_num
#

  
# refine the peak of the merged cluster based on occruency of cluster across multiple files 
def refinePeak (name_str,bed12_df):
    each_df = bed12_df[bed12_df['name'] == name_str].copy()
    pos = each_df.iloc[:,8].copy()
    # find the most frequenct position as refine PAS 
    counts = Counter(pos.to_list())
    _mf_pos = max(counts,key = counts.get)
    _start = _mf_pos - 1
    _cs = min(pos.to_list())
    _ce = max(pos.to_list())
    saminf_list = each_df.iloc[0,3].split(',')
    sam_dict = {}
    for _each in saminf_list:
        _k,_v = _each.split(':')
        sam_dict[_k] = _v
    _sam_id = ':'.join(sam_dict.keys())
    _name = _sam_id + '|' + each_df.iloc[0,0] + ':' + str(_cs) + '-' + str(_ce)
    #cl_da = {'chrom':each_df['chrom'][0],'start':_start,'end':_mf_pos,'name':_name,'score':each_df['score'][0],'strand':each_df['strand'][0]}
    #peak_df = pd.DataFrame(cl_da, index = ['chrom', 'start', 'end', 'name', 'score','strand'])
    peak = [each_df.iloc[0,0],_start,_mf_pos,_name,each_df.iloc[0,4],each_df.iloc[0,5]]
    return peak
    
    
# merge cluster from multiple samples, the cluster file from each samples should be seperated with "," 
def mergeCluster (cluster_tsv_files,output_cluster,read_cut,sam_num_cut,dis_cut,header,sam_lab):
    file_list = cluster_tsv_files.split(',')
    print(f'### merging cluster from {len(file_list)} files')
    if sam_lab != None:
        lab_list = sam_lab.split(',')
        if len(file_list) != len(lab_list):
            print(f'### warning number of files {len(file_list)} and number of labels {len(lab_list)} does not match')
            lab_list = None
    else:
        lab_list = None
    #
    i = 1
    df_list = []
    for _file in file_list:
        c_df = pd.read_csv(_file,sep = '\t',header = header)
        if lab_list != None:
            each_sam = lab_list[i - 1]
        else:
            each_sam = 'S' + str(i)
        # check if the file have head 
        if c_df.columns[3] != 'name':
            c_df.iloc[:,3] = each_sam + ':' + c_df.iloc[:,3]
            if read_cut != None:
                c_df = c_df.loc[c_df.iloc[:,4] > read_cut].copy()
        else:
            c_df['name'] = each_sam + ':' + c_df['name']
            if read_cut != None:
                c_df = c_df.loc[c_df['score'] > read_cut].copy()
        c_df = c_df.iloc[:,0:6].copy()
        df_list.append(c_df)
        i += 1
    union_df = pd.concat(df_list)
    union_bed = BedTool.from_dataframe(union_df)
    union_bed = union_bed.sort()
    merge_bed = union_bed.merge(c = '4,5,6',d = dis_cut, o = 'collapse,sum,distinct',s = True)
    merge_df = merge_bed.to_dataframe()
    sam_nums = merge_df['name'].apply(calSamNum)
    filtered_df = merge_df[sam_nums > sam_num_cut].copy()
    filtered_bed = BedTool.from_dataframe(filtered_df)
    # intesect filtered_df with union_df
    bed_wo = filtered_bed.intersect(union_bed,s = True,wa = True, wb = True)
    wo_df = bed_wo.to_dataframe()
    # to speed up, process the bed sperated by chrom 
    chrs = wo_df['chrom'].unique()
    output_list = []
    for _chr in chrs:
        chr_wo_df = wo_df[wo_df['chrom'] == _chr].copy()
        chr_cuid = chr_wo_df['name'].unique()
        print(f'### processing {len(chr_cuid)} cluster from chromosome {_chr}')
        peak_list = list(map(lambda _c: refinePeak(_c,chr_wo_df), chr_cuid))
        output_list += peak_list
        #break
    if output_cluster != None:
        #filtered_df.to_csv(output_cluster, index = False, sep = '\t')
        with open (output_cluster,'w') as w:
            for cl in output_list:
                w.write('\t'.join(str(_e) for _e in cl) + '\n')
    return filtered_df

if len(sys.argv) < 4:
    print("Usage:python mergeCluster.py [input_file1,input_file2,input_file3...] [output_file] [read_cut] [sample_number_cut] [distance_cut] [sample_label]")
    sys.exit(1)
else:
    if sys.argv[6] != 'none':
        header = 'infer'
    else:
        header = None
    sam_lab = None
    read_cut = int(sys.argv[3])
    sam_num_cut = int(sys.argv[4])
    dis_cut = int(sys.argv[5])
    if len(sys.argv) > 7:
        sam_lab = sys.argv[7]        
    mergeCluster(sys.argv[1],sys.argv[2],read_cut,sam_num_cut,dis_cut,header,sam_lab)
    
    
