# callCluster
# python=3.6-3.10
# Bin Zhang
# Data: Sep 14, 2021
# Last modification: Nov, 14, 2023
# version: 1.0
# ENV: 
# @ pip install cigar
# @ pip install portion
# @ pip install pysam
# @ pip install progressbar

import os
import re
import sys
import cigar
import pysam
import datetime
import numpy as np
import pandas as pd
from multiprocessing import Pool
#import portion as P
#from collections import Counter
#from collections import defaultdict
#from progressbar import ProgressBar, Bar, ETA

cluster_width = 24
# input file format should be bam to bed format
def input2dict (input_file,input_format,if_stranded):
    #if_stranded = int(if_stranded)
    # to determin strand 
    strand_dict = {'+':'-','-':'+'}
    locus_covs = {'+' : {},'-' : {}, '?' : {}}
    if 'bed' in input_format:
        print('[INFO] The input: {} is in bed format'.format(input_file))
        with open(input_file,'r') as bed:
            for _line in bed:
                eles = re.split('[\t\ ]',_line.rstrip())
                _chr = eles[0]
                _strand = eles[5]
                _locus = eles[2]
                if if_stranded == 2:
                    _strand = strand_dict[_strand]
                if _strand == '-':
                    _locus = eles[1]
#                
                #print(_chr,_strand,_locus)
                try:
                    locus_covs[_strand][_chr][_locus] += 1
                except:
                    try:
                        locus_covs[_strand][_chr].setdefault(_locus,1)
                    except:
                        locus_covs[_strand].setdefault(_chr,{})
                        locus_covs[_strand][_chr].setdefault(_locus,1)
                                 
    elif 'bam' in input_format:
        print('[INFO] The input: {} is in bam format'.format(input_file))
        for read in pysam.AlignmentFile(input_file, "rb"):
            # skip secondary alignment and PCR duplicates 
            if (read.flag & 256) | (read.flag & 512) | (read.flag & 1024):
                continue
            _chr = read.reference_name
            if read.flag & 16:
                _strand = '-'
            else:
                _strand = '+'        
            _end = read.reference_start
            cigars = read.cigartuples
            for ec in cigars:
                if ec[0] == 0 or ec[0] == 2:
                    _end += ec[1]
                elif ec[0] == 3:
                    _end += ec[1]
                #elif ec[0] == 1:
                #    _end -= ec[1]
            _locus = _end
            if if_stranded == 2:
                _strand = strand_dict[_strand]
            if _strand == '-':
                _locus = read.reference_start
            #print("bed",_chr,read.reference_start,_end,read.query_name,'255',_strand)
            _locus = str(_locus)           
            try:
                locus_covs[_strand][_chr][_locus] += 1
            except:
                try:
                    locus_covs[_strand][_chr].setdefault(_locus,1)
                except:
                    locus_covs[_strand].setdefault(_chr,{})
                    locus_covs[_strand][_chr].setdefault(_locus,1)    
    else:
        print('[INFO] undetermined format {}'.format(input_format))
    return locus_covs
    #print(locus_covs['-']['chr1']['3020262'])
    #for _strand in locus_covs:
    #    for _chr in locus_covs[_strand]:
    #        for _locus in sorted (locus_covs[_strand][_chr]):
    #            print(_chr,_locus,_strand)

# convert dict key to int  
def key2int (unsorted_dict):
    int_key = []
    for k in unsorted_dict:
        int_key.append(int(k))
    return int_key

# find peak 
def findPeak (cov_dict):
    global cluster_width
    _n = 0
    _p = 0
    _max = 0
    for _i in cov_dict:
        if cov_dict[_i] > _max:
            _p = _i
            _max = cov_dict[_i]
    loci = key2int(cov_dict)
    _l = loci[-1]
    _r = loci[0]
    _peak = int(_p)
    for _locus in loci:
        if abs(_locus - _peak) < cluster_width + 1:
            _n += cov_dict[str(_locus)]
            if _locus < _l:
                _l = _locus
            if _locus > _r:
                _r = _locus
            cov_dict.pop(str(_locus))
    res = [cov_dict,[_l - 1,_r,_peak,_n]]
    return res

#"""
# call cluster from the dict 
def callCluser (input_file,input_format,if_stranded,out_cluster,out_dis):
    start_time = datetime.datetime.now()
    locus_covs = input2dict(input_file,input_format,if_stranded)
    _tc = 0
    _tc1 = 0
    cluster_dict = {'+' : {}, '-' : {}}
    singles_dict = {'+' : {}, '-' : {}}
    with open(out_dis, 'w') as w:
        _sep = '\t'
        for _strand in locus_covs:
            for _chr in locus_covs[_strand]:
                loci = key2int(locus_covs[_strand][_chr])
                _tmp = 0
                _start = 0
                ### use the mark to indicate if the before one loci is far than 24 ...  
                _marker = 1
                _last = 0
                cluster_dict[_strand].setdefault(_chr,{})
                singles_dict[_strand].setdefault(_chr,{})
                for _locus in sorted (loci):
                    _last = _locus
                    _num = locus_covs[_strand][_chr][str(_locus)]
                    _tc += _num
                    if _tmp == 0:
                        _tmp = _locus
                    else:
                        _d = _locus - _tmp
                        #_id = _chr + '|' + _strand + '|NA|NA|NA'
                        _dis_line = _sep.join(str(e) for e in [_chr,_tmp,_locus,_d]) + '\n'
                        w.write(_dis_line)
                        if _d < cluster_width + 1:
                            if _start == 0:
                                _start = _tmp
                            try:
                                cluster_dict[_strand][_chr][str(_start)].append(_locus)
                            except:
                                cluster_dict[_strand][_chr][str(_start)] = [_locus]
                            _marker = 0   
                        else:
                            _start = 0
                            if _marker == 0:
                                _marker = 1
                            else:
                                singles_dict[_strand][_chr][str(_tmp)] = locus_covs[_strand][_chr][str(_tmp)]
                                _tc1 += locus_covs[_strand][_chr][str(_tmp)]
                        _tmp = _locus
                if _start == 0:
                    singles_dict[_strand][_chr][str(_last)] = locus_covs[_strand][_chr][str(_last)]
                    _tc1 += locus_covs[_strand][_chr][str(_last)]    
        print('[INFO] there are {} reads in total'.format(_tc))
        w.close()
        #print('[INFO] there are {} reads in single locus'.format(_tc1))
    with open(out_cluster,'w') as w:
        _sep = '\t'
        print('[INFO] writing single locus cluster to {}'.format(out_cluster))
        for _strand in singles_dict:
            for _chr in singles_dict[_strand]:
                for _s in singles_dict[_strand][_chr]:
                    _start = int(_s)
                    _n = singles_dict[_strand][_chr][_s]
                    _id = '{}|{}|NA|{}|{}'.format(_chr,_strand,_s,_n)
                    _c_line = _sep.join(str(e) for e in [_chr,_start - 1,_start,_id,_n,_strand,_start - 1,_start]) + '\n'
                    w.write(_c_line)
        print('[INFO] {} reads are from single locus'.format(_tc1))
        _tc2 = 0
        print('[INFO] writing multiple locus cluster to: {}'.format(out_cluster))
        for _strand in cluster_dict:
            for _chr in cluster_dict[_strand]:
                for _s in cluster_dict[_strand][_chr]:
                    _start = int(_s)
                    _tc2 += locus_covs[_strand][_chr][_s]
                    if singles_dict[_strand][_chr].get(_s,-1) != -1:
                        print('[INFO] duplicated counting: {}'.format(_start))
                    for _e in cluster_dict[_strand][_chr][_s]:
                        if _e == _start:
                            print('[INFO] locus {} is duplicated')
                        _tc2 += locus_covs[_strand][_chr][str(_e)]
                        if singles_dict[_strand][_chr].get(str(_e),-1) != -1:
                            print('[INFO] duplicated counting: {}'.format(_e))
                    if cluster_dict[_strand][_chr][_s][-1] - _start < cluster_width + 1: 
                        # cluster width no larger than 24nt 
                        _p = _s
                        _n = locus_covs[_strand][_chr][_s]
                        for _e in cluster_dict[_strand][_chr][_s]:
                            _n += locus_covs[_strand][_chr][str(_e)]
                            if locus_covs[_strand][_chr][str(_e)] > locus_covs[_strand][_chr][_p]:
                                _p = str(_e)
                        #_id = _chr + '|' + _strand + '|NA|' + _p + '|' + _n
                        _id = '{}|{}|NA|{}|{}'.format(_chr,_strand,_p,_n)
                        _peak = int(_p)
                        _c_line = _sep.join(str(e) for e in [_chr,_start,cluster_dict[_strand][_chr][_s][-1],_id,_n,_strand,_peak - 1,_peak]) + '\n'
                        w.write(_c_line)
                    else: 
                        # cluster width larger than 24nt 
                       cov_dict = {}
                       cov_dict[_s] = locus_covs[_strand][_chr][_s]
                       for _e in cluster_dict[_strand][_chr][_s]:
                           cov_dict[str(_e)] = locus_covs[_strand][_chr][str(_e)]
                       while(len(cov_dict.keys()) > 0):
                           #print(cov_dict)
                           peaks = findPeak(cov_dict)
                           cov_dict = peaks[0]
                           #_id = _chr + '|' + _strand + '|NA|' + peaks[1][2] + '|' + peaks[1][3]
                           #print(peaks[0])
                           #print(peaks[1])
                           _id = '{}|{}|NA|{}|{}'.format(_chr,_strand,peaks[1][2],peaks[1][3])
                           _c_line = _sep.join(str(e) for e in [_chr,peaks[1][0],peaks[1][1],_id,peaks[1][3],_strand,peaks[1][2] - 1,peaks[1][2]]) + '\n'
                           w.write(_c_line)                 
        print('[INFO] {} reads are from multiple locus'.format(_tc2))
        end_time = datetime.datetime.now()
        print('[INFO] finished and the running time is {}'.format(end_time - start_time))
#"""
if len(sys.argv) < 4:
    print("Usage:python callCluster.py input_file input_format if_stranded out_cluster out_dis")
    sys.exit(1)
else:
    if_stranded = int(sys.argv[3])
    #input2dict(sys.argv[1],sys.argv[2],if_stranded)
    callCluser(sys.argv[1],sys.argv[2],if_stranded,sys.argv[4],sys.argv[5])
    
