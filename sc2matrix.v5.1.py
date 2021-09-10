# sc2matrix
# python=3.9
# Bin Zhang
# Adapted from Juexiao Zhou's script 
# Data: Sep 8, 2021
# version: 5.1
# ENV: 
# @ pip install cigar
# @ pip install portion
# @ pip install pysam
# @ pip install progressbar

import argparse
import os
import cigar
import pandas as pd
import portion as P
import datetime
from multiprocessing import Pool
import numpy as np
from collections import Counter
from collections import defaultdict
import pysam
from progressbar import ProgressBar, Bar, ETA


if __name__ == '__main__':
    
    #"""
    parser = argparse.ArgumentParser(description='sc2matrix')
    parser.add_argument('--bam' ,help='input alignment file in bam format')
    parser.add_argument('--ref', help='reference annotation in bed format')
    parser.add_argument('--barcode', help='barocde tsv file in the format from cellranger')
    parser.add_argument('--overlap_ratio', default='0', help='ratio of overlapped region among the min length of read and annotation')
    parser.add_argument('--overlap_nt', default='1', help='1 nucleotide')
    parser.add_argument('--t',default='1', help='numeber of thread')
    parser.add_argument('--s',default='1', help='strandless: 0, forward: 1, reverse: 2')
    parser.add_argument('--seed',default='3', help='length of seed for searching overlaps')
    parser.add_argument('--output',default='matrix.mtx', help='output file')
    parser.add_argument('--if_split',default='yes', help='if split the junction read into multiple blocks')
    parser.add_argument('--count_end',default='no', help='only count the 5 or 3 end of the read, incompatible with if_split yes')
    parser.add_argument('--model',default='normal', help='normal or noindex. normal model requires bam index')
    #parser.add_argument('--test_id',default='A01018:86:HG7GNDSXY:4:2533:24370:9251', help='test id')
    args = parser.parse_args()
        
    sam = args.bam #'test.scseq.1w.bam'
    ref = args.ref #'final.filtered.merged.cluster.cut2sam.refind.win48.bed'
    barcode = args.barcode #'barcodes.tsv'
    overlap_ratio_threshold = float(args.overlap_ratio) #0.01
    overlap_nt_threshold = int(args.overlap_nt)
    thread = int(args.t)
    consider_strand = int(args.s)
    seed_len = int(args.seed)
    output = args.output
    if_split = args.if_split
    count_end = args.count_end
    model = args.model
    
    print('[INFO] loaded count module')
    start_time = datetime.datetime.now()
    print('[INFO] parsing reference annotation file')

    ref_data = {}
    ref_region = {}
    _ref_id = 0
    with open(ref,'r') as f:
        for line in f:
            line = line.rstrip('\n').split('\t')
            _ref_id += 1
            _chr = line[0]
            _start = int(line[1])
            _end = int(line[2])
            _strand = line[5]
            _interval = P.closed(_start, _end)
            _start_seed = str(_start)[0:seed_len]
            _end_seed = str(_end)[0:seed_len]
            if ref_region.get(_chr,-1) == -1:
                ref_region[_chr] = [_start,_end]
            elif _start < ref_region[_chr][0]:
                ref_region[_chr][0] = _start
            elif _end > ref_region[_chr][1]:
                ref_region[_chr][1] = _end
            # determin strand 
            _strand_dict = {'+':'-','-':'+'}
            if consider_strand == 0:
                _strand = '?'
            elif consider_strand == 2:
                _strand = _strand_dict[_strand]
            if ref_data.get(_chr,-1) == -1:
                ref_data.setdefault(_chr,{'+' : {}, '-' : {}, '?' : {}})   
            # processed all possible seeds even though start and end have different seeds
            for _seed in range(int(_start_seed),int(_end_seed)+1):
                try:
                    ref_data[_chr][_strand]['{}'.format(_seed)].append({'start':_start,'end':_end,'id':_ref_id, 'interval':_interval})
                except:
                    ref_data[_chr][_strand]['{}'.format(_seed)]=[{'start':_start,'end':_end,'id':_ref_id, 'interval':_interval}]
    #print(ref_data['chr1']['-']['100'])
    print('[INFO] loading barcode')
    _cb_id = 0
    cb_dict = {}
    with open(barcode,'r') as f:
        for line in f:
            _cb_id += 1
            _cb = line.rstrip('\n')
            cb_dict[_cb] = _cb_id
            
    def len_interval(interval):
        #print(interval)
        out = 0
        for a in list(interval):
            out += a.upper - a.lower
        return out

    def len_union(interval1, interval2):
        return len_interval(interval1.union(interval2))

    def len_intersection(interval1, interval2):
        return len_interval(interval1.intersection(interval2))

    def len_difference(interval1, interval2):
        return len_interval(interval1.difference(interval2))

    def overlap_ratio(interval1, interval2):
        #return len_intersection(interval1,interval2)/len_union(interval1,interval2)
        return len_intersection(interval1,interval2)/min(len_interval(interval1),len_interval(interval2))
    #
    def cigar2interval(cigar_list):
        a = P.closed(cigar_list[0][0],cigar_list[0][1])
        if len(cigar_list)>1:
            for i in range(1,len(cigar_list)):
                a = a|P.closed(cigar_list[i][0],cigar_list[i][1])
        return a
            
    def find_overlap(_cig_interval, datas):
        out = []
        for data in datas:
            if _cig_interval.overlaps(data['interval']):
                #out.append([len_interval(_cig_interval),len_interval(data['interval']),overlap_ratio(_cig_interval,data['interval']),data,len_intersection(_cig_interval, data['interval'])])
                #out.append([len_interval(_cig_interval),len_interval(data['interval']),len_intersection(_cig_interval, data['interval']),data['id']])
                out.append([data['id'],overlap_ratio(_cig_interval,data['interval']),len_intersection(_cig_interval, data['interval'])])
                break
        if len(out)!=0:
            return out
        else:
            return False

    #_overlap = find_overlap(cigar2interval([[10024476,10024477]]),ref_data['chr1']['-']['100'])
    #print("testing overlap",_overlap[0])
    def get_barcode(tags):
        _CB = 'NA'
        _UB = 'NA'
        for eacht in tags:
            if(eacht[0] == "CB"):
                _CB = eacht[1]
            elif(eacht[0] == "UB"):
                _UB = eacht[1]
        return([_CB, _UB])
     
    # process a single read    
    def process_single_read(readstr):
        global chr_ref_data
        global seed_len
        global if_split
        global count_end
        #print(readstr)
        read = readstr.split('\t')
        _strand = read[0]
        _start = int(read[1]) - 1
        _cb = read[2]
        _UB = read[3]
        _cigar = read[4]
        seeds = []
        _seed = read[1][0:seed_len]
        cig_processed = []
        cig_list = list(cigar.Cigar(_cigar).items())
        ### if split junction read into multiple blocks 
        if if_split == 'yes':
            seeds.append(_seed)
            for ec in cig_list:
                if ec[1] == 'M':
                    _end = _start + ec[0]
                    cig_processed.append([_start,_end])
                    # if the junction jumped to the region with another seed 
                    if(str(_start)[0:seed_len] != _seed):
                        _seed = str(_start)[0:seed_len]
                        seeds.append(_seed)
                    # if the region start and end has different seed 
                    if(str(_end)[0:seed_len] != _seed):
                        _seed = str(_end)[0:seed_len]
                        seeds.append(_seed)
                    _start += ec[0]
                elif ec[1] == 'N': # junction indicate in the cigar 
                    _start += ec[0]
        else: 
            _end = _start
            for ec in cig_list:
                if ec[1] == 'M' or ec[1] == 'N':
                    _end = _end + ec[0]
            # only count 5'end or 3'end 
            if '5' in count_end:
                #print('[INFO] only counting the read 5p end')
                if _strand == '-':
                    cig_processed.append([_end,_end + 1])
                    seeds.append(str(_end)[0:seed_len])
                else:
                    cig_processed.append([_start,_start + 1])
                    seeds.append(_seed)
            elif '3' in count_end:
                #print('[INFO only counting the read 3p end]')
                if _strand == '+':
                    cig_processed.append([_end,_end + 1])
                    seeds.append(str(_end)[0:seed_len])
                else:
                    cig_processed.append([_start,_start + 1])
                    seeds.append(_seed)
            else:
                cig_processed.append([_start,_end])
                s1 = int(str(_start)[0:seed_len])
                s2 = int(str(_end)[0:seed_len]) + 1
                for _es in range(s1,s2):
                    seeds.append(str(_es))
        #
        cig_interval = cigar2interval(cig_processed)
        _outinfo = 'NA'
        for _seed in seeds:
            #if(_seed == '100'):
                #print(_UB,cig_interval,chr_ref_data[_strand][_seed])
            try:
                _overlap = find_overlap(cig_interval,chr_ref_data[_strand][_seed])
                #if [_UB == 'TATCTCTTCTGA']:
                    #print(cig_interval,chr_ref_data[_strand][_seed])
            except:
                _overlap = False
            if _overlap != False:
                _overlap = _overlap[0]
                _outinfo = '{}\t{}\t{}\t{}\t{}'.format(_overlap[0],_cb,_UB,_overlap[2],_overlap[1])
        return _outinfo

    #     
    def processed_data_to_MEX(datal):
        global overlap_ratio_threshold
        global overlap_nt_threshold
        count_dict = defaultdict(dict)
        mex_out = []
        _hit_read_num = 0
        _hit_umi_num = 0
        for _line in datal:
            if _line != 'NA':
                #print(_line)
                _hit_read_num += 1
                inf = _line.rstrip('\n').split('\t')
                _feature_id = inf[0]
                _barcode_id = inf[1]
                _umi = inf[2]
                _overlap_nt = int(inf[-2])
                _overlap_ratio = float(inf[-1])
                if _overlap_ratio > overlap_ratio_threshold and _overlap_nt >= overlap_nt_threshold:
                    try:
                        count_dict[_feature_id][_barcode_id][_umi] = _overlap_nt
                    except:
                        count_dict[_feature_id][_barcode_id] = {_umi : _overlap_nt}
        for _feature_id in count_dict:
            for _barcode_id in count_dict[_feature_id]:
                _count = len(count_dict[_feature_id][_barcode_id])
                _out = '{}\t{}\t{}\n'.format(_feature_id,_barcode_id,_count)
                _hit_umi_num += _count
                mex_out.append(_out)
        return [mex_out,[_hit_read_num,_hit_umi_num]]
    #     
    def write_mex(mex):
        mex_line = mex[0]
        with open(output,'a+') as w:
            if True:
                for _line in mex_line:
                    w.write(_line)
                
    print('[INFO] parsing bam file {}'.format(sam))
    _total_read_num = 0
    _not_skipped_num = 0
    processed_data = {}
#    samfile = pysam.AlignmentFile(sam, "rb", threads = thread) 
    if 'noindex' in model:
        print('[INFO] running no index model')
        all_reads = {}
        bam_start_time = datetime.datetime.now()
        for read in pysam.AlignmentFile(sam, "rb", threads = thread):
            _total_read_num += 1
            if (read.flag & 256) | (read.flag & 512) | (read.flag & 1024):
                continue
            barcodes = get_barcode(read.get_tags())
            # skip reads without cell barcode or UMI tags 
            if barcodes[0] == 'NA' or barcodes[1] == 'NA':
                continue
            elif cb_dict.get(barcodes[0],-1) == -1: # skip reads with cell barcode beyond the list from the provided barcode tsv file 
                continue
            else:
                _cb = cb_dict[barcodes[0]]
            _not_skipped_num += 1
            if read.flag & 16:
                _strand = '-'
            else:
                _strand = '+'        
            if consider_strand == 0:
                _strand = '?'
            _sep = '\t'
            _chr = read.reference_name
            all_reads.setdefault(_chr,[]).append(_sep.join(str(e) for e in [_strand,read.reference_start,_cb,barcodes[1],read.cigarstring]))
        
        bam_end_time = datetime.datetime.now()
        print('[INFO] processed {} reads and {} reads were skipped'.format(_not_skipped_num,_total_read_num - _not_skipped_num))
        print('[INFO] Time used in parsing bam file:{}'.format(bam_end_time - bam_start_time))
        print('[INFO] overlapping with reference annotation')
        for _chr in all_reads:
            print('[INFO] processing {}'.format(_chr))
            chr_start_time = datetime.datetime.now()
            # skip chrM
            if ref_data.get(_chr,-1) == -1:
                continue
            chr_ref_data = ref_data[_chr]
            chr_reads = all_reads[_chr]
            with Pool(thread) as p:
                widgets = [Bar(),ETA()]
                pbar = ProgressBar(widgets=widgets,maxval=len(chr_reads))
                chr_processed_data = list(pbar(p.map(process_single_read,chr_reads)))

            chr_end_time = datetime.datetime.now()
            processed_data[_chr] = chr_processed_data
            print("[INFO] Time used: {}".format(chr_end_time - chr_start_time))
    elif 'normal' in model:
        samfile = pysam.AlignmentFile(sam, "rb", threads = thread)
        for _chr in ref_region:
            _chr_start = ref_region[_chr][0]
            _chr_end = ref_region[_chr][1]
            chr_ref_data = ref_data[_chr]
            chr_start_time = datetime.datetime.now()
            #print(chr_ref_data['-']['100'])
            print('[INFO] processing alignment in {}:{}-{}'.format(_chr,_chr_start,_chr_end))
            inter = samfile.fetch(_chr,_chr_start,_chr_end)
            chr_reads = []
            for read in inter:
                _total_read_num += 1

                # skip secondary alignment and PCR duplicates 
                if (read.flag & 256) | (read.flag & 512) | (read.flag & 1024):
                    continue
                barcodes = get_barcode(read.get_tags())
                # skip reads without cell barcode or UMI tags 
                if barcodes[0] == 'NA' or barcodes[1] == 'NA':
                    continue
                elif cb_dict.get(barcodes[0],-1) == -1: # skip reads with cell barcode beyond the list from the provided barcode tsv file 
                    continue
                else:
                    _cb = cb_dict[barcodes[0]]
            
                _not_skipped_num += 1
                if read.flag & 16:
                    _strand = '-'
                else:
                    _strand = '+'        
                if consider_strand == 0:
                    _strand = '?'
                _sep = '\t'
                chr_reads.append(_sep.join(str(e) for e in [_strand,read.reference_start,_cb,barcodes[1],read.cigarstring]))

            print('[INFO] overlapping with reference annotation')
            with Pool(thread) as p:
                widgets = [Bar(),ETA()]
                pbar = ProgressBar(widgets=widgets,maxval=len(chr_reads))
                chr_processed_data = list(pbar(p.map(process_single_read,chr_reads)))

            chr_end_time = datetime.datetime.now()
            print('[INFO] processed {} reads'.format(_total_read_num))
            processed_data[_chr] = chr_processed_data
            print("[INFO] Time used: {}".format(chr_end_time - chr_start_time))
            #break
        print('[INFO] finally processed {} reads and {} reads were skipped'.format(_total_read_num,_total_read_num - _not_skipped_num))
    else:
        print('[INFO] Undetermined model option: {}'.format(model))
    
    
    print('[INFO] converting processed data to Market Exchange (MEX) format matrix')
    with Pool(thread) as p:
        mex_out = p.map(processed_data_to_MEX,processed_data.values())
    _line_num = 0
    _anno_read_num = 0
    _anno_umi_num = 0
    for _each_mex in mex_out:
        _line_num += len(_each_mex[0])
        _anno_read_num += _each_mex[1][0]
        _anno_umi_num += _each_mex[1][1]
    print('[INFO] finally {} reads and {} UMI were overlapped with reference annnotation'.format(_anno_read_num,_anno_umi_num))    
    print('[INFO] output to: {} in MEX Format with features: {} and cells: {}'.format(output,_ref_id,_cb_id))
    _head_line = '%%MatrixMarket matrix coordinate integer general\n%meta_data: { software_version: S3ESeq-beta}\n'
    with open(output,'w') as w:
        w.write(_head_line)
        w.write('{}\t{}\t{}\n'.format(_ref_id,_cb_id,_line_num))

    with Pool(thread) as p:
        p.map(write_mex, mex_out)    
    end_time = datetime.datetime.now()
    print("[INFO] Total Time used: {}".format(end_time - start_time))
    print('[INFO] finished')
    