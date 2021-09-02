# sc2matrix
# python=3.9
# Juexiao Zhou
# Data: Sep 1, 2021
# version: 3.3
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
import pysam
from progressbar import ProgressBar, Bar, ETA
    
if __name__ == '__main__':
    
    #"""
    parser = argparse.ArgumentParser(description='sc2matrix')
    parser.add_argument('--bam' ,help='path/to/sample.bam')
    parser.add_argument('--ref', help='path/to/sample_ref.bed')
    parser.add_argument('--barcode', help='path/to/sample_barcode.txt')
    parser.add_argument('--overlap_ratio', default='0.1', help='0.1')
    parser.add_argument('--t',default='1', help='numeber of thread')
    parser.add_argument('--s',default='1', help='ignore strand: 0, positive: 1, negative: 2')
    parser.add_argument('--head',default='2', help='length of seed for searching')
    args = parser.parse_args()

    sam = args.bam #'test.scseq.1w.bam'
    ref = args.ref #'final.filtered.merged.cluster.cut2sam.refind.win48.bed'
    barcode = args.barcode #'barcodes.tsv'
    overlap_ratio_threshold = float(args.overlap_ratio) #0.1
    thread = int(args.t)
    consider_strand = int(args.s)
    head_len = int(args.head)
    #"""

    """
    sam = 'test.1m.bam'
    ref = 'final.filtered.merged.cluster.cut2sam.refind.win48.bed'
    barcode = 'barcodes.tsv'
    overlap_ratio_threshold = 0.1
    thread = 112
    consider_strand = 1
    head_len = 3
    #"""

    print('[INFO] loaded count module')
    start_time = datetime.datetime.now()

    os.system('rm -r ./tmp')
    os.system('mkdir -p ./tmp')
    os.system('rm {}.processed'.format(sam))
    os.system('rm {}.matrix'.format(sam))
    print('[INFO] clear legacy data')
    print('[INFO] mkdir ./tmp')

    print('[INFO] parsing ref file')
    ref_data = {'+':{},'-':{},'?':{}}
    _refinfos = []
    with open(ref,'r') as f:
        for line in f:
            line = line.rstrip('\n').split('\t')
            _chr = line[0]
            _start = int(line[1])
            _end = int(line[2])
            _strand = line[5]
            _info = '{}:{}-{}'.format(_chr,_start,_end)
            _refinfos.append(_info)
            _interval = P.closed(_start, _end)
            _start_head = str(_start)[0:head_len]
            _end_head = str(_end)[0:head_len]
            for _head in range(int(_start_head),int(_end_head)+1):
                if consider_strand == 1:
                    try:
                        ref_data[_strand][_chr][_head].append({'start':_start,'end':_end,'info':_info, 'interval':_interval})
                    except:
                        try:
                            ref_data[_strand][_chr][_head]=[{'start':_start,'end':_end,'info':_info, 'interval':_interval}]
                        except:
                            ref_data[_strand][_chr]={'{}'.format(_head):[{'start':_start,'end':_end,'info':_info, 'interval':_interval}]}
                elif consider_strand == 0:
                    try:
                        ref_data['?'][_chr][_head].append({'start':_start,'end':_end,'info':_info, 'interval':_interval})
                    except:
                        try:
                            ref_data['?'][_chr][_head]=[{'start':_start,'end':_end,'info':_info, 'interval':_interval}]
                        except:
                            ref_data['?'][_chr]={'{}'.format(_head):[{'start':_start,'end':_end,'info':_info, 'interval':_interval}]}
                elif consider_strand == 2:
                    _strand_dict = {'+':'-','-':'+'}
                    try:
                        ref_data[_strand_dict[_strand]][_chr][_head].append({'start':_start,'end':_end,'info':_info, 'interval':_interval})
                    except:
                        try:
                            ref_data[_strand_dict[_strand]][_chr][_head]=[{'start':_start,'end':_end,'info':_info, 'interval':_interval}]
                        except:
                            ref_data[_strand_dict[_strand]][_chr]={'{}'.format(_head):[{'start':_start,'end':_end,'info':_info, 'interval':_interval}]}

    print('[INFO] parsed ref file')

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
        return len_intersection(interval1,interval2)/len_union(interval1,interval2)

    def cigar2interval(cigar_list):
        a = P.closed(cigar_list[0][0],cigar_list[0][1])
        if len(cigar_list)>1:
            for i in range(1,len(cigar_list)):
                a = a|P.closed(cigar_list[i][0],cigar_list[i][1])
        return a

    def process_cigar(pos, cig):

        # cig = '19M17725N72M'
        # process_cigar(3001708, cig)
        # OUT: [[3001708, 3001727], [3019452, 3019524]]

        cigar_list = list(cigar.Cigar(cig).items())
        out_list = []
        current_pos = pos
        for c in cigar_list:
            if c[1]=='M':
                out_list.append([current_pos,current_pos+c[0]])
                current_pos+=c[0]
            else:
                current_pos+=c[0]
        return out_list

    def find_overlap(_cig_interval, datas):
        out = []
        for data in datas:
            if _cig_interval.overlaps(data['interval']):
                out.append([len_interval(_cig_interval),len_interval(data['interval']),overlap_ratio(_cig_interval,data['interval']),data, len_intersection(_cig_interval, data['interval'])])
                #print(len_interval(_cig_interval))
                #print(len_interval(data['interval']))
                #print(overlap_ratio(_cig_interval,data['interval']))
                break
        if len(out)!=0:
            return out
        else:
            return False

    def process_single_line_from_sam(line):

        global ref_data
        global sam
        global consider_strand

        raw_line = line.rstrip('\n')
        line = line.rstrip('\n').split('\t')
        _chr = line[2]
        _start = int(line[3])
        _start_head = str(_start)[0:head_len]
        #_strand = '+' #TODO
        _strand = int(line[1])

        if _strand & 16:
            _strand = '-'
        else:
            _strand = '+'
            
        if consider_strand == 0:
            _strand = '?'
        
        _cig = line[5]
        _cig_processed = process_cigar(_start, _cig)
        _cig_interval = cigar2interval(_cig_processed)
        _CB = 'NA'
        _UB = 'NA'
        for x in line:
            if 'CB' in x:
                _CB = x
            elif 'UB' in x:
                _UB = x
        try:
            _overlap = find_overlap(_cig_interval,ref_data[_strand][_chr][_start_head])
        except:
            _overlap = False
        if _overlap!=False:
        #if True:
            _overlap = _overlap[0]
            _outinfo = '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(_CB, _UB,_overlap[3]['info'],_overlap[0],_overlap[1],_overlap[4], _overlap[2])
            #return line, _chr, _start, _cig, _CB, _UB, _cig_processed, _cig_interval, _overlap
            #print([line, _chr, _start, _cig, _CB, _UB, _cig_processed, _cig_interval, _overlap])
        else:
            _outinfo = '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(_CB, _UB,'NA','NA','NA','NA','NA')
        return _outinfo
        #with open('{}.processed'.format(sam),'a+') as w:
        #    w.write(_outinfo)
    
    print('[INFO] start parsing bam file')
    lines = []
    for line in pysam.AlignmentFile(sam,'rb'):
        if (line.flag & 256)|(line.flag & 512)|(line.flag & 1024):
            continue
        else:
            lines.append(line.tostring())

    #with open(sam,'r') as f:
    #    lines = f.readlines()
    with Pool(thread) as p:
        widgets = [Bar(),ETA()]
        pbar = ProgressBar(widgets=widgets,maxval=len(lines))
        processed_data = list(pbar(p.imap(process_single_line_from_sam,lines)))
        #processed_data = p.map(process_single_line_from_sam, lines)
        #p.map(process_single_line_from_sam, lines)

    #process_single_line_from_sam(_line)
    #process_single_line_from_sam(_line, ref_data)

    #print('[INFO] saving processed sam file to {}.processed'.format(sam))
    #with open('{}.processed'.format(sam),'w') as w:
    #    for line in processed_data:
    #        w.write(line)
    #print('[INFO] saved processed sam file to {}.processed'.format(sam))
    print('[INFO] finished porocessing bam file')
    
    print('[INFO] loading barcode')
    barcodes = []
    with open(barcode,'r') as f:
        for line in f:
            barcodes.append(line.rstrip('\n'))
    barcodes_dict = {}
    for barcode in barcodes:
        barcodes_dict[barcode]=0

    def sam2count(line):
        line = line.rstrip('\n').split('\t')
        _CB = line[0]
        _UB = line[1]
        _CB = _CB.split(':')[-1]
        _UB = _UB.split(':')[-1]
        _info = line[-5]
        if _info!='NA':
            _overlap_len = int(line[-2])
            _overlap_ratio = float(line[-1])
            if _overlap_ratio > overlap_ratio_threshold:
                with open('./tmp/{}.txt'.format(_info),'a+') as w:
                    w.write('{}_{}\n'.format(_CB,_UB))

    print('[INFO] loading processed sam')
    #sam_processed = open('{}.processed'.format(sam),'r')
    lines = []
    for line in processed_data:
        if 'CB' in line and 'UB' in line:
            lines.append(line)

    print('[INFO] start writing tmp files under ./tmp')
    with Pool(thread) as p:
        p.map(sam2count, lines)
        #widgets = [Bar(),ETA()]
        #pbar = ProgressBar(widgets=widgets,maxval=len(lines))
        #pbar(p.imap(sam2count,lines))

    def write_matrix_unit(_info):
        
        global barcodes_dict
        
        with open('{}.matrix'.format(sam),'a+') as w:
            if True:
            #    _info = 'chr1:3630763-3630811'
            #for _info in _refinfos:
                if os.path.isfile('./tmp/{}.txt'.format(_info)):
                    with open('./tmp/{}.txt'.format(_info),'r') as f:
                        lines = list(set(f.readlines()))
                        lines = [x.split('_')[0] for x in lines]
                        tmp_barcodes_dicts = barcodes_dict.copy()
                        for line in lines:
                            try:
                                tmp_barcodes_dicts[line]+=1
                                #print(line)
                            except:
                                pass
                        w.write('{}\t{}\n'.format(_info, '\t'.join([str(x) for x in tmp_barcodes_dicts.values()])))
                else:
                    w.write('{}\t{}\n'.format(_info, '\t'.join([str(x) for x in barcodes_dict.values()])))
    
    print('[INFO] start bulding matrix with size rows: {}, cols: {}'.format(len(_refinfos),len(barcodes)))
    with open('{}.matrix'.format(sam),'w') as w:
        w.write('{}\t{}\n'.format('region', '\t'.join(barcodes)))
        
    with Pool(thread) as p:
        #p.map(sam2count, lines)
        #widgets = [Bar(),ETA()]
        #pbar = ProgressBar(widgets=widgets,maxval=len(lines))
        #pbar(p.imap(write_matrix_unit,_refinfos))
        p.map(write_matrix_unit, _refinfos)
    
    print('[INFO] matrix saved to {}.matrix'.format(sam))

    #os.system('rm -rf ./tmp')
    print('[INFO] ./tmp deleted')

    end_time = datetime.datetime.now()
    print("[INFO] Time used: {}".format(end_time - start_time))
    print('[INFO] finished')    