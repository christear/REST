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
#import resource


if __name__ == '__main__':
    
    #"""
    parser = argparse.ArgumentParser(description='sc2matrix')
    parser.add_argument('--bam' ,help='input alignment file in bam format')
    parser.add_argument('--ref', help='reference annotation in bed format')
    parser.add_argument('--barcode', help='barocde tsv file in the format from cellranger')
    parser.add_argument('--overlap_ratio', default='0', help='ratio of overlapped region among the union of read and annotation')
    parser.add_argument('--overlap_nt', default='1', help='1 nucleotide')
    parser.add_argument('--t',default='1', help='numeber of thread')
    parser.add_argument('--s',default='1', help='ignore strand: 0, positive: 1, negative: 2')
    parser.add_argument('--seed',default='3', help='length of seed for searching overlaps')
    parser.add_argument('--output',default='out.matrix', help='output file')
    parser.add_argument('--if_split',default='yes', help='if split the junction read into multiple blocks')
    parser.add_argument('--count_end',default='no', help='only count the 5 or 3 end of the read, incompatible with if_split yes')
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
    tmpdir = '_tmp_' + output
#    test_id = args.test_id
    if_split = args.if_split
    count_end = args.count_end
    #"""
    #process = psutil.Process(os.getpid())
    """
    sam = 'test.1m.bam'
    ref = 'final.filtered.merged.cluster.cut2sam.refind.win48.bed'
    barcode = 'barcodes.tsv'
    overlap_ratio_threshold = 0.1
    thread = 112
    consider_strand = 1
    seed_len = 3
    #"""

    print('[INFO] loaded count module')
    start_time = datetime.datetime.now()

    os.system('if [ -d {} ]; then rm -r {}; fi'.format(tmpdir,tmpdir))
    os.system('mkdir -p {}'.format(tmpdir))
    os.system('if [ -f {} ]; then rm {}; fi'.format(output,output))
    print('[INFO] clear legacy data')
    print('[INFO] creating temporary dir:',tmpdir)
    print('[INFO] parsing reference annotation file')
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
            _start_head = str(_start)[0:seed_len]
            _end_head = str(_end)[0:seed_len]
            # processed all possible seeds even though start and end have different seeds
            for _head in range(int(_start_head),int(_end_head)+1):
                if consider_strand == 1:
                    try:
                        ref_data[_strand][_chr]['{}'.format(_head)].append({'start':_start,'end':_end,'info':_info, 'interval':_interval})
                    except:
                        try:
                            ref_data[_strand][_chr]['{}'.format(_head)]=[{'start':_start,'end':_end,'info':_info, 'interval':_interval}]
                        except:
                            ref_data[_strand][_chr]={'{}'.format(_head):[{'start':_start,'end':_end,'info':_info, 'interval':_interval}]}
                elif consider_strand == 0:
                    try:
                        ref_data['?'][_chr]['{}'.format(_head)].append({'start':_start,'end':_end,'info':_info, 'interval':_interval})
                    except:
                        try:
                            ref_data['?'][_chr]['{}'.format(_head)]=[{'start':_start,'end':_end,'info':_info, 'interval':_interval}]
                        except:
                            ref_data['?'][_chr]={'{}'.format(_head):[{'start':_start,'end':_end,'info':_info, 'interval':_interval}]}
                elif consider_strand == 2:
                    _strand_dict = {'+':'-','-':'+'}
                    try:
                        ref_data[_strand_dict[_strand]][_chr]['{}'.format(_head)].append({'start':_start,'end':_end,'info':_info, 'interval':_interval})
                    except:
                        try:
                            ref_data[_strand_dict[_strand]][_chr]['{}'.format(_head)]=[{'start':_start,'end':_end,'info':_info, 'interval':_interval}]
                        except:
                            ref_data[_strand_dict[_strand]][_chr]={'{}'.format(_head):[{'start':_start,'end':_end,'info':_info, 'interval':_interval}]}

    print('[INFO] parsed reference annotation file')

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
        
    # process a single read    
    def process_single_read(readstr):
        global ref_data
        global seed_len
        global if_split
        global count_end
#        print(readstr)
        read = readstr.split('\t')
        _chr = read[0]
        _strand = read[1]
        #_start = read[2]
        _start = int(read[2])
        _CB = read[3]
        _UB = read[4]
        _cigar = read[5]
        seeds = []
        _seed = read[2][0:seed_len]
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
        _outinfo = '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(_CB,_UB,'NA','NA','NA','NA','NA')
        for _seed in seeds:
            try:
                _overlap = find_overlap(cig_interval,ref_data[_strand][_chr][_seed])
            except:
                _overlap = False
            if _overlap != False:
                _overlap = _overlap[0]
                _outinfo = '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(_CB,_UB,_overlap[3]['info'],_overlap[0],_overlap[1],_overlap[4], _overlap[2])
        #with open('{}.processed'.format(sam),'a+') as w:
        #    w.write(_outinfo)
        return _outinfo
    # process a group of reads 
    def process_read_group(readl):
        outs = []
        for read in readl:
            _outinfo = process_single_read(read)
            outs.append(_outinfo)
        return outs
#   
    def get_barcode(tags):
        _CB = 'NA'
        _UB = 'NA'
        for eacht in tags:
            if(eacht[0] == "CB"):
                _CB = "CB:" + eacht[1]
            elif(eacht[0] == "UB"):
                _UB = "UB:" + eacht[1]
        return(_CB + '\t' + _UB)

    print('[INFO] start parsing bam file')
    #lines = []
    lines = {}
    readbam_start_time = datetime.datetime.now()
    _total_read_num = 0
#    _anno_match_num = 0
    for line in pysam.AlignmentFile(sam,'rb',threads = thread):
        _total_read_num += 1
        if (line.flag & 256)|(line.flag & 512)|(line.flag & 1024):
            continue
        else:
            _barcode = get_barcode(line.get_tags())
            # skip reads without cell barcode or UMI information 
            #if barcodes[0] == 'NA' or barcodes[1] == 'NA':
            if 'NA' in _barcode:
                continue        
            if line.flag & 16:
                _strand = '-'
            else:
                _strand = '+'        
            if consider_strand == 0:
                _strand = '?'
            
            _key = _total_read_num/500000
            _sep = '\t'
            #lines.setdefault(_key, []).append([line.reference_name,_strand,line.reference_start,barcodes,line.cigarstring])
            lines.setdefault(_key, []).append(_sep.join(str(e) for e in [line.reference_name,_strand,line.reference_start,_barcode,line.cigarstring]))
    readbam_end_time = datetime.datetime.now()
    print('[INFO]',_total_read_num,'reads have been porocessed from the bam file')
    print("[INFO] Time used for processing bam file: {}".format(readbam_end_time - readbam_start_time))

    print('[INFO] overlapping with reference annotation')
    anno_start_time = datetime.datetime.now()
    with Pool(thread) as p:
        widgets = [Bar(),ETA()]
        pbar = ProgressBar(widgets=widgets,maxval=len(lines))
        processed_data = list(pbar(p.map(process_read_group,lines.values())))
        print("[INFO] processed ",len(processed_data)," reads")
    anno_end_time = datetime.datetime.now()
    print('[INFO] finished overlapping with reference annotation ')
    print("[INFO] Time used for annotation: {}".format(anno_end_time - anno_start_time))    
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
        _n = 0
        if _info!='NA':
            _overlap_len = int(line[-2])
            _overlap_ratio = float(line[-1])
            if _overlap_ratio > overlap_ratio_threshold and _overlap_len >= overlap_nt_threshold:
                with open('{}/{}.txt'.format(tmpdir,_info),'a+') as w:
                    w.write('{}_{}\n'.format(_CB,_UB))
                    _n += 1
        return(_n)

    print('[INFO] loading processed results')
    lines = []
    _hits_num = 0
    for aline in processed_data:
        line = aline[0]
        if 'CB' in line and 'UB' in line:
            lines.append(line)
        if "NA" not in line:
            _hits_num += 1
    #
    print('[INFO]',_hits_num,' reads are overlapped with reference annotation')
    write_start_time = datetime.datetime.now()
    print('[INFO] start writing tmp files under {}'.format(tmpdir))
    with Pool(thread) as p:
        overlap_nums = p.map(sam2count, lines)
    _over_num = 0
    for n in overlap_nums:
        _over_num += n
    print('[INFO]',_over_num,' UMI are assigned to annotation based on threshhold',overlap_ratio_threshold,'and',overlap_nt_threshold)
    def write_matrix_unit(_info):
        
        global barcodes_dict
        
        with open(output,'a+') as w:
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
    
    print('[INFO] start bulding matrix with rows: {}, cols: {}'.format(len(_refinfos),len(barcodes)))
    with open(output,'w') as w:
        w.write('{}\t{}\n'.format('region', '\t'.join(barcodes)))
        
    with Pool(thread) as p:
        p.map(write_matrix_unit, _refinfos)
    
    print('[INFO] matrix saved to {}'.format(output))
    write_end_time = datetime.datetime.now()
    print("[INFO] Time used for writing: {}".format(write_end_time - write_start_time))
    os.system('rm -r {}'.format(tmpdir))  
    print('[INFO] {} deleted'.format(tmpdir))
    end_time = datetime.datetime.now()
    print("[INFO] Time used: {}".format(end_time - start_time))
#    print('[INFO] The maximum memory usage is {}'.format(memmax)) 
    print('[INFO] finished')