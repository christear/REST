# utils for PAS analysis 
# python=3.6-3.10
# Bin Zhang
# Data: Nov 12, 2023
# version: 1.0
# ENV: 

import re
import sys
import pandas as pd
#from seq_utils import kmer2seq
import pybedtools
from pybedtools import BedTool,example_filename
#
from gtfparse import read_gtf
import portion as P

# connect kmer to sequence: wrote by myself 
def kmer2seq2 (kmer_seq):
    kmer_list = re.split(' ',kmer_seq)
    _cs = ''
    for i in range(0,len(kmer_list) - 1):
        _cs += kmer_list[i][0]
    _cs += kmer_list[len(kmer_list) - 1]
    return _cs
# 
def toUpper (seq):
    us = seq.upper()
    return us
    
# copied from DNABERT/motif/motif_utils.py 
def seq2kmer(seq, k):
    """
    Convert original sequence to kmers
    
    Arguments:
    seq -- str, original sequence.
    k -- int, kmer of length k specified.
    
    Returns:
    kmers -- str, kmers separated by space

    """
    kmer = [seq[x:x+k] for x in range(len(seq)+1-k)]
    kmers = " ".join(kmer)
    return kmers
    
def kmer2seq(kmers):
    """
    Convert kmers to original sequence
    
    Arguments:
    kmers -- str, kmers separated by space.
    
    Returns:
    seq -- str, original sequence.

    """
    kmers_list = kmers.split(" ")
    bases = [kmer[0] for kmer in kmers_list[0:-1]]
    bases.append(kmers_list[-1])
    seq = "".join(bases)
    assert len(seq) == len(kmers_list) + len(kmers_list[0]) - 1
    return seq
    
# determine internal primming (IP)
# 7 consecutive A or more than 8 A among 10 nucleotide 
def determineIP (seq):
    _subseq = seq[100:110]
    #_subseq = seq
    if re.match(r'.*A{7}',_subseq):
        _ip = 1
    else:
        _na = 0
        for _n in _subseq:
            if _n == "A":
                _na += 1
        if _na > 7:
            _ip = 1
        else:
            _ip = 0
    return _ip

### determine PAS motif
# motif_list is an array of PAS motifs ranked by frequency, human and mouse may be different from the third  
# k: only check the top k motifs  
def determineMotif (seq,motif_list,k):
    # only using the upstream sequences of a PAS 
    _subseq = seq[0:104]
    _nm = 0
    for _i in range(0,k):
        if motif_list[_i] in _subseq:
            _nm += 1
            #_nm = 1
    return _nm
   
### relabel the training data based on PAS motif and internal primming 
def relabelPred (pred_txt,motif_file,output):
    pred_df = pd.read_csv(pred_txt,sep = '\t')
    cseq_list = pred_df["sequence"].apply(kmer2seq)
    ip_res = cseq_list.apply(determineIP)
    motif_df = pd.read_csv(motif_file,sep = '\t',header = None)
    motif_list = motif_df[0].values.tolist()
    # the top 3 PAS motif
    motif_res = cseq_list.apply(determineMotif,motif_list = motif_list,k = 3)
    # all possible PAS motif 
    motif_res2 = cseq_list.apply(determineMotif,motif_list = motif_list,k = len(motif_list))
    relabel = []
    #n = 0
    for i in range(0,len(motif_res)):
        # relabel PAS in FN without PAS motif and with internal primming as negative
        if pred_df["label"][i] == 1 and pred_df["pred_01"][i] == 0 and motif_res2[i] == 0 and ip_res[i] == 1:
            relabel.append(0)
        # relabel predicted PAS not from ground truth (FP)
        elif pred_df["label"][i] == 0 and pred_df["pred_01"][i] == 1:
            # relabel FP with PAS motif and without internal primming as positive
            if motif_res[i] > 0 and ip_res[i] == 0:
                relabel.append(1)
            # relabel FP without PAS motif or with internal primming as negative
            else:
                relabel.append(0)
        else:
            #n += 1
            relabel.append(pred_df["label"][i])
    out_df = pred_df
    out_df['relabel'] = relabel
    if output != None:
        #repred.to_csv(output, index=False, sep = "\t")
        out_df.to_csv(output, index = False, sep = '\t')
    #pred_df['relabel'] = relabel
    return out_df
    
### filter chr
def filter_chr(chr):
    r = 0
    if re.match(r"chr[1-9XY]",chr):
        r = 1
    return r
    
### get strand
def rev_strand (strand):
    str_dict = {"+":"-","-":"+"}
    return str_dict[strand]
    
### get peak
def get_peak (cluster_row):
    #values = cluster_row.values.tolist()[0]
    values = cluster_row.tolist()
    peak = values[7]
    if values[5] == "-":
        peak = values[6]
    return peak
    
### cluster to bed object 
def cluster2bed (cluster,flank_win,strand,if_filter_chr):
    cluster_df = pd.read_csv(cluster,sep = '\t',header = None)
    # only keep cluster from chr1-22 and chrX and chrY, remove chrM
    if if_filter_chr == True:
        fc = cluster_df[0].apply(filter_chr)
        cluster_df = cluster_df[fc == 1].copy()
    # adjust strand based on sequencing data 
    if strand == 2:
        cluster_df[5] = cluster_df[5].apply(rev_strand)
    # get peak from each cluster
    peaks = cluster_df.apply(get_peak,axis = 1)
    bed_df =  cluster_df.loc[:,[0,1,2,3,4,5]].copy()
    if flank_win > 1:
        bed_df[1] = peaks - flank_win
        bed_df[2] = peaks + flank_win
    elif flank_win == 1:
        bed_df[1] = peaks - flank_win
        bed_df[2] = peaks
    else:
        print(f'### Error: flank_win {flank_win} is not correct')
    # create bed from data.frame
    bed = BedTool.from_dataframe(bed_df)
    return bed
    
### add labels to bed object based on annotation (bedformat)
# entries has intersection with annotation or not will label with: true/1 and false/0
def label (x):
    if x:
        l = 1
    else:
        l = 0
    return l

def addLabel2bed(bed,annotation):
    bed_df = bed.to_dataframe()
    anno_bed = bed.intersect(annotation,s = True)
    anno_bed_df = anno_bed.to_dataframe()
    ids = bed_df['name']
    bed_df['anno'] = ids.isin(anno_bed_df['name'])
    bed_df['label'] = bed_df['anno'].apply(label)
    return bed_df
    
##### functions for annotatePAS
# extract last exons from each transcript for each gene and merge the overlapped ones 
def collapseExon(gene_id,exon_df):
    gene_df = exon_df[exon_df['gene_id'] == gene_id].copy()
    _chr = gene_df['seqname'].iloc[0]
    _strand = gene_df['strand'].iloc[0]
    _gene_name = gene_df['gene_name'].iloc[0]
    #txn_uid = gene_df['transcript_id'].unique()
    regs = P.empty()
    for i,row in gene_df.iterrows():
        regs = regs | P.closed(row['start'],row['end'])
    # format bed output
    j = 1
    out_bed = []
    if len(regs) == 1:
        each_r = regs[0]
        out_bed.append([_chr,each_r.lower,each_r.upper,':'.join(str(_e) for _e in [gene_id,_gene_name,'last_exon',len(regs)]),j,_strand])
    else:
        for each_r in regs:
            if j == 1:
                if _strand == '+':
                    out_bed.append([_chr,each_r.lower,each_r.upper,':'.join(str(_e) for _e in [gene_id,_gene_name,'first_exon',len(regs)]),j,_strand])
                else:
                    out_bed.append([_chr,each_r.lower,each_r.upper,':'.join(str(_e) for _e in [gene_id,_gene_name,'last_exon',len(regs)]),j,_strand])
            elif j == len(regs):
                if _strand == '+':
                    out_bed.append([_chr,each_r.lower,each_r.upper,':'.join(str(_e) for _e in [gene_id,_gene_name,'last_exon',len(regs)]),j,_strand])
                else:
                    out_bed.append([_chr,each_r.lower,each_r.upper,':'.join(str(_e) for _e in [gene_id,_gene_name,'first_exon',len(regs)]),j,_strand])
            else:
                out_bed.append([_chr,each_r.lower,each_r.upper,':'.join(str(_e) for _e in [gene_id,_gene_name,'mid_exon',len(regs)]),j,_strand])
            j += 1
    # 
    out_bed = pd.DataFrame(out_bed)
    return out_bed

# extend bed, the reason not using slopBed is due to it requires chr.length file 
def extendBedDf (bed_line,ext_len,direction):
    values = bed_line.tolist()
    #print(values)
    if direction == 'forward':
        if values[5] == '+':
            values[2] = values[2] + ext_len
        else:
            values[1] = max(values[1] - ext_len,1)
    elif direction == 'reverse':
        if values[5] == '+':
            values[1] = max(values[1] - ext_len,1)
        else:
            values[2] = values[2] + ext_len
    else:
        print(f'### direction {direction} is not determined')
    #values = pd.DataFrame(values)
    return values
    
# get neighboring region of a gene 
def neighborReg (bed_line,reg_len,direction):
    values = bed_line.tolist()
    #print(values)
    if direction == 'forward':
        if values[5] == '+':
            _start = values[2]
            _end = values[2] + reg_len
        else:
            _end = values[1]
            _start = max(values[1] - reg_len,1)
    elif direction == 'reverse':
        if values[5] == '+':
            _end = values[1]
            _start = max(values[1] - reg_len,1)
        else:
            _start = values[2]
            _end = values[2] + reg_len
    else:
        print(f'### direction {direction} is not determined')
    #values = pd.DataFrame(values)
    values[1] = _start
    values[2] = _end
    return values
    
##
# split the genome into 7 distinct categories, including: terminal_exon, ups_exon, intron, ncRNA, ups_tss1k, extend_pas10k, intergenic 
def annotatePAS(input_pas,anno_gtf,output,header,distance):
    print(f'### reading PAS from {input_pas}')
    pas_df = pd.read_csv(input_pas,sep = '\t',header = header)
    # default distance is 50bp
    pas_win50 = pas_df.copy()
    #distance = 24
    pas_win50.iloc[:,1] = pas_win50.iloc[:,2] - distance
    pas_win50.iloc[:,2] = pas_win50.iloc[:,2] + distance
    # pas_bed = BedTool.from_dataframe(pas_df)
    print(f'### loading annotaton from {anno_gtf}')
    gtf_df = read_gtf(anno_gtf)
    gene_df = gtf_df[gtf_df['feature'] == 'gene'].copy()
    gene_bed_df = gene_df[['seqname','start','end']].copy()
    gene_bed_df['name'] = gene_df['gene_id'] + ':' + gene_df['gene_name']
    gene_bed_df['gene_type'] = gene_df['gene_type']
    gene_bed_df['strand'] = gene_df['strand']
    
    # 10k regions downstream of a gene 
    pas_10k = gene_bed_df.apply(neighborReg,reg_len = 10000,direction = 'forward',axis = 1)
    pas_10k = pd.DataFrame(pas_10k.tolist())
    pas_10k.columns = gene_bed_df.columns
    # 1k regions upstream of a gene 
    tss_1k = gene_bed_df.apply(neighborReg,reg_len = 1000,direction = 'reverse',axis = 1)
    tss_1k = pd.DataFrame(tss_1k.tolist())
    tss_1k.columns = gene_bed_df.columns
    
    # genic regions
    pcg_df = gene_bed_df[gene_bed_df['gene_type'] == 'protein_coding'].copy()
    ncg_df = gene_bed_df[gene_bed_df['gene_type'] != 'protein_coding'].copy()
    # 
    pcg_exon_df = gtf_df[(gtf_df['feature'] == 'exon') & (gtf_df['gene_type'] == 'protein_coding')].copy()
    # annotated transcirpt end 
    pcg_txn_df = gtf_df[(gtf_df['feature'] == 'transcript') & (gtf_df['gene_type'] == 'protein_coding')].copy()
    pcg_txn_bed_df = pcg_txn_df[['seqname','start','end']].copy()
    pcg_txn_bed_df['name'] = pcg_txn_df['gene_id'] + ':' + pcg_txn_df['gene_name'] + ':' + pcg_txn_df['transcript_id']
    pcg_txn_bed_df['transcript_type'] = pcg_txn_df['transcript_type']
    pcg_txn_bed_df['strand'] = pcg_txn_df['strand']
    pcg_txn1 = pcg_txn_bed_df.apply(neighborReg,reg_len = 1, direction = 'forward',axis = 1)
    pcg_txn1_df = pd.DataFrame(pcg_txn1.tolist())
    # to speed up, split dataframe with chromosome 
    #chrs = pcg_exon_df['seqname'].unique()
    chrs = (pas_df.iloc[:,0]).unique()
    all_output = {}
    for _chr in chrs:
        print(f'### processing chromosome {_chr}')
        #chr_pas_bed = BedTool.from_dataframe(pas_df[pas_df.iloc[:,0] == _chr])
        chr_pas_bed = BedTool.from_dataframe(pas_win50[pas_win50.iloc[:,0] == _chr])
        # collapse exon for protein coding genes 
        #chr_pcg_exon_df = pcg_exon_df[pcg_exon_df['seqname'] == _chr].copy()
        chr_pcg_exon_df = pcg_exon_df[pcg_exon_df.iloc[:,0] == _chr].copy()
        chr_gids = chr_pcg_exon_df['gene_id'].unique()
        bed_list = list(map(lambda _gid: collapseExon(_gid,chr_pcg_exon_df),chr_gids))
        collapse_exon_df = pd.concat(bed_list)
        # PAS overlapped with exon of protein coding gene 
        cexon_bed = BedTool.from_dataframe(collapse_exon_df)
        pas_cexon = chr_pas_bed.intersect(cexon_bed,s = True,wa = True, wb = True)
        pas_cexon_df = pas_cexon.to_dataframe()
        txn1_bed = BedTool.from_dataframe(pcg_txn1_df[pcg_txn1_df.iloc[:,0] == _chr])
        pas_txn1 = chr_pas_bed.intersect(txn1_bed,s = True,wa = True, wb = True)
        pas_txn1_df = pas_txn1.to_dataframe()
        last_exon_pas = pas_cexon_df[pas_cexon_df['blockCount'].str.contains('last_exon')].copy()
        # PAS overlapped with protein_coding gene but not exons
        pcg_bed = BedTool.from_dataframe(pcg_df[pcg_df.iloc[:,0] == _chr])
        pas_pcgnoexon = chr_pas_bed.intersect(cexon_bed,v = True,s = True).intersect(pcg_bed,s = True,wa = True,wb = True)
        pas_pcgnoexon_df = pas_pcgnoexon.to_dataframe()
        # PAS overlapped with ncRNA 
        ncg_bed = BedTool.from_dataframe(ncg_df[ncg_df.iloc[:,0] == _chr])
        pas_ncg = chr_pas_bed.intersect(pcg_bed,s = True,v = True).intersect(ncg_bed,s = True,wa = True, wb = True)
        pas_ncg_df = pas_ncg.to_dataframe()
        # PAS at the 10k region downstream of a gene 
        pas10k_bed = BedTool.from_dataframe(pas_10k[pas_10k.iloc[:,0] == _chr])
        pas_dns10k = chr_pas_bed.intersect(pcg_bed,s = True,v = True).intersect(ncg_bed,s = True,v = True).intersect(pas10k_bed,s = True,wa = True,wb = True)
        pas_dns10k_df = pas_dns10k.to_dataframe()
        # PAS at the 1k region upstream of TSS antisense 
        tss1k_bed = BedTool.from_dataframe(tss_1k[tss_1k.iloc[:,0] == _chr])
        pas_tss1k = chr_pas_bed.intersect(pcg_bed,s = True,v = True).intersect(ncg_bed,s = True,v = True).intersect(pas10k_bed,s = True,v = True).intersect(tss1k_bed,wa = True,wb = True)
        pas_tss1k_df = pas_tss1k.to_dataframe()
        #
        #print('### Summarizing results')
        chr_out_put = []
        #chr_pas = chr_pas_bed.to_dataframe()
        chr_pas = pas_df[pas_df.iloc[:,0] == _chr].copy()
        for index,row in chr_pas.iterrows():
            #_pas = row['name']
            out_line = row.values.tolist()
            _pas = out_line[3]
            _type = 'NA'
            if 'name' in last_exon_pas.columns and _pas in last_exon_pas['name'].values:
                each = last_exon_pas[last_exon_pas['name'] == _pas].copy()
                _type = 'last_exon'
            elif 'name' in pas_txn1_df.columns and _pas in pas_txn1_df['name'].values:
                each = pas_txn1_df[pas_txn1_df['name'] == _pas].copy()
                _type = 'intron_anno'
            elif 'name' in pas_pcgnoexon_df.columns and _pas in pas_pcgnoexon_df['name'].values:
                each = pas_pcgnoexon_df[pas_pcgnoexon_df['name'] == _pas].copy()
                _type = 'intron_unanno'
            elif 'name' in pas_cexon_df.columns and _pas in pas_cexon_df['name'].values:
                each = pas_cexon_df[pas_cexon_df['name'] == _pas].copy()
                _type = 'ups_exon'
            elif 'name' in pas_ncg_df.columns and _pas in pas_ncg_df['name'].values:
                each = pas_ncg_df[pas_ncg_df['name'] == _pas].copy()
                _type = 'ncRNA'
            elif 'name' in pas_dns10k_df.columns and _pas in pas_dns10k_df['name'].values:
                each = pas_dns10k_df[pas_dns10k_df['name'] == _pas].copy()
                _type = 'extend_10k'
            elif 'name' in pas_tss1k_df.columns and _pas in pas_tss1k_df['name'].values:
                each = pas_tss1k_df[pas_tss1k_df['name'] == _pas].copy()
                _type = 'tss_1k'
            else:
                out_line += ['intergenic',-1,-1,'NA']
                chr_out_put.append(out_line)
            if _type != 'NA':
                for i in range(each.shape[0]):
                    out_line2 = out_line + [_type]
                    out_line2 += each.iloc[i,:][7:10].tolist()
                    chr_out_put.append(out_line2)
        all_output[_chr] = chr_out_put
    #
    print(f'### writing the output to {output}')
    with open(output,'w') as w:
        for _chr in all_output:
            for each_line in all_output[_chr]:
                w.write('\t'.join(str(_e) for _e in each_line) + '\n')
    print('### annotating PAS Done')

    

    
    
    
    
    
    