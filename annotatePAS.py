# annotate PAS with gene annotation gtf
# python=3.6-3.10
# Bin Zhang
# Data: Nov 26, 2023
# version: 1.0
# ENV:

import os
import sys
import numpy as np
import pandas as pd
import pybedtools
from pybedtools import BedTool
#
from gtfparse import read_gtf
import portion as P

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
    
# split the genome into 7 distinct categories, including: terminal_exon, ups_exon, intron, ncRNA, ups_tss1k, extend_pas10k, intergenic 
def annotatePAS(input_pas,anno_gtf,output,header):
    print(f'### reading PAS from {input_pas}')
    pas_df = pd.read_csv(input_pas,sep = '\t',header = header)
    # default distance is 50bp
    pas_win50 = pas_df.copy()
    pas_win50.iloc[:,1] = pas_win50.iloc[:,2] - 25
    pas_win50.iloc[:,2] = pas_win50.iloc[:,2] + 25
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
            _pas = row['name']
            out_line = row.values.tolist()
            _type = 'NA'
            if 'name' in last_exon.columns and _pas in last_exon_pas['name'].values:
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
    print('### Done')
    
if len(sys.argv) < 4:
    print('python annotatePAS.py input_pas annotation_gtf output header')
    sys.exit(1)
else:
    if sys.argv[4] != 'none':
        header = 'infer'
    else:
        header = None
    annotatePAS(sys.argv[1],sys.argv[2],sys.argv[3],header)
    
        
    
        
    
    